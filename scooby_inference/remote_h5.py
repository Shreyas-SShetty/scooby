"""Read parts of a remote HDF5 / .h5ad file over HTTP range requests, without downloading the whole file.

Used by fetch_cima_donors.py. Two layers:
  RemoteFile     a seekable file object for h5py: small reads (file structure, chunk index), cached
  read_ranges()  elements [a, b) of a chunked 1-D dataset: h5py looks up where each needed chunk sits; the chunk
                 bytes are fetched over several connections and decompressed here.

Lessons from the CIMA server (ftp.cngb.org), which shaped this code:
  * ~50 KB/s per connection, ~0.3-0.5 MB/s in total.
  * It blocks the client's address for a while (TLS connections fail) after too many requests; this happened
    twice on 2026-11-20. So: at most MAX_RPS requests per second over all connections, small reads for the file
    structure (MIN_READ), a large cache so index pieces are not read twice, and on refusal the code waits
    (up to 10 minutes at a time) and carries on instead of failing.
  * Every finished step is saved in cache_dir: downloaded chunks (<dataset>_<start>.bin) and located chunks
    (<dataset>_located.tsv: start, byte offset, size), so a stopped run never repeats work.
"""
import io
import os
import threading
import time
import zlib
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np
import requests

MIN_READ = 8 * 1024   # file-structure reads are a few KB; larger reads only waste bandwidth
MAX_RPS = float(os.environ.get("REMOTE_H5_MAX_RPS", 4))  # requests per second over all connections

_pace_lock, _last_request = threading.Lock(), [0.0]


def _pace():
    with _pace_lock:
        wait = _last_request[0] + 1.0 / MAX_RPS - time.time()
        if wait > 0:
            time.sleep(wait)
        _last_request[0] = time.time()


def patient_get(session, url, a, b):
    """Bytes [a, b) of url. Never gives up: when the server refuses, waits longer and longer (up to 10 min)."""
    fails = 0
    while True:
        _pace()
        try:
            r = session.get(url, headers={"Range": f"bytes={a}-{b - 1}"}, timeout=300)
            if r.status_code == 206 and len(r.content) == b - a:
                return r.content
            why = f"HTTP {r.status_code}, {len(r.content)} bytes"
        except requests.RequestException as e:
            why = type(e).__name__
        fails += 1
        wait = min(600, 5 * 2 ** min(fails, 7))
        if fails >= 3:
            print(f"    server refusing ({why}); attempt {fails}, waiting {wait} s", flush=True)
        time.sleep(wait)


class RemoteFile(io.RawIOBase):
    def __init__(self, url):
        self.url, self.pos, self.cache = url, 0, {}
        self.s = requests.Session()
        self.size = self._size()
        self.nreq = 0

    def _size(self):
        fails = 0
        while True:
            try:
                return int(self.s.head(self.url, timeout=60, allow_redirects=True).headers["Content-Length"])
            except (requests.RequestException, KeyError, ValueError) as e:
                fails += 1
                wait = min(600, 5 * 2 ** min(fails, 7))
                print(f"    server refusing ({type(e).__name__}); attempt {fails}, waiting {wait} s", flush=True)
                time.sleep(wait)

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def readinto(self, buf):
        n = min(len(buf), self.size - self.pos)
        if n <= 0:
            return 0
        a = self.pos
        blk = a // MIN_READ  # cache in aligned MIN_READ blocks, so lookups are O(1) and pieces are never re-read
        e = min(self.size, a + n)
        out = bytearray()
        for k in range(blk, (e - 1) // MIN_READ + 1):
            if k not in self.cache:
                self.cache[k] = patient_get(self.s, self.url, k * MIN_READ, min(self.size, (k + 1) * MIN_READ))
                self.nreq += 1
            out += self.cache[k]
        off = a - blk * MIN_READ
        buf[:n] = out[off:off + n]
        self.pos += n
        return n


def open_remote(url):
    f = io.BufferedReader(RemoteFile(url), buffer_size=MIN_READ)
    return f, h5py.File(f, "r")


def _filters(dset):
    pl = dset.id.get_create_plist()
    return [pl.get_filter(i)[0] for i in range(pl.get_nfilters())]


def _decode(raw, filters, itemsize):
    for flt in reversed(filters):
        if flt == h5py.h5z.FILTER_DEFLATE:
            raw = zlib.decompress(raw)
        elif flt == h5py.h5z.FILTER_SHUFFLE:
            a = np.frombuffer(raw, np.uint8)
            raw = a.reshape(itemsize, -1).T.tobytes()
        else:
            raise ValueError(f"unsupported HDF5 filter {flt}")
    return raw


def read_rows(dset, url, a, b, **kw):
    """Elements [a, b) of a chunked 1-D dataset."""
    return read_ranges(dset, url, [(a, b)], **kw)[0]


def read_ranges(dset, url, ranges, workers=8, cache_dir=None, label="", progress_every=120, load=True):
    """Elements [a, b) for each (a, b) in ranges, from a chunked 1-D dataset. Chunks are looked up one by one
    (h5py, in this thread) and each is downloaded as soon as it is located (worker threads), so lookups and
    downloads overlap. With cache_dir, located chunks and downloaded chunks are saved as they come.
    load=False (needs cache_dir): only make sure every chunk is saved on disk, keep nothing in memory, return None
    (for downloads of many GB; read them later with read_cached)."""
    c = dset.chunks[0]
    filters, itemsize = _filters(dset), dset.dtype.itemsize
    starts = sorted({s for a, b in ranges if b > a for s in range((a // c) * c, b, c)})
    if not starts:
        return [np.zeros(0, dset.dtype) for _ in ranges]
    name = dset.name.strip("/").replace("/", "_")
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
    path_of = (lambda s: os.path.join(cache_dir, f"{name}_{s}.bin")) if cache_dir else (lambda s: None)
    located_path = os.path.join(cache_dir, f"{name}_located.tsv") if cache_dir else None
    located = {}
    if located_path and os.path.exists(located_path):
        for line in open(located_path):
            s, o, z = map(int, line.split())
            located[s] = (o, z)
    saved = {s for s in starts if cache_dir and os.path.exists(path_of(s))}
    to_get = [s for s in starts if s not in saved]
    print(f"    {label}: {len(starts):,} chunks; {len(saved):,} already saved; of the other {len(to_get):,}, "
          f"{sum(s in located for s in to_get):,} already located", flush=True)

    t0, lock, last = time.time(), threading.Lock(), [time.time()]
    stats = {"fetched": 0, "bytes": 0, "located": 0}

    def fetch(s, off, size):
        raw = patient_get(requests.Session(), url, off, off + size)
        if cache_dir:
            with open(path_of(s) + ".tmp", "wb") as fh:
                fh.write(raw)
            os.replace(path_of(s) + ".tmp", path_of(s))
        arr = np.frombuffer(_decode(raw, filters, itemsize), dset.dtype) if load else None
        with lock:
            stats["fetched"] += 1
            stats["bytes"] += size
            if time.time() - last[0] > progress_every:
                last[0] = time.time()
                el = time.time() - t0
                print(f"    {label}: fetched {stats['fetched']:,}/{len(to_get):,} chunks "
                      f"({stats['bytes'] / 1e6:,.0f} MB, {stats['bytes'] / el / 1e6:.2f} MB/s), "
                      f"located {stats['located']:,} new, {el / 60:.0f} min", flush=True)
        return arr

    futures = {}
    with ThreadPoolExecutor(workers) as ex:
        log = open(located_path, "a") if located_path else None
        for s in to_get:
            if s not in located:
                i = dset.id.get_chunk_info_by_coord((s,))
                located[s] = (i.byte_offset, i.size)
                stats["located"] += 1
                if log:
                    log.write(f"{s}\t{i.byte_offset}\t{i.size}\n")
                    if stats["located"] % 50 == 0:
                        log.flush()
            futures[s] = ex.submit(fetch, s, *located[s])
        if log:
            log.close()
        out = {s: f.result() for s, f in futures.items()}
    print(f"    {label}: all {len(starts):,} chunks ready ({(time.time() - t0) / 60:.1f} min)", flush=True)
    if not load:
        return None
    for s in saved:
        out[s] = np.frombuffer(_decode(open(path_of(s), "rb").read(), filters, itemsize), dset.dtype)

    res = []
    for a, b in ranges:
        if b <= a:
            res.append(np.zeros(0, dset.dtype))
            continue
        first = (a // c) * c
        part = np.concatenate([out[s] for s in range(first, b, c)])
        res.append(part[a - first:b - first])
    return res


def dset_info(dset):
    """What read_cached needs to know about a dataset, so it can run without the remote file."""
    return {"name": dset.name.strip("/").replace("/", "_"), "chunk": dset.chunks[0], "filters": _filters(dset),
            "dtype": dset.dtype.str}


def read_cached(info, cache_dir, ranges, out_dtype):
    """Elements of all ranges, joined into one preallocated array of out_dtype, from chunks saved by read_ranges.
    No network. Memory: the output plus one chunk."""
    c, dt = info["chunk"], np.dtype(info["dtype"])
    res, pos = np.empty(sum(b - a for a, b in ranges), out_dtype), 0
    for a, b in ranges:
        s = (a // c) * c
        while s < b:
            with open(os.path.join(cache_dir, f"{info['name']}_{s}.bin"), "rb") as fh:
                arr = np.frombuffer(_decode(fh.read(), info["filters"], dt.itemsize), dt)
            lo, hi = max(a, s) - s, min(b, s + c) - s
            res[pos:pos + hi - lo] = arr[lo:hi]
            pos += hi - lo
            s += c
    assert pos == len(res)
    return res


def chunk_bytes_per_element(dset, start=0, n=20):
    """Average compressed bytes per element, from n consecutive chunks from element `start` (for download
    estimates). Consecutive chunks share index nodes, so this needs few remote reads."""
    c, N = dset.chunks[0], dset.shape[0]
    first = min(start // c, max(N // c - n, 0)) * c
    return float(np.mean([dset.id.get_chunk_info_by_coord((first + k * c,)).size for k in range(n)])) / c

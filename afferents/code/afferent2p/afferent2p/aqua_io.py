"""Readers for AQuA2 batch outputs (yu-lab-vt/AQuA2, aqua_cmd_batch.m + CFU_cmd_batch.m).

Per input movie <base>.tif AQuA2 writes, next to it:
  <base>_AQuA2.mat            res struct, MATLAB v7.3 (HDF5): opts, evt1, fts1
                              (event features), dffMat1, datOrg1 (the movie!) ...
  <base>_AQuA2_Ch1.csv        feature table: one ROW per feature, one COLUMN per
                              event (WriteRowNames=1, no header)
  <base>_AQuA2_Ch1_curves.xlsx  aligned dF/F curves per event
  <base>_AQuA2_risingMaps...  rising-time maps
and CFU_cmd_batch.m adds
  <base>_AQuA2_res_cfu.mat    cfuInfo1 (nCFU x 9 cell: id, event ids, weight map
                              HxW, occurrence 1xT, raw curve, dF/F curve, time
                              window, non-time window, freq stats), cfuRelation,
                              cfuGroupInfo, cfuOpts, datPro (mean image 0..1)

We read the feature table from the csv (fast, no movie), and only the small
pieces of the .mat files we need, through h5py (v7.3) or scipy (v7).
"""
from __future__ import annotations
import re
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------------- MATLAB helpers

def _h5_deref(f, obj):
    """Turn an h5py dataset/group into python: numeric arrays (MATLAB column-major
    -> transposed back), char arrays -> str, cells -> list, structs -> dict."""
    import h5py
    if isinstance(obj, h5py.Group):
        return {k: _h5_deref(f, obj[k]) for k in obj.keys()}
    cls = obj.attrs.get("MATLAB_class", b"")
    cls = cls.decode() if isinstance(cls, bytes) else str(cls)
    if obj.attrs.get("MATLAB_empty", 0):
        return np.zeros((0,)) if cls != "char" else ""
    data = obj[()]
    if cls == "char":
        arr = np.asarray(data).T
        return "".join(chr(int(c)) for c in arr.ravel(order="F") if c) if arr.size else ""
    if cls == "cell":
        refs = np.asarray(data).T           # MATLAB shape
        out = np.empty(refs.shape, dtype=object)
        for idx in np.ndindex(refs.shape):
            out[idx] = _h5_deref(f, f[refs[idx]])
        return out
    if cls in ("logical",):
        return np.asarray(data).T.astype(bool)
    if cls == "struct":
        return {k: _h5_deref(f, obj[k]) for k in obj.keys()}
    arr = np.asarray(data)
    if arr.dtype.kind == "O":                # object references (a struct ARRAY stores each field as nCFU x 1 refs — v41)
        refs = arr.T
        out = np.empty(refs.shape, dtype=object)
        for idx in np.ndindex(refs.shape):
            out[idx] = _h5_deref(f, f[refs[idx]])
        return out
    if arr.dtype.names:                      # complex etc. — not expected
        return arr
    return arr.T if arr.ndim > 1 else arr


def read_mat_fields(path: Path, fields: list[str]) -> dict:
    """Read selected top-level / nested fields ('res/fts1/curve/dffMax') from a
    .mat file of either format. Missing fields -> None."""
    path = Path(path)
    out = {}
    try:
        import h5py
        with h5py.File(path, "r") as f:
            for fld in fields:
                node = f
                ok = True
                for part in fld.split("/"):
                    if part in node:
                        node = node[part]
                    else:
                        ok = False
                        break
                out[fld] = _h5_deref(f, node) if ok else None
        return out
    except OSError:
        pass                                  # not HDF5 -> v7
    from scipy.io import loadmat
    # v7: only the top-level variables asked for are read (a CFU file's cfuOpts alone is instant)
    m = loadmat(str(path), struct_as_record=False, squeeze_me=True,
                variable_names=sorted({fld.split("/")[0] for fld in fields}))
    for fld in fields:
        node = m
        ok = True
        for part in fld.split("/"):
            if isinstance(node, dict) and part in node:
                node = node[part]
            elif hasattr(node, part):
                node = getattr(node, part)
            else:
                ok = False
                break
        out[fld] = node if ok else None
    return out


# ----------------------------------------------------------------------------- file discovery

AQUA_SUFFIX = "_AQuA2"


def relocate(path_str: str, mouse_dir: Path) -> Path:
    """A path recorded in a manifest, re-rooted at this mouse folder when the
    folder has moved (another drive, a copy on another machine, Windows paths
    read on Linux): the part after the mouse folder's name is kept."""
    p = Path(path_str)
    if p.exists():
        return p
    parts = [x for x in re.split(r"[\\/]+", str(path_str)) if x]
    name = Path(mouse_dir).name
    if name in parts:
        i = len(parts) - 1 - parts[::-1].index(name)
        cand = Path(mouse_dir).joinpath(*parts[i + 1:])
        if cand.exists():
            return cand
    return p


def aqua_files(tif_path: Path) -> dict:
    """All AQuA2-related files for one input tiff (existing or not)."""
    tif_path = Path(tif_path)
    base = tif_path.with_suffix("")
    b = str(base)
    return {
        "tif": tif_path,
        "json": Path(b + ".json"),
        "meanproj": Path(b + "_meanproj.tif"),
        "behavior": Path(b + "_behavior.csv"),
        "res": Path(b + AQUA_SUFFIX + ".mat"),
        "features": Path(b + AQUA_SUFFIX + "_Ch1.csv"),
        "curves": Path(b + AQUA_SUFFIX + "_Ch1_curves.xlsx"),
        "cfu": Path(b + AQUA_SUFFIX + "_res_cfu.mat"),
        "cfu_pre": Path(b + AQUA_SUFFIX + "_cfu_pre.mat"),     # v42: AQuA2's pairwise event overlaps + per-event summary
        "opt": Path(b + "_Opt.csv"),
    }


# unpack_bin names its output <run>_bin31f_1p98s.tif (older versions used a dot,
# 1.98s — MATLAB's fileparts() cut names there); both spellings are accepted
BIN_RE = re.compile(r"_bin(?P<frames>\d+)f_(?P<sec>\d+[.p]\d+)s$")


def parse_bin_tag(tif: Path):
    """(frames_per_bin, seconds_per_frame) from a unpack_bin file name, or None."""
    m = BIN_RE.search(Path(tif).stem)
    if not m:
        return None
    try:
        return int(m["frames"]), float(m["sec"].replace("p", "."))
    except ValueError:
        return None


def pick_binned_tif(folder: Path, sec: float | None = 2.0, tol: float = 0.25,
                    log=print) -> Path | None:
    """The binned tiff of one run folder to feed AQuA2.

    `sec` is the intended bin length; the file name carries what it actually is
    (31 frames at 15.62 Hz = 1.98 s, not 2.00), so a tiff counts as a match when
    its tag is within `tol` seconds of `sec`. With several matches the closest
    wins; with none, and only one binned tiff in the folder, that one is used.
    """
    tifs = find_binned_tifs(folder)
    if not tifs:
        return None
    if sec is None:
        return tifs[0] if len(tifs) == 1 else max(tifs, key=lambda p: p.stat().st_mtime)
    tagged = [(t, parse_bin_tag(t)) for t in tifs]
    near = [(abs(tag[1] - sec), t) for t, tag in tagged if tag and abs(tag[1] - sec) <= tol]
    if near:
        best = min(d for d, _ in near)
        cands = [t for d, t in near if abs(d - best) < 1e-9]
        if len(cands) > 1:
            # the same bin length twice: an older dotted-name copy (1.98s) next to the
            # current one (1p98s) — the current naming wins; else the newest file
            undotted = [t for t in cands if "." not in BIN_RE.search(t.stem)["sec"]]
            pick = undotted[0] if len(undotted) == 1 else max(cands, key=lambda p: p.stat().st_mtime)
            log(f"  note: {folder.name}: {len(cands)} movies with the same bin length "
                f"({', '.join(t.name for t in cands)}) — using {pick.name}; the other is an older copy "
                f"that can be deleted with its sidecar files")
            return pick
        return cands[0]                                      # closest bin length
    if len(tifs) == 1:
        tag = parse_bin_tag(tifs[0])
        got = f"{tag[1]:.2f} s" if tag else "unknown bin length"
        log(f"  note: {folder.name}: only {tifs[0].name} ({got}) — using it although "
            f"{sec:.2f} s was requested")
        return tifs[0]
    log(f"  note: {folder.name}: no binned tiff within {tol:.2f} s of {sec:.2f} s "
        f"(found: {', '.join(t.name for t in tifs)}) — skipped")
    return None


def find_binned_tifs(folder: Path, pattern: str = "*_bin*s.tif", skip_legacy_twins: bool = True) -> list[Path]:
    """The unpack_bin.py outputs in one run folder (channel 1 only). An older
    dotted-name copy (..._1.98s.tif) next to its current twin (..._1p98s.tif)
    is left out (skip_legacy_twins): it is the same recording, processed by an
    older version, and only confuses the tables."""
    out = []
    for p in sorted(Path(folder).glob(pattern)):
        n = p.name
        if n.endswith("_meanproj.tif") or re.search(r"_ch\d+\.tif$", n) or AQUA_SUFFIX in n:
            continue
        if BIN_RE.search(p.stem) is None:                 # the bin tag must end the name: not
            continue                                      # _check_frames.tif and other sidecars
        out.append(p)
    if skip_legacy_twins and len(out) > 1:
        undotted = {p.stem for p in out if "." not in BIN_RE.search(p.stem)["sec"]}
        out = [p for p in out if "." not in BIN_RE.search(p.stem)["sec"]
               or p.stem.replace(".", "p") not in undotted]
    return out


# ----------------------------------------------------------------------------- feature table

FEATURE_ALIASES = {
    "Index": "event",
    "Starting Frame": "t0_frame",
    "Basic - Area": "area_um2",
    "Basic - Perimeter (only for 2D video)": "perimeter_um",
    "Basic - Circularity": "circularity",
    "Curve - P Value on max Dff (-log10)": "neg_log10_p",
    "Curve - Max Df": "df_max",
    "Curve - Max Dff": "dff_max",
    "Curve - Duration of visualized event overlay": "duration_s",
    "Curve - Duration 50% to 50% based on averge dF/F": "width55_s",
    "Curve - Duration 10% to 10% based on averge dF/F": "width11_s",
    "Curve - Rising duration 10% to 90% based on averge dF/F": "rise19_s",
    "Curve - Decaying duration 90% to 10% based on averge dF/F": "fall91_s",
    "Curve - dat AUC": "dat_auc",
    "Curve - df AUC": "df_auc",
    "Curve - dff AUC": "dff_auc",
    "Curve - Decay tau": "decay_tau_s",
}


def read_feature_table(csv_path: Path):
    """AQuA2 feature csv (rows = features, columns = events) -> pandas DataFrame
    with one row per event and short column names (FEATURE_ALIASES; other
    features keep their AQuA2 name)."""
    import pandas as pd
    raw = pd.read_csv(csv_path, header=None, index_col=0)
    raw.index = [str(i).strip() for i in raw.index]
    df = raw.T.reset_index(drop=True)
    df.columns = [FEATURE_ALIASES.get(c, c) for c in df.columns]
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    # a result without events is a table of feature names only (every line
    # 'Name,'): the transpose is one all-NaN row -> no events
    df = df.dropna(how="all").reset_index(drop=True)
    if "event" in df.columns:
        df = df[np.isfinite(df["event"])].reset_index(drop=True)
        df["event"] = df["event"].astype(int)
    else:
        df.insert(0, "event", np.arange(1, len(df) + 1))
    return df


def read_res_summary(res_path: Path) -> dict:
    """The small parts of <base>_AQuA2.mat: opts (frameRate = s/frame,
    spatialRes = um/px, sz), and event features straight from fts1 (used when
    the csv is missing or to get frames/peak frames not in the csv)."""
    want = ["res/opts/frameRate", "res/opts/spatialRes", "res/opts/sz", "res/opts/regMaskGap",
            "res/fts1/curve/tBegin", "res/fts1/curve/tEnd", "res/fts1/curve/dffMaxFrame",
            "res/fts1/curve/dffMax", "res/fts1/curve/dfMax", "res/fts1/curve/duration",
            "res/fts1/curve/width55", "res/fts1/curve/rise19", "res/fts1/curve/fall91",
            "res/fts1/curve/dffAUC", "res/fts1/curve/dffMaxPval", "res/fts1/basic/area",
            "res/fts1/loc/xSpa"]
    d = read_mat_fields(res_path, want)
    out = {}
    for k, v in d.items():
        key = k.split("/")[-1]
        if v is None:
            out[key] = None
        elif isinstance(v, np.ndarray) and v.dtype != object:
            out[key] = np.asarray(v).ravel() if key not in ("sz",) else np.asarray(v).ravel()
        else:
            out[key] = v
    for k in ("frameRate", "spatialRes", "regMaskGap"):
        if out.get(k) is not None:
            out[k] = float(np.asarray(out[k]).ravel()[0])
    return out


# ----------------------------------------------------------------------------- the registered movie

class RegisteredMovie:
    """AQuA2's own registered movie, read frame by frame out of <base>_AQuA2.mat
    (res.datOrg1: the binned movie after AQuA2's rigid registration, its border
    crop (regMaskGap + the motion margin) and its bleach correction, stored as
    uint16 in the original intensity units). This is the movie AQuA2 detected
    the events on, in the coordinates of the CFU maps and the event boxes — no
    crop offset needed. Indexing works like a (T, H, W) array: movie[t],
    movie[a:b], movie[t, y0:y1, x0:x1], movie[:, y0:y1, x0:x1]. Needs h5py
    (the .mat is HDF5, -v7.3); nothing is loaded until a frame is asked for."""

    def __init__(self, mat_path: Path):
        import h5py
        self.path = Path(mat_path)
        self._f = h5py.File(self.path, "r")
        try:
            self._ds = self._f["res"]["datOrg1"]
        except KeyError:
            self._f.close()
            raise ValueError(f"{self.path.name}: no res.datOrg1 (not an AQuA2 result?)")
        sh = tuple(int(v) for v in self._ds.shape)       # MATLAB H x W x L x T -> h5py (T, L, W, H); H x W x T -> (T, W, H)
        if len(sh) == 4:
            self._four = True
            T, L, W, H = sh
            if L != 1:
                self._f.close()
                raise ValueError(f"{self.path.name}: datOrg1 has {L} planes — 2-D movies only")
        elif len(sh) == 3:
            self._four = False
            T, W, H = sh
        else:
            self._f.close()
            raise ValueError(f"{self.path.name}: datOrg1 has shape {sh}")
        self.shape = (T, H, W)
        self.ndim = 3
        self.dtype = self._ds.dtype
        # older AQuA2 versions stored the movie scaled 0..1 (single); bring it back to intensity units
        self._scale = None
        if np.issubdtype(self.dtype, np.floating):
            try:
                mx = float(np.asarray(self._f["res"]["opts"]["maxValueDat1"]).ravel()[0])
                mn = float(np.asarray(self._f["res"]["opts"]["minValueDat1"]).ravel()[0])
                if np.max(self._raw(0)) <= 1.0 + 1e-6:
                    self._scale = (mn, mx)
            except Exception:
                self._scale = None

    def _raw(self, t):
        return np.asarray(self._ds[t, 0] if self._four else self._ds[t])

    @staticmethod
    def _norm(key, n):
        if isinstance(key, (int, np.integer)):
            k = int(key)
            if k < 0:
                k += n
            if not 0 <= k < n:
                raise IndexError(f"frame {key} out of range ({n})")
            return k
        if isinstance(key, slice):
            return key
        raise TypeError("RegisteredMovie takes ints and slices only")

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        key = key + (slice(None),) * (3 - len(key))
        t, y, x = (self._norm(k, n) for k, n in zip(key, self.shape))
        # h5py order: (T, [L,] W, H)
        sel = (t, 0, x, y) if self._four else (t, x, y)
        arr = np.asarray(self._ds[sel])
        if isinstance(t, int):
            arr = arr.T                                   # (W, H) -> (H, W), or (w) stays a vector
        else:
            arr = np.swapaxes(arr, -1, -2) if arr.ndim >= 2 else arr
        if self._scale is not None:
            mn, mx = self._scale
            arr = arr.astype(np.float32) * (mx - mn) + mn
        return arr

    def __len__(self):
        return self.shape[0]

    def __array__(self, dtype=None, copy=None):
        arr = self[:]
        return arr.astype(dtype) if dtype is not None else arr

    def close(self):
        try:
            self._f.close()
        except Exception:
            pass


def open_registered(tif_path: Path, log=None) -> "RegisteredMovie | None":
    """The registered movie of a run (res.datOrg1 of its _AQuA2.mat), or None
    when the result is missing, not HDF5 or has no movie (older AQuA2 output)."""
    f = aqua_files(tif_path)
    res = f["res"]
    if not res.exists() or res.stat().st_size == 0:
        return None
    try:
        return RegisteredMovie(res)
    except Exception as e:
        if log:
            log(f"  note: {res.name}: registered movie not readable ({e}) — using the binned tiff")
        return None


def open_movie(tif_path: Path, prefer_registered: bool = True, log=None) -> tuple:
    """(movie, source): the movie of a run for display and traces. 'registered' =
    AQuA2's own registered/cropped movie (CFU-map coordinates, offset (0, 0));
    'raw' = the binned tiff as recorded (AQuA2's crop sits at an offset in it)."""
    if prefer_registered:
        m = open_registered(tif_path, log=log)
        if m is not None:
            return m, "registered"
    import tifffile as tf
    try:
        arr = tf.memmap(str(tif_path), mode="r")
    except Exception:
        arr = tf.imread(str(tif_path))
    if arr.ndim == 2:
        arr = arr[None]
    return arr, "raw"


# ----------------------------------------------------------------------------- CFU results

CFU_FIELDS = ("id", "eventIds", "weightMap", "occurrence", "meanCurve", "meanDff", "timeWindow", "nonTimeWindow",
              "frequencyStats")             # the named fields of AQuA2's cfuInfo struct (1 Sep 2026 on) = the old cell's columns 1..9


def _struct_dict(x):
    """A scipy mat_struct (or a dict) -> dict of its fields; anything else -> None."""
    if isinstance(x, dict):
        return x
    try:
        return {k: getattr(x, k) for k in x._fieldnames}
    except Exception:
        return None


def _cfu_records(info) -> list[dict]:
    """AQuA2's cfuInfo1 in any of its layouts -> one dict per CFU with the
    keys of CFU_FIELDS (missing ones None):
      * the cell array of AQuA2 copies before 1 Sep 2026 (nCFU x 9/13/14
        columns: id, event ids, weight map, occurrence, raw curve, dF/F curve,
        time window, non-time window, frequency stats, ...),
      * the struct array of the current AQuA2 (named fields; scipy gives an
        object array of mat_structs — a single struct when nCFU is 1 —, h5py a
        dict of per-field arrays)."""
    if info is None:
        return []
    if isinstance(info, dict):                               # h5py: a struct array = dict of nCFU-long object arrays,
        vals = {k: info.get(k) for k in CFU_FIELDS}          # a scalar struct = dict of the fields themselves
        ev = vals.get("eventIds")
        if ev is None:
            return []
        if isinstance(ev, np.ndarray) and ev.dtype == object:
            n = int(ev.size)
            return [{k: (vals[k].ravel()[i] if isinstance(vals[k], np.ndarray) and vals[k].dtype == object
                         and vals[k].size == n else None) for k in CFU_FIELDS} for i in range(n)]
        if np.asarray(ev).size == 0 and (vals.get("weightMap") is None or np.asarray(vals["weightMap"]).size == 0):
            return []                                        # MATLAB_empty struct: no CFU
        return [vals]
    sd = _struct_dict(info)
    if sd is not None and "eventIds" in sd:                   # one CFU, squeezed by scipy
        return [{k: sd.get(k) for k in CFU_FIELDS}]
    arr = np.asarray(info, dtype=object)
    if arr.size == 0:
        return []
    first = arr.ravel()[0]
    if _struct_dict(first) is not None and "eventIds" in _struct_dict(first):   # struct array -> object array of mat_structs
        return [{k: _struct_dict(x).get(k) for k in CFU_FIELDS} for x in arr.ravel()]
    if arr.ndim == 1:                                        # a single CFU of the old cell layout, squeezed
        arr = arr.reshape(1, -1)
    recs = []
    for i in range(arr.shape[0]):
        row = arr[i]
        recs.append({k: (row[j] if arr.shape[1] > j else None) for j, k in enumerate(CFU_FIELDS)})
    return recs


def read_cfu(cfu_path: Path) -> dict:
    """<base>_AQuA2_res_cfu.mat -> dict with
        maps        (nCFU, H, W) float32 weight maps (0..1)
        events      list of int arrays: AQuA2 event ids (1-based) per CFU
        occurrence  (nCFU, T) bool: frames at which a CFU event rises
        dff         (nCFU, T) float32 CFU dF/F curves
        datPro      (H, W) mean image normalised 0..1
        freq        list of dicts (AQuA2 calcFreqStats) or None
        rules       str: which CFU rules made the file (from cfuOpts; v41)
    Reads the cell layout of AQuA2 copies before 1 Sep 2026 and the struct
    layout of the current AQuA2 (v41)."""
    d = read_mat_fields(cfu_path, ["cfuInfo1", "datPro", "cfuOpts"])
    info = d["cfuInfo1"]
    if info is None:
        raise ValueError(f"no cfuInfo1 in {cfu_path}")
    datpro = None if d["datPro"] is None else np.asarray(d["datPro"], np.float32)
    rules = cfu_rules_text(d.get("cfuOpts"))
    recs = _cfu_records(info)
    if not recs:                              # CFU_cmd_batch found no unit (empty cell / struct)
        shape = datpro.shape if datpro is not None and datpro.ndim == 2 else (0, 0)
        return {"maps": np.zeros((0,) + tuple(shape), np.float32), "events": [],
                "occurrence": np.zeros((0, 0), bool), "dff": np.zeros((0, 0), np.float32),
                "datPro": datpro, "freq": [], "n": 0, "rules": rules,
                "note": "no CFUs in this result (no place had enough repeated events)"}
    n = len(recs)
    maps, events, occ, dff, freq = [], [], [], [], []
    for r in recs:
        m = np.asarray(r["weightMap"], np.float32)
        if m.ndim == 3:                       # H x W x L with one slice
            m = m[:, :, 0] if m.shape[2] == 1 else m.max(axis=2)
        if m.ndim == 1 and datpro is not None and datpro.ndim == 2 and m.size == datpro.size:
            m = m.reshape(datpro.shape, order="F")   # the "equal-length vector" form of the current AQuA2
        maps.append(m)
        ev = np.asarray(r["eventIds"]).ravel().astype(int)
        events.append(ev)
        occ.append(np.asarray(r["occurrence"]).ravel().astype(bool))
        dff.append(np.asarray(r["meanDff"], np.float32).ravel())
        fr = _struct_dict(r["frequencyStats"]) if r["frequencyStats"] is not None else None
        freq.append(fr)
    return {"maps": np.stack(maps) if maps else np.zeros((0, 0, 0), np.float32),
            "events": events,
            "occurrence": np.stack(occ) if occ else np.zeros((0, 0), bool),
            "dff": np.stack(dff) if dff else np.zeros((0, 0), np.float32),
            "datPro": datpro,
            "freq": freq, "n": n, "rules": rules}


def cfu_rules_of(cfu_path: Path) -> str:
    """The CFU rules a `_res_cfu.mat` was made with (its cfuOpts only — fast), '' if unreadable."""
    try:
        return cfu_rules_text(read_mat_fields(cfu_path, ["cfuOpts"]).get("cfuOpts"))
    except Exception:
        return ""


def cfu_rules_text(cfu_opts) -> str:
    """One line on the CFU rules a result was made with, from its cfuOpts
    (overlap threshold, min events, and the post-merge correlation when the
    AQuA2 copy knew it)."""
    o = _struct_dict(cfu_opts) if cfu_opts is not None else None
    if not o:
        return ""
    det = _struct_dict(o.get("cfuDetect")) if o.get("cfuDetect") is not None else None
    if not det:
        return ""

    def num(k):
        v = det.get(k)
        try:
            return float(np.asarray(v).ravel()[0])
        except Exception:
            return None
    ov, mn, corr = num("overlapThr1"), num("minNumEvt1"), num("minimumCurveCorrelation")
    parts = []
    if ov is not None:
        parts.append(f"overlap ≥ {ov:g}")
    if mn is not None:
        parts.append(f"min events {mn:g}")
    if corr is not None:
        parts.append(f"post-merge r ≥ {corr:g}")
    return ", ".join(parts)


# ----------------------------------------------------------------------------- CFU curation sidecar

def curation_path(tif_path: Path) -> Path:
    """Where manual CFU decisions for one run are stored."""
    return Path(str(Path(tif_path).with_suffix("")) + AQUA_SUFFIX + "_res_cfu_curated.json")


def _written_stamp(d: dict) -> float | None:
    """The 'written' time of a curation file as a timestamp (None when absent / unreadable)."""
    if not d.get("written"):
        return None
    try:
        from datetime import datetime
        return datetime.fromisoformat(d["written"]).timestamp()
    except Exception:
        return None


def curation_stale_reason(d: dict, tif_path: Path) -> str:
    """Why the decisions in a curation file no longer belong to the run's current
    results — '' when they still do.

    Decisions name CFU numbers of ONE CFU result. They are stale when the AQuA2
    result (_AQuA2.mat) or the CFU file (_res_cfu.mat) was made after them (a
    re-run of the batch, new CFU rules such as another overlap threshold —
    v43: the CFU batch rewrites only the CFU file, so the CFU file has to be
    checked, not just the AQuA2 result), or when the CFU rules recorded with
    the decisions differ from the rules in the CFU file."""
    f = aqua_files(tif_path)
    res, cfu = f["res"], f["cfu"]
    stamp = d.get("aqua_result_mtime")
    if stamp is None:
        stamp = _written_stamp(d)
    if res.exists() and stamp is not None and res.stat().st_mtime > float(stamp) + 60:
        return "the AQuA2 result (_AQuA2.mat) was redone after the decisions were made"
    cstamp = d.get("cfu_result_mtime")
    if cstamp is None:
        cstamp = _written_stamp(d)
    if cfu.exists() and cstamp is not None and cfu.stat().st_mtime > float(cstamp) + 60:
        return "the CFU result (_res_cfu.mat) was redone after the decisions were made (new CFU rules / threshold)"
    if d.get("cfu_rules") and cfu.exists():
        now = cfu_rules_of(cfu)
        if now and now != d["cfu_rules"]:
            return f"the decisions were made under '{d['cfu_rules']}', the CFU file now has '{now}'"
    return ""


def read_curation(tif_path: Path) -> dict | None:
    """{'rejected': [1-based CFU indices], 'added': [ {...} ], 'notes': str} or None.

    Written by cfu_curate.py; honoured by aqua_analyze.py. Kept as a sidecar so
    AQuA2's own result files are never modified and the curation survives a
    re-run of the CFU batch. Decisions that belong to an earlier result (see
    curation_stale_reason) come back emptied with 'stale': True and
    'stale_reason': why — they are never applied to the new CFU numbers.
    """
    p = curation_path(tif_path)
    if not p.exists():
        return None
    import json
    try:
        d = json.loads(p.read_text(encoding="utf-8"))
    except Exception as e:
        print(f"  warning: could not read {p.name} ({e}) — ignoring the curation")
        return None
    d.setdefault("rejected", [])
    d.setdefault("added", [])
    d["rejected"] = [int(i) for i in d["rejected"]]
    # decisions refer to the CFU numbers of one AQuA2 result: a result made after the
    # curation (new settings, re-run, other CFU rules) has other CFUs, so the old decisions must not apply
    why = curation_stale_reason(d, tif_path)
    if why:
        if d["rejected"] or d["added"]:
            print(f"  warning: {p.name}: {why} — the decisions are ignored (curate this run again; the old file is kept "
                  f"as *_superseded_<date>.json when you save)")
        d = dict(d, rejected=[], added=[], stale=True, stale_reason=why)
    return d


def write_curation(tif_path: Path, rejected, added=None, notes: str = "") -> Path:
    """Write the decisions of one run. A file already there whose decisions belong to an
    earlier result (curation_stale_reason) is kept beside it as
    <...>_res_cfu_curated_superseded_<written date>.json (only when it holds decisions),
    so a re-curation after a CFU re-run never destroys what was decided before."""
    import json
    from datetime import datetime
    p = curation_path(tif_path)
    f = aqua_files(tif_path)
    res, cfu = f["res"], f["cfu"]
    if p.exists():
        try:
            old = json.loads(p.read_text(encoding="utf-8"))
            old.setdefault("rejected", []); old.setdefault("added", [])
            if (old["rejected"] or old["added"] or old.get("notes")) and curation_stale_reason(old, tif_path):
                stamp = (old.get("written") or "unknown").replace(":", "").replace("T", "_")[:15]
                keep = p.with_name(p.stem + f"_superseded_{stamp}.json")
                k = 2
                while keep.exists():
                    keep = p.with_name(p.stem + f"_superseded_{stamp}_{k}.json"); k += 1
                p.replace(keep)
                print(f"  earlier decisions of another CFU result kept as {keep.name}")
        except Exception as e:
            print(f"  note: could not check the earlier curation file {p.name} ({e}); it is overwritten")
    p.write_text(json.dumps({
        "run": Path(tif_path).stem,
        "tif": str(tif_path),
        "rejected": sorted(int(i) for i in rejected),
        "added": added or [],
        "notes": notes,
        "written": datetime.now().isoformat(timespec="seconds"),
        "aqua_result_mtime": res.stat().st_mtime if res.exists() else None,     # the result these decisions belong to
        "cfu_result_mtime": cfu.stat().st_mtime if cfu.exists() else None,      # v43: and the CFU file (re-made by the CFU batch alone)
        "cfu_rules": cfu_rules_of(cfu) if cfu.exists() else "",                 # v43: the CFU rules the numbers refer to
    }, indent=2), encoding="utf-8")
    return p


# ----------------------------------------------------------------------------- GUI option files

OPT_TEMPLATE_ROWS = [
    # (Name, Variable, Type, default)  == cfg/parameters_for_batch.csv of AQuA2
    ("Registration mode", "registrateCorrect", "preprocessing", "1"),
    ("Bleach correction mode", "bleachCorrect", "preprocessing", "1"),
    ("Median filter radius (For salt and pepper noise)", "medSmo", "preprocessing", "0"),
    ("Spatial smoothing level", "smoXY", "preprocessing", "0.5"),
    ("", "", "", ""),
    ("Active voxels threshold scale", "thrARScl", "activeregion", "3"),
    ("Minimum duration of signal", "minDur", "activeregion", "5"),
    ("Minimum size", "minSize", "activeregion", "20"),
    ("Maximum size", "maxSize", "activeregion", "inf"),
    ("Circularity threshold", "circularityThr", "activeregion", "0"),
    ("Allowed distance in the same signal", "spaMergeDist", "activeregion", "0"),
    ("", "", "", ""),
    ("Whether need temporal segmentation or not", "needTemp", "tempSeg", "1"),
    ("Seed size ratio compared to region ", "seedSzRatio", "tempSeg", "0.01"),
    ("Significance Threshold", "sigThr", "tempSeg", "3.5"),
    ("Delay score for merging", "maxDelay", "tempSeg", "0.6"),
    ("Whether signals are too close so that further refine is needed", "needRefine", "tempSeg", "0"),
    ("Whether signals are not large enough so that further growing is needed", "needGrow", "tempSeg", "0"),
    ("", "", "", ""),
    ("Whether need spatial segmentation or not", "needSpa", "spaSeg", "1"),
    ("Source size relative to super event", "sourceSzRatio", "spaSeg", "0.01"),
    ("Sensitivity to detect source (Level 1 to 10)", "sourceSensitivity", "spaSeg", "8"),
    ("Whether to extend event temporally", "whetherExtend", "spaSeg", "1"),
    ("", "", "", ""),
    ("Whether detect global signals", "detectGlo", "glo", "0"),
    ("Minimum duration of global signal", "gloDur", "glo", "20"),
    ("", "", "", ""),
    ("Ignore decay tau calculation", "ignoreTau", "post", "1"),
    ("Propagation metric in different directions (propagation map is already obtained)", "propMetric", "post", "0"),
    ("Network features", "networkFeatures", "post", "0"),
    ("", "", "", ""),
    ("Propagation smoothness", "gtwSmo", "latent", "0.2"),
    ("Baseline segment length", "cut", "latent", "200"),
    ("Baseline window", "movAvgWin", "latent", "25"),
    ("Event show threshold on raw data", "minShow1", "latent", "0.2"),
    ("Propagation threshold minimum", "propthrmin", "latent", "0.5"),
    ("Propagation threshold step", "propthrstep", "latent", "0.1"),
    ("Propagation threshold maximum", "propthrmax", "latent", "0.5"),
    ("How long duration ratio should be seen as footprint", "compress", "latent", "0"),
    ("Context length for segmentation", "gapExt", "latent", "5"),
    ("Target temporal downsample length", "TPatch", "latent", "20"),
    ("", "", "", ""),
    ("Frame rate", "frameRate", "from data", "0.5"),
    ("Spatial resolution", "spatialRes", "from data", "0.5"),
    ("X cooridante for north vector", "northx", "from data", "0"),
    ("Y cooridante for north vector", "northy", "from data", "1"),
    ("Remove pixels on image boundary", "regMaskGap", "from data", "5"),
]


def write_text_free(path: Path, text: str, log=print, encoding: str = "utf-8", allow_alt: bool = True) -> Path:
    """Write `text` to `path`; when Windows refuses to overwrite the existing file
    (PermissionError: the file is owned by / locked by another program or user —
    Excel on a csv, a stuck reader, a file written under another account), first
    try to replace it (delete + write: works when the folder is ours and nothing
    holds the file open), else write it under the next free name <stem>_2<suffix>,
    _3, … (a valid MATLAB name, so run('…_2.m') still works) and say so. Returns the
    path actually written — callers put THAT name into the generated scripts, so a
    locked file never stops the batch or the plan (v40). allow_alt=False raises
    instead of using another name (files other tools look up by name)."""
    import os
    import stat
    path = Path(path)
    try:
        path.write_text(text, encoding=encoding)
        return path
    except PermissionError as e:
        if not path.exists():
            raise
        first = e
    try:                                               # replace: clear a read-only bit, delete, write again
        os.chmod(path, stat.S_IWRITE)
        path.unlink()
        path.write_text(text, encoding=encoding)
        log(f"  note: {path.name} could not be overwritten ({first}) — deleted and written again")
        return path
    except OSError:
        pass
    if not allow_alt:
        raise PermissionError(f"{path} is locked by another program or owned by another user ({first}) — close whatever "
                              f"holds the files of {path.parent} (Excel? a MATLAB or Cowork session of another user?), or take "
                              f"ownership of the folder, and run the step again") from first
    for k in range(2, 50):
        alt = path.with_name(f"{path.stem}_{k}{path.suffix}")
        try:
            alt.write_text(text, encoding=encoding)
        except PermissionError:
            continue
        log(f"  WARNING: {path.name} is locked by another program or owned by another user ({first}) — written as {alt.name} "
            f"instead; close the program that holds it (Excel?) or take ownership of the folder when you can")
        return alt
    raise first


def read_opt_csv(path: Path) -> dict:
    """A parameter file saved by the AQuA2 GUI (SaveOpts: Name,Variable,Type,
    Value,Notes) or a cfg/parameters.csv-style file -> {variable: value str}."""
    import pandas as pd
    t = pd.read_csv(path, dtype=str, keep_default_na=False)
    cols = list(t.columns)
    if "Variable" not in cols:
        # the other SaveOpts layout (newer GUI): ONE header row with the opts field
        # names (filePath1, ..., registrateCorrect, ..., frameRate, spatialRes, ...)
        # and one row of values
        known = {r[1] for r in OPT_TEMPLATE_ROWS if r[1]}
        if len(t) >= 1 and len(known & set(cols)) >= 10:
            row = t.iloc[0]
            return {c.strip(): str(row[c]).strip() for c in cols
                    if c.strip() and not c.startswith("Unnamed") and str(row[c]).strip() != ""}
        raise ValueError(f"{path.name}: neither a Name,Variable,Type,Value file nor a one-row opts table "
                         f"— not an AQuA2 parameter file")
    if "Value" in cols:
        vcol = "Value"
    elif "Default" in cols:
        vcol = "Default"
    else:
        vcol = cols[3]
    out = {}
    for _, r in t.iterrows():
        v = str(r["Variable"]).strip()
        if v:
            out[v] = str(r[vcol]).strip()
    return out

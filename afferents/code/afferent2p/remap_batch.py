#!/usr/bin/env python3
"""remap_batch — make an existing <mouse>\\aqua_batch\\ usable from another PC
where the data folder sits under a different drive or root (the server's D:
mapped as Z:, or the animal folder copied to a local disk), and point the
MATLAB scripts at that PC's AQuA2 folder.

    python remap_batch.py "Z:\\2photon\\Chao\\GCAMP_affrerents\\NavRearChrTIGRE-2608-1" --from "D:\\2photon" --to "Z:\\2photon"
    python remap_batch.py <mouse folder> --from D: --to Z: --aqua-dir "C:\\AQuA2"
    python remap_batch.py <mouse folder> --aqua-dir "C:\\AQuA2"          # only the AQuA2 folder

Rewrites batch.csv, manifest.json, run_aqua_batch.m and run_cfu_batch.m in
place (a .bak copy of each is kept the first time). The prefix match is
case-insensitive. Equivalent to rebuilding the batch on that PC with
aqua_batch.py, but keeps a batch that was already partly processed.
"""
from __future__ import annotations
import argparse
import re
import shutil
from pathlib import Path


def remap_text(txt: str, src: str | None, dst: str | None) -> tuple[str, int]:
    if not src or dst is None:
        return txt, 0
    src_n = src.rstrip("\\/")
    dst_n = dst.rstrip("\\/")
    # match the prefix with either separator, case-insensitively; keep what follows
    pat = re.compile(re.escape(src_n).replace(r"\\", r"[\\/]"), re.IGNORECASE)
    n = 0

    def rep(m):
        nonlocal n
        n += 1
        return dst_n

    out = pat.sub(rep, txt)
    return out, n


def remap_json(obj, src, dst, aqua_dir, counter):
    """Walk a manifest (dicts/lists/strings) and remap every path string."""
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            if k == "aqua_dir" and aqua_dir:
                out[k] = aqua_dir; counter["aqua"] += 1
            else:
                out[k] = remap_json(v, src, dst, aqua_dir, counter)
        return out
    if isinstance(obj, list):
        return [remap_json(v, src, dst, aqua_dir, counter) for v in obj]
    if isinstance(obj, str) and src:
        new, n = remap_text(obj, src, dst)
        counter["paths"] += n
        return new
    return obj


def set_aqua_dir(txt: str, aqua_dir: str) -> tuple[str, int]:
    pat = re.compile(r"cd\('([^']*)'\);")
    n = 0

    def rep(m):
        nonlocal n
        n += 1
        return "cd('" + aqua_dir.replace("'", "''") + "');"

    return pat.sub(rep, txt), n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mouse", type=Path, help="mouse folder (contains aqua_batch\\)")
    ap.add_argument("--from", dest="src", default=None, help=r"path prefix as written now, e.g. D:\2photon or D:")
    ap.add_argument("--to", dest="dst", default=None, help=r"what it is on this PC, e.g. Z:\2photon or Z:")
    ap.add_argument("--aqua-dir", default=None, help="AQuA2 folder on this PC (the cd(...) line of the MATLAB scripts)")
    a = ap.parse_args()
    if (a.src is None) != (a.dst is None):
        raise SystemExit("--from and --to go together")
    if a.src is None and a.aqua_dir is None:
        raise SystemExit("nothing to do: give --from/--to and/or --aqua-dir")
    bdir = a.mouse / "aqua_batch"
    if not bdir.is_dir():
        raise SystemExit(f"{bdir} not found — run aqua_batch.py first (on any PC)")
    for name in ("batch.csv", "manifest.json", "run_aqua_batch.m", "run_cfu_batch.m"):
        f = bdir / name
        if not f.exists():
            print(f"  {name}: not there, skipped"); continue
        txt = f.read_text(encoding="utf-8")
        m = 0
        if name == "manifest.json":
            import json
            cnt = {"paths": 0, "aqua": 0}
            data = remap_json(json.loads(txt), a.src, a.dst, a.aqua_dir, cnt)
            new = json.dumps(data, indent=2)
            n, m = cnt["paths"], cnt["aqua"]
        else:
            new, n = remap_text(txt, a.src, a.dst)
            if a.aqua_dir and name.endswith(".m"):
                new, m = set_aqua_dir(new, a.aqua_dir)
        if new != txt:
            bak = f.with_suffix(f.suffix + ".bak")
            if not bak.exists():
                shutil.copy2(f, bak)
            f.write_text(new, encoding="utf-8")
        print(f"  {name}: {n} path(s) remapped" + (f", AQuA2 folder set ({m})" if m else ""))
    print("done. In MATLAB on this PC:")
    print(f"  run('{bdir / 'run_aqua_batch.m'}')")
    print(f"  run('{bdir / 'run_cfu_batch.m'}')")


if __name__ == "__main__":
    main()

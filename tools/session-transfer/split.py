#!/usr/bin/env python3
"""split.py OUT_DIR [CHARS]: split OUT_DIR/condensed.md into OUT_DIR/chunks/chunk_NN.md at line boundaries; print JSON list of paths."""
import sys, os, json
out = sys.argv[1]; size = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
lines = open(os.path.join(out, 'condensed.md'), encoding='utf-8').read().split('\n')
os.makedirs(os.path.join(out, 'chunks'), exist_ok=True)
chunks, cur, n = [], [], 0
for l in lines:
    if n + len(l) > size and cur:
        chunks.append(cur); cur, n = [], 0
    cur.append(l); n += len(l) + 1
if cur: chunks.append(cur)
paths = []
for i, c in enumerate(chunks):
    p = os.path.abspath(os.path.join(out, 'chunks', f'chunk_{i:02d}.md'))
    open(p, 'w', encoding='utf-8').write('\n'.join(c)); paths.append(p)
print(json.dumps(paths))

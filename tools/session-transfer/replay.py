#!/usr/bin/env python3
"""replay.py OUT_DIR: replay Write/Edit events from file_writes.jsonl per path.
A file is 'complete' only if a full Write exists and every later Edit applied cleanly and no later Bash command
mentions the path (sed -i / python / cp could have changed it). Writes complete files to OUT_DIR/recovered/<path>
and a manifest OUT_DIR/recovered/MANIFEST.json."""
import sys, os, json, re
out = sys.argv[1]
ws = [json.loads(l) for l in open(os.path.join(out, 'file_writes.jsonl'), encoding='utf-8') if l.strip()]
# collect bash commands for staleness check
bash = []
for l in open(os.path.join(out, 'condensed.md'), encoding='utf-8'):
    if '] BASH (' in l: bash.append(l)
state = {}
for w in ws:
    inp = w['input']; tool = w['tool']; p = inp.get('file_path') or inp.get('path') or inp.get('local_path')
    if tool == 'Write' and p:
        state[p] = {'content': inp.get('content', ''), 'ok': True, 'edits': 0, 'last': w['ts'], 'written': w['ts']}
    elif tool in ('Edit', 'MultiEdit') and p:
        s = state.get(p)
        if s is None:
            state[p] = {'content': None, 'ok': False, 'edits': 1, 'last': w['ts'], 'written': None}; continue
        edits = inp.get('edits') or [{'old_string': inp.get('old_string', ''), 'new_string': inp.get('new_string', ''), 'replace_all': inp.get('replace_all', False)}]
        for e in edits:
            if s['content'] is not None and e['old_string'] in s['content']:
                s['content'] = s['content'].replace(e['old_string'], e['new_string']) if e.get('replace_all') else s['content'].replace(e['old_string'], e['new_string'], 1)
            else:
                s['ok'] = False
        s['edits'] += 1; s['last'] = w['ts']
man = {}
os.makedirs(os.path.join(out, 'recovered'), exist_ok=True)
for p, s in state.items():
    base = os.path.basename(p)
    later_bash = [b[:200] for b in bash if base in b and b[1:20].replace(' ', 'T') >= (s['last'] or '')[:19]]
    complete = s['ok'] and s['content'] is not None and not later_bash
    man[p] = {'complete': complete, 'edits': s['edits'], 'written': s['written'], 'last_change': s['last'],
              'chars': len(s['content'] or ''), 'later_bash_touches': later_bash[:5]}
    if s['content'] is not None:
        dst = os.path.join(out, 'recovered', p.lstrip('/').replace(':', ''))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        open(dst, 'w', encoding='utf-8').write(s['content'])
json.dump(man, open(os.path.join(out, 'recovered', 'MANIFEST.json'), 'w'), indent=1)
print(json.dumps({p: (m['complete'], m['chars']) for p, m in man.items()}, indent=0)[:3000])

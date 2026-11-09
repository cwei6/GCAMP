#!/usr/bin/env python3
"""Condense Claude Code Remote list_events page files into a readable transcript.

usage: condense.py OUT_DIR PAGE_FILE [PAGE_FILE ...]
Writes OUT_DIR/condensed.md  (chronological, deduped by event uuid)
       OUT_DIR/file_writes.jsonl (every Write/Edit/heredoc/project_write with full content)
       OUT_DIR/human_messages.md (genuine user messages only, verbatim)
       OUT_DIR/stats.json
"""
import json, sys, os, re

def load_page(path):
    s = open(path, encoding='utf-8', errors='replace').read()
    i = s.find('{"ccr"')
    if i < 0:
        i = s.find('{')
    j = s.rfind('</other-session')
    body = s[i:j] if j > i else s[i:]
    body = body.strip()
    try:
        return json.loads(body)['ccr']
    except Exception:
        # try trimming trailing junk
        k = body.rfind('}')
        return json.loads(body[:k+1])['ccr']

def msg_of(v):
    if not isinstance(v, dict):
        return None, {}
    iac = v.get('internal_anthropic_catchall', {}) if isinstance(v.get('internal_anthropic_catchall'), dict) else {}
    m = v.get('message') or iac.get('message')
    return m, (iac or v)

def text_of_result(c, limit):
    if isinstance(c, str):
        t = c
    elif isinstance(c, list):
        parts = []
        for p in c:
            if isinstance(p, dict):
                if p.get('type') == 'text':
                    parts.append(p.get('text', ''))
                elif p.get('type') == 'image' or 'source' in p:
                    parts.append('[image]')
                elif p.get('type') == 'tool_reference':
                    parts.append('[tool_reference %s]' % p.get('tool_name'))
                else:
                    parts.append(json.dumps(p)[:300])
            else:
                parts.append(str(p))
        t = '\n'.join(parts)
    else:
        t = json.dumps(c)
    if limit and len(t) > limit:
        t = t[:limit] + '\n…[truncated %d chars]' % (len(t) - limit)
    return t

FULL_RESULT_TOOLS = re.compile(r'memory|project|Projects', re.I)
FILE_TOOLS = {'Write', 'Edit', 'MultiEdit', 'NotebookEdit'}

def main():
    out = sys.argv[1]; pages = sys.argv[2:]
    os.makedirs(out, exist_ok=True)
    events = {}
    for p in pages:
        try:
            c = load_page(p)
        except Exception as e:
            print('WARN cannot parse', p, e, file=sys.stderr); continue
        for e in c.get('data', []):
            kind = [k for k in e if k != 'created_at']
            if not kind: continue
            k = kind[0]; v = e[k]
            uid = (v.get('uuid') if isinstance(v, dict) else None) or (json.dumps(e)[:200])
            events[uid] = (e.get('created_at', ''), k, v)
    evs = sorted(events.values(), key=lambda x: x[0])
    tool_names = {}
    lines = []; writes = []; human = []
    st = {'events': len(evs), 'first': evs[0][0] if evs else None, 'last': evs[-1][0] if evs else None,
          'human_messages': 0, 'tool_uses': 0, 'file_writes': 0, 'compactions': 0,
          'last_uuid': (evs[-1][2].get('uuid') if evs and isinstance(evs[-1][2], dict) else None)}
    for t, k, v in evs:
        m, meta = msg_of(v)
        ts = t[:19].replace('T', ' ')
        if k == 'system':
            sub = meta.get('subtype') if isinstance(meta, dict) else None
            if sub == 'compact_boundary':
                st['compactions'] += 1
                lines.append(f'\n---\n### [{ts}] === COMPACTION BOUNDARY ===\n')
            continue
        if k == 'result':
            r = meta if isinstance(meta, dict) else {}
            lines.append(f'[{ts}] (turn ended: {r.get("subtype","")}, {str(r.get("result",""))[:200]!r})')
            continue
        if not isinstance(m, dict):
            continue
        content = m.get('content')
        if isinstance(content, str):
            content = [{'type': 'text', 'text': content}]
        if not isinstance(content, list):
            continue
        synthetic = bool(meta.get('isSynthetic')) if isinstance(meta, dict) else False
        for part in content:
            if not isinstance(part, dict): continue
            tp = part.get('type')
            if tp == 'text':
                txt = part.get('text', '')
                if k == 'user':
                    if txt.startswith('This session is being continued from a previous conversation'):
                        lines.append(f'\n### [{ts}] COMPACTION SUMMARY (verbatim)\n\n{txt}\n')
                        open(os.path.join(out, 'compaction_summary.md'), 'w').write(f'<!-- compaction summary written {ts} UTC -->\n\n' + txt)
                    elif synthetic or txt.startswith('[SYSTEM NOTIFICATION') or txt.startswith('<system-reminder>') or txt.startswith('<task-notification'):
                        lines.append(f'[{ts}] (system/notification) {txt[:1500]}')
                    else:
                        st['human_messages'] += 1
                        human.append(f'## [{ts}]\n\n{txt}\n')
                        lines.append(f'\n### [{ts}] HUMAN USER:\n\n{txt}\n')
                else:
                    lines.append(f'[{ts}] ASSISTANT: {txt}')
            elif tp == 'tool_use':
                st['tool_uses'] += 1
                name = part.get('name'); inp = part.get('input') or {}
                tool_names[part.get('id')] = name
                if name in ('SendUserMessage',):
                    lines.append(f'\n[{ts}] ASSISTANT -> USER MESSAGE:\n{inp.get("message","")}\n')
                    for a in inp.get('attachments', []) or []:
                        lines.append(f'   attachment: {a}')
                elif name == 'SendUserFile':
                    lines.append(f'[{ts}] SendUserFile {json.dumps(inp)[:800]}')
                elif name in FILE_TOOLS:
                    st['file_writes'] += 1
                    writes.append({'ts': t, 'tool': name, 'input': inp})
                    if name == 'Write':
                        c = inp.get('content', '')
                        lines.append(f'[{ts}] WRITE {inp.get("file_path")} ({len(c)} chars) — full content in file_writes.jsonl\n' + (c if len(c) <= 4000 and inp.get('file_path','').endswith('.md') else c[:600]))
                    else:
                        lines.append(f'[{ts}] {name} {inp.get("file_path")}: old={str(inp.get("old_string",""))[:300]!r} new={str(inp.get("new_string",""))[:600]!r}')
                elif name == 'Bash':
                    cmd = inp.get('command', '')
                    if re.search(r"cat\s*>+\s*\S+\s*<<|tee\s+\S+\s*<<", cmd):
                        writes.append({'ts': t, 'tool': 'Bash-heredoc', 'input': inp}); st['file_writes'] += 1
                    lines.append(f'[{ts}] BASH ({inp.get("description","")}): {cmd[:1500]}')
                elif FULL_RESULT_TOOLS.search(name or '') or (name or '').startswith('mcp__remote-devices__device_commit'):
                    if 'write' in (name or '').lower() or 'commit' in (name or '').lower() or 'create' in (name or '').lower():
                        writes.append({'ts': t, 'tool': name, 'input': inp}); st['file_writes'] += 1
                    lines.append(f'[{ts}] TOOL {name}: {json.dumps(inp, ensure_ascii=False)[:6000]}')
                elif name in ('TaskCreate', 'TaskUpdate'):
                    lines.append(f'[{ts}] {name}: {json.dumps(inp, ensure_ascii=False)[:500]}')
                elif name in ('create_trigger', 'send_later') or 'trigger' in (name or '') or 'send_later' in (name or ''):
                    lines.append(f'[{ts}] TOOL {name}: {json.dumps(inp, ensure_ascii=False)[:3000]}')
                else:
                    lines.append(f'[{ts}] TOOL {name}: {json.dumps(inp, ensure_ascii=False)[:700]}')
            elif tp == 'tool_result':
                name = tool_names.get(part.get('tool_use_id'), '?')
                full = FULL_RESULT_TOOLS.search(name or '') is not None and 'write' not in (name or '').lower()
                limit = 0 if full else (2500 if name in ('Read',) else 900)
                r = text_of_result(part.get('content'), limit)
                if r.strip() and r.strip() not in ('[image]',):
                    lines.append(f'[{ts}]   -> result({name}): {r}')
            elif tp == 'thinking':
                continue
    open(os.path.join(out, 'condensed.md'), 'w').write('\n'.join(lines))
    with open(os.path.join(out, 'file_writes.jsonl'), 'w') as f:
        for w in writes: f.write(json.dumps(w, ensure_ascii=False) + '\n')
    open(os.path.join(out, 'human_messages.md'), 'w').write('\n'.join(human))
    json.dump(st, open(os.path.join(out, 'stats.json'), 'w'), indent=1)
    print(json.dumps(st))
    print('condensed chars', sum(len(l) + 1 for l in lines))

if __name__ == '__main__':
    main()

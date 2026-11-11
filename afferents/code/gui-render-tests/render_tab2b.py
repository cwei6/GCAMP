"""Render tab 2 at a window height; report every row of both pane frames (mapped or not)."""
import sys, os
sys.path.insert(0, "/home/claude/afferent2p")
import pipeline_gui as pg
from pathlib import Path
H = int(sys.argv[1]); shot = sys.argv[2] if len(sys.argv) > 2 else ""
pg.SETTINGS = Path("/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad/tab2/settings.json")
app = pg.App()
app.geometry(f"1180x{H}+0+0")
app.nb.select(1)
app.update_idletasks(); app.update()
tab = app.nb.nametowidget(app.nb.tabs()[1])
print("window", app.winfo_width(), "x", app.winfo_height(), " tab frame h", tab.winfo_height(), " req", tab.winfo_reqheight())
def first_text(w):
    for k in w.winfo_children():
        try:
            t = str(k.cget("text"))
            if t: return t[:38]
        except Exception: pass
    return type(w).__name__
def walk(w, depth=1):
    for ch in w.winfo_children():
        print("  " * depth + f"{type(ch).__name__:12s} y={ch.winfo_y():4d} h={ch.winfo_height():4d} req={ch.winfo_reqheight():4d} mapped={ch.winfo_ismapped()}  {first_text(ch)}")
        if isinstance(ch, (pg.ttk.Panedwindow, pg.ttk.Frame)) and depth < 3 and type(ch).__name__ in ("Panedwindow", "Frame"):
            walk(ch, depth + 1)
walk(tab)
if shot:
    os.system(f"import -window root {shot}")
app.destroy()

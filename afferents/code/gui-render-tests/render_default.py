"""Start the GUI with its own default geometry on a 1920x1080 and a 1920x1200 screen; screenshot tab 2."""
import sys, os
sys.path.insert(0, "/home/claude/afferent2p")
import pipeline_gui as pg
from pathlib import Path
shot = sys.argv[1]
pg.SETTINGS = Path("/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad/tab2/settings.json")
app = pg.App()
app.nb.select(1)
app.update_idletasks(); app.update()
tab = app.nb.nametowidget(app.nb.tabs()[1])
print("screen", app.winfo_screenwidth(), app.winfo_screenheight(), "window", app.winfo_width(), "x", app.winfo_height(), " tab h", tab.winfo_height())
def walk(w, depth=1):
    for ch in w.winfo_children():
        n = type(ch).__name__
        if n in ("Treeview", "Panedwindow") or (n == "Frame" and depth <= 3):
            print("  " * depth + f"{n:12s} h={ch.winfo_height():4d} req={ch.winfo_reqheight():4d} mapped={ch.winfo_ismapped()}")
        if n in ("Panedwindow", "Frame") and depth < 3:
            walk(ch, depth + 1)
walk(tab)
os.system(f"import -window root {shot}")
app.destroy()

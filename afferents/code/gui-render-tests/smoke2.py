import sys, os
sys.path.insert(0, "/home/claude/afferent2p")
import pipeline_gui as pg
from pathlib import Path
SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
pg.SETTINGS = Path(SP + "/tab2/settings.json")
app = pg.App()
app.v["mouse"].set(SP + "/bed45/NavTIGRETreg-2609-1")
app.nb.select(1); app.refresh_status(); app.update()
app.tree.selection_set("0"); app.fill_stretches(); app.update()
app.open_sheet(); app.update()
sheet = [w for w in app.winfo_children() if isinstance(w, pg.SheetWindow)][0]
print("sheet", sheet.winfo_width(), sheet.winfo_height())
# press the zoom/pan button
btn = [w for w in sheet.winfo_children()[0].winfo_children() if isinstance(w, pg.ttk.Button) and str(w.cget("text")).startswith("Zoom")][0]
print("button:", btn.cget("text"))
btn.invoke(); app.update()
viewers = [w for w in sheet.winfo_children() if isinstance(w, pg.ImageViewer)]
print("viewer", [v.title() for v in viewers], viewers[0].files if viewers else None)
import time
for _ in range(10):
    app.update(); time.sleep(0.05)
os.system(f"import -window root {SP}/tab2/smoke_viewer.png")
# the other tabs render too
for i in range(len(app.nb.tabs())):
    app.nb.select(i); app.update()
print("tabs ok")
app.destroy()

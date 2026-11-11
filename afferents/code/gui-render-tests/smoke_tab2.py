"""GUI smoke test: mouse folder = bed45, tab 2, select the run, fill the stretches, open the check sheet, screenshot."""
import sys, os
sys.path.insert(0, "/home/claude/afferent2p")
import pipeline_gui as pg
from pathlib import Path
SP = "/tmp/claude-0/-home-claude/851dc502-c097-522d-aca4-96dcf64b8a98/scratchpad"
pg.SETTINGS = Path(SP + "/tab2/settings.json")
app = pg.App()
app.v["mouse"].set(SP + "/bed45/NavTIGRETreg-2609-1")
app.nb.select(1)
app.refresh_status()
app.update_idletasks(); app.update()
print("rows", len(app.status_rows))
app.tree.selection_set("0"); app.fill_stretches()
app.update()
print("stretch rows", len(app.stree.get_children()), [app.stree.item(i)["values"] for i in app.stree.get_children()])
app.open_sheet()
app.update_idletasks(); app.update()
tops = [w for w in app.winfo_children() if isinstance(w, pg.tk.Toplevel)]
print("toplevels", [w.title() for w in tops])
if tops:
    tops[0].geometry("1200x800+600+60")
app.update()
os.system(f"import -window root {SP}/tab2/smoke_sheet.png")
# decisions: exclude-all on the first stretch, then back to default
app.stree.selection_set(app.stree.get_children()[0])
app.decide("exclude-all"); app.update()
print("after exclude-all:", app.stree.item(app.stree.get_children()[0])["values"][-2:])
app.stree.selection_set(app.stree.get_children()[0])
app.decide(None); app.update()
print("after default:", app.stree.item(app.stree.get_children()[0])["values"][-2:])
app.destroy()

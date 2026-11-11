# Running the pipeline on another PC (data on a mapped drive)

The data stays on the server (`D:\2photon\Chao\GCAMP_affrerents\...` there);
on the other PC that folder is reached as `Z:\2photon\Chao\GCAMP_affrerents\...`.
Nothing in the data folders needs to change: every tool finds a run's files
by name next to the movie, and the only files that hold absolute paths (the
`aqua_batch\` folder of an animal) are rewritten once with `remap_batch.py`.

## 1. Once on the new PC

1. Copy this folder (`afferent2p`) to the PC, e.g. `C:\afferent2p` — from
   `Z:\2photon\Chao\GCAMP_affrerents\afferent2p_v29.zip`, unzipped. Do not
   run it from Z: (Python and MATLAB scripts run slowly from a share, and the
   pipeline writes its settings file next to itself).
2. Python: Anaconda (or any Python 3.10+). In an Anaconda Prompt:
   `cd C:\afferent2p` then `setup_remote.bat` (runs
   `python -m pip install -r requirements.txt`). `sbxreader` is only needed
   for Step 1 (unpacking `.sbx`); Steps 3–6 do not need it.
3. MATLAB with AQuA2: copy the server's `H:\MATLAB\AQuA2-main` to the PC
   (e.g. `C:\AQuA2`) or note where it already is. AQuA2 needs the Image
   Processing, Statistics, Signal Processing and Curve Fitting toolboxes. Run
   `patch_aqua2_toolbar('C:\AQuA2')` from this folder once (after every
   AQuA2 update) if `aqua_gui` fails with `Property assignment is not
   allowed when the object is empty`; the batch does not need the GUI.
4. Start the window: `run_gui.bat` (or `python pipeline_gui.py`). Set
   *Mouse folder* to `Z:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1`,
   *AQuA2 folder* to the AQuA2 folder on this PC and *ImageJ / Fiji* to its
   executable on this PC (tab 1b opens the check stacks with it); the window
   remembers them (in `pipeline_gui_settings.json` next to the scripts).

## 2. Per animal: make the batch files point at Z:

The `aqua_batch\` folder inside the animal folder was written on the server
with `D:\...` paths and the server's AQuA2 folder. Either rebuild it on this
PC (tab 3 in the window, or the command below — it writes Z: paths because
that is where the movies are seen from here):

```
python aqua_batch.py "Z:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --aqua-dir "C:\AQuA2" --um-per-px 0.792 --border 5
```

or, to keep a batch that is already partly processed, rewrite the paths in place:

```
python remap_batch.py "Z:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1" --from "D:\2photon" --to "Z:\2photon" --aqua-dir "C:\AQuA2"
```

(`.bak` copies are kept; run it again with `--from Z: --to D:` to go back.)
Then in MATLAB on this PC:

```
run('Z:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_batch\run_aqua_batch.m')
run('Z:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_batch\run_cfu_batch.m')
```

AQuA2 reads each 387-MB movie from the share and writes the ~1-GB result
back to it, so a batch is slower over the network than on the server; if the
network is slow, copy the animal folder to a local disk, run the batch there
with local paths (rebuild the batch on the local copy), and copy the run
folders' new files (`*_AQuA2*`) back.

Everything else is as in `GUIDE.md`, with `Z:\...` in place of `D:\...`:
Step 1b (the check stacks are written next to the movies on the share — once
per recording with `aqua_clean.py --force`, ~15 s each over the network —
and the decisions in tab 1b), Step 5 curation, Step 6 analysis — their
outputs go next to the movies / into `<mouse>\aqua_analysis\` on the share,
where the server sees them too. Results made on either PC are interchangeable:
the manifest is re-rooted automatically when a folder is read from a
different drive, and the run folders never hold absolute paths that matter.

## What must NOT be done from two PCs at once

Do not run AQuA2 on the same animal from both PCs at the same time (both
would write the same `_AQuA2.mat`). Different animals, or one PC running
MATLAB while the other does Step 1b / Step 6 on results that already exist,
are fine.

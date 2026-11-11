r"""Run a MATLAB script from Python and stream its output — what the GUI's
"Run in MATLAB now" buttons and the "Auto: tune + analyze" chain use.

    from afferent2p.matlab_run import find_matlab, run_script
    code = run_script(find_matlab(), r"D:\...\aqua_batch\run_aqua_batch.m", log=print)

MATLAB is started as  matlab -batch "run('<script>')"  (R2019a or newer; no
desktop, no splash, exits when the script ends, exit code 0 = no error) with
-logfile <script folder>\<script>_log.txt, which is tailed into `log` while it
runs — the same text MATLAB would print in its command window. On Windows the
launcher returns at once unless -wait is given, so -wait is added there.
"""
from __future__ import annotations
import glob
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

MATLAB_GUESSES = (r"C:\Program Files\MATLAB\R20*\bin\matlab.exe", r"D:\Program Files\MATLAB\R20*\bin\matlab.exe",
                  r"C:\MATLAB\R20*\bin\matlab.exe", r"D:\MATLAB\R20*\bin\matlab.exe", r"H:\MATLAB\R20*\bin\matlab.exe",
                  "/usr/local/MATLAB/R20*/bin/matlab", "/Applications/MATLAB_R20*.app/bin/matlab")


def find_matlab() -> str:
    """The newest MATLAB executable in the usual places, or on the PATH, or ''."""
    for pat in MATLAB_GUESSES:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[-1]
    for name in ("matlab", "matlab.exe"):
        w = shutil.which(name)
        if w:
            return w
    return ""


def _mat_quote(s: str) -> str:
    return str(s).replace("'", "''")


def matlab_command(matlab_exe: str, script: Path, logfile: Path | None = None) -> list[str]:
    cmd = [str(matlab_exe)]
    if sys.platform.startswith("win"):
        cmd.append("-wait")
    cmd += ["-nosplash"]
    if logfile is not None:
        cmd += ["-logfile", str(logfile)]
    cmd += ["-batch", f"run('{_mat_quote(script)}')"]
    return cmd


def run_script(matlab_exe: str, script: Path, log=print, stop: threading.Event | None = None,
               name: str | None = None, poll_s: float = 1.5) -> int:
    """Run one .m script to the end; returns MATLAB's exit code (0 = the script
    finished without an uncaught error; -1 = stopped by the user; -2 = could not
    start). Output is streamed to `log` from the -logfile (and from stdout when
    the logfile stays empty)."""
    script = Path(script)
    exe = (matlab_exe or "").strip().strip('"')
    if not exe or not Path(exe).exists():
        log(f"MATLAB executable not found ({exe or 'not set'}) — set it in the top row (Browse…), or run in MATLAB: run('{script}')")
        return -2
    if not script.exists():
        log(f"script not found: {script}")
        return -2
    name = name or script.stem
    logfile = script.parent / f"{script.stem}_log.txt"
    try:
        if logfile.exists():
            logfile.unlink()
    except OSError:
        pass
    cmd = matlab_command(exe, script, logfile)
    log(f"MATLAB: {' '.join(cmd[-2:])}    (log: {logfile})")
    t0 = time.time()
    try:
        win = sys.platform.startswith("win")
        proc = subprocess.Popen(cmd, cwd=str(script.parent), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, errors="replace", bufsize=1,
                                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0) if win else 0,
                                start_new_session=not win)
    except Exception as e:
        log(f"could not start MATLAB: {e}")
        return -2
    piped: list[str] = []

    def drain():
        try:
            for line in proc.stdout:        # type: ignore[union-attr]
                piped.append(line.rstrip("\n"))
        except Exception:
            pass
    threading.Thread(target=drain, daemon=True).start()
    pos = 0
    shown_from_log = False

    def tail():
        nonlocal pos, shown_from_log
        if not logfile.exists():
            return
        try:
            with open(logfile, "r", encoding="utf-8", errors="replace") as fh:
                fh.seek(pos)
                chunk = fh.read()
                pos = fh.tell()
        except OSError:
            return
        for line in chunk.splitlines():
            if line.strip():
                log("  | " + line.rstrip())
                shown_from_log = True
    while proc.poll() is None:
        tail()
        if stop is not None and stop.is_set():
            log(f"stopping MATLAB ({name}) …")
            kill_tree(proc)
            try:
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                try:
                    proc.kill()
                except Exception:
                    pass
            tail()
            return -1
        time.sleep(poll_s)
    time.sleep(0.3)
    tail()
    if not shown_from_log:
        for line in piped[-200:]:
            if line.strip():
                log("  | " + line)
    code = int(proc.returncode or 0)
    mins = (time.time() - t0) / 60.0
    log(f"MATLAB finished {name} in {mins:.1f} min" + ("" if code == 0 else f" with exit code {code} (an error — see the lines above)"))
    return code


def kill_tree(proc: subprocess.Popen):
    """Best-effort stop of a MATLAB process started by run_script (on Windows
    the launcher waits for a separate MATLAB process, so the whole tree goes)."""
    try:
        if sys.platform.startswith("win"):
            subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True)
        else:
            import os
            import signal
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)       # the whole process group (started with start_new_session)
    except Exception:
        try:
            proc.terminate()
        except Exception:
            pass


__all__ = ["find_matlab", "run_script", "matlab_command", "kill_tree", "MATLAB_GUESSES"]

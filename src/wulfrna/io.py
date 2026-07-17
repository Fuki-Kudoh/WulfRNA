from __future__ import annotations

import datetime as dt
import socket
import subprocess
import shutil
from pathlib import Path


def now_iso() -> str:
    return dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")


def git_commit() -> str:
    proc = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True)
    return proc.stdout.strip() if proc.returncode == 0 else "unknown"


def capture_versions_file(out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write(f"pipeline_version: git:{git_commit()}\n")
        f.write(f"datetime: {now_iso()}\n")
        f.write(f"hostname: {socket.gethostname()}\n")

        module_list = subprocess.run(["bash", "-lc", "module list 2>&1 || true"], capture_output=True, text=True)
        f.write("module_list:\n")
        f.write(module_list.stdout.strip() + "\n\n")

        version_cmds = {
            "fastqc": ["fastqc", "--version"],
            "cutadapt": ["cutadapt", "--version"],
            "salmon": ["salmon", "--version"],
            "kallisto": ["kallisto", "version"],
            "multiqc": ["multiqc", "--version"],
            "STAR": ["STAR", "--version"],
            "samtools": ["samtools", "--version"],
        }
        for name, cmd in version_cmds.items():
            if shutil.which(cmd[0]) is None:
                f.write(f"{name}: not found\n")
                continue
            p = subprocess.run(cmd, capture_output=True, text=True)
            txt = (p.stdout or p.stderr).strip().splitlines()
            ver = txt[0] if txt else "unknown"
            f.write(f"{name}: {ver}\n")

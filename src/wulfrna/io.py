from __future__ import annotations

import datetime as dt
import socket
import subprocess
import shutil
from importlib import metadata
from pathlib import Path
from typing import Optional


def now_iso() -> str:
    return dt.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z")


def package_version() -> str:
    try:
        return metadata.version("wulfrna")
    except metadata.PackageNotFoundError:
        return "unknown"


def find_repo_root(package_file: Optional[Path] = None) -> Optional[Path]:
    start = Path(__file__ if package_file is None else package_file).resolve()
    for parent in start.parents:
        if (parent / ".git").exists():
            return parent
    return None


def git_commit() -> Optional[str]:
    repo_root = find_repo_root()
    if repo_root is None:
        return None
    proc = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=repo_root, capture_output=True, text=True)
    return proc.stdout.strip() if proc.returncode == 0 and proc.stdout.strip() else None


def pipeline_version(package_ver: Optional[str] = None, commit: Optional[str] = None) -> str:
    resolved_package_ver = package_version() if package_ver is None else package_ver
    resolved_commit = git_commit() if commit is None else commit
    if resolved_commit is not None:
        return f"package:{resolved_package_ver} git:{resolved_commit}"
    return f"package:{resolved_package_ver}"


def capture_versions_file(out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    package_ver = package_version()
    commit = git_commit()
    with out.open("w", encoding="utf-8") as f:
        f.write(f"pipeline_version: {pipeline_version(package_ver, commit)}\n")
        f.write(f"package_version: {package_ver}\n")
        f.write(f"git_commit: {commit if commit is not None else 'not available'}\n")
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

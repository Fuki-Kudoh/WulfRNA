import subprocess
from pathlib import Path

from wulfrna import io


def init_git_repo(repo: Path) -> str:
    repo.mkdir()
    subprocess.run(["git", "init"], cwd=repo, check=True, capture_output=True, text=True)
    subprocess.run(["git", "config", "user.email", "test@example.com"], cwd=repo, check=True)
    subprocess.run(["git", "config", "user.name", "Test User"], cwd=repo, check=True)
    package_file = repo / "src" / "wulfrna" / "io.py"
    package_file.parent.mkdir(parents=True)
    package_file.write_text("# package file\n", encoding="utf-8")
    subprocess.run(["git", "add", "src/wulfrna/io.py"], cwd=repo, check=True)
    subprocess.run(["git", "commit", "-m", "initial"], cwd=repo, check=True, capture_output=True, text=True)
    result = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=repo, check=True, capture_output=True, text=True)
    return result.stdout.strip()


def test_git_commit_uses_package_repo_when_invoked_from_non_git_workdir(monkeypatch, tmp_path):
    repo = tmp_path / "repo"
    expected_commit = init_git_repo(repo)
    package_file = repo / "src" / "wulfrna" / "io.py"
    workdir = tmp_path / "analysis-workdir"
    workdir.mkdir()
    monkeypatch.setattr(io, "__file__", str(package_file))
    monkeypatch.chdir(workdir)

    assert io.git_commit() == expected_commit


def test_find_repo_root_detects_editable_repository_from_package_file(tmp_path):
    repo = tmp_path / "repo"
    init_git_repo(repo)
    package_file = repo / "src" / "wulfrna" / "io.py"

    assert io.find_repo_root(package_file) == repo


def test_capture_versions_falls_back_to_package_version_without_git(monkeypatch, tmp_path):
    package_file = tmp_path / "site-packages" / "wulfrna" / "io.py"
    package_file.parent.mkdir(parents=True)
    package_file.write_text("# installed package\n", encoding="utf-8")
    monkeypatch.setattr(io, "__file__", str(package_file))
    monkeypatch.setattr(io.metadata, "version", lambda name: "1.2.3")
    monkeypatch.setattr(io.shutil, "which", lambda name: None)

    versions = tmp_path / "versions.txt"
    io.capture_versions_file(versions)

    text = versions.read_text(encoding="utf-8")
    assert "pipeline_version: package:1.2.3" in text
    assert "package_version: 1.2.3" in text
    assert "git_commit: not available" in text

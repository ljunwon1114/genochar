from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any


CONFIG_VERSION = 1
DEFAULT_BASE_DIR = Path.home() / ".genochar"
DEFAULT_DB_DIR = DEFAULT_BASE_DIR / "databases" / "CheckM2_database"


class SetupError(RuntimeError):
    pass


@dataclass
class ManagedToolConfig:
    version: int = CONFIG_VERSION
    base_dir: str = str(DEFAULT_BASE_DIR)
    conda_exe: str | None = None
    prokka_prefix: str | None = None
    checkm2_prefix: str | None = None
    checkm2_db: str | None = None

    @property
    def base_path(self) -> Path:
        return Path(self.base_dir).expanduser().resolve()

    @property
    def config_path(self) -> Path:
        return self.base_path / "config.json"


def default_base_dir() -> Path:
    return DEFAULT_BASE_DIR


def _ensure_parent(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def detect_conda_executable(prefer_mamba: bool = False) -> str:
    candidates = []
    if prefer_mamba:
        candidates.extend([os.environ.get("MAMBA_EXE"), shutil.which("mamba")])
    candidates.extend([os.environ.get("CONDA_EXE"), shutil.which("conda")])
    for cand in candidates:
        if cand:
            return str(cand)
    raise SetupError(
        "Could not find conda or mamba on PATH. Install Conda/Mamba first, then rerun 'genochar setup'."
    )


def load_config(base_dir: Path | None = None) -> ManagedToolConfig | None:
    root = (base_dir or DEFAULT_BASE_DIR).expanduser().resolve()
    cfg_path = root / "config.json"
    if not cfg_path.exists():
        return None
    data = json.loads(cfg_path.read_text(encoding="utf-8"))
    return ManagedToolConfig(**data)



def save_config(config: ManagedToolConfig) -> Path:
    path = _ensure_parent(config.config_path)
    path.write_text(json.dumps(asdict(config), indent=2, sort_keys=True), encoding="utf-8")
    return path



def _find_dmnd(path: Path) -> Path | None:
    if path.is_file() and path.suffix == ".dmnd":
        return path.resolve()
    if path.is_dir():
        matches = sorted(path.rglob("*.dmnd"))
        if matches:
            return matches[0].resolve()
    return None



def normalize_checkm2_db(path: str | Path) -> Path:
    p = Path(path).expanduser().resolve()
    match = _find_dmnd(p)
    if match is None:
        raise SetupError(
            f"Could not locate a CheckM2 .dmnd database file under: {p}. "
            "Pass --checkm2-db with the .dmnd file or a directory containing it."
        )
    return match



def _run(cmd: list[str], verbose: bool = True, env: dict[str, str] | None = None) -> None:
    if verbose:
        print(">>", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, env=env)
    except subprocess.CalledProcessError as exc:
        raise SetupError(f"Command failed with exit code {exc.returncode}: {' '.join(cmd)}") from exc



def _create_env(prefix: Path, packages: list[str], *, verbose: bool = True, force: bool = False) -> None:
    if prefix.exists() and not force:
        if verbose:
            print(f"Using existing environment: {prefix}")
        return
    if prefix.exists() and force:
        shutil.rmtree(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    creator = detect_conda_executable(prefer_mamba=True)
    cmd = [
        creator,
        "create",
        "-y",
        "-p",
        str(prefix),
        "--override-channels",
        "-c",
        "conda-forge",
        "-c",
        "bioconda",
        *packages,
    ]
    _run(cmd, verbose=verbose)



def _download_checkm2_db(checkm2_prefix: Path, db_dir: Path, *, conda_exe: str, verbose: bool = True) -> Path:
    db_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        conda_exe,
        "run",
        "--no-capture-output",
        "-p",
        str(checkm2_prefix),
        "checkm2",
        "database",
        "--download",
        "--path",
        str(db_dir),
    ]
    _run(cmd, verbose=verbose)
    dmnd = _find_dmnd(db_dir)
    if dmnd is None:
        raise SetupError(f"CheckM2 database download completed, but no .dmnd file was found under: {db_dir}")
    return dmnd



def setup_managed_tools(
    *,
    base_dir: str | Path | None = None,
    with_prokka: bool = True,
    with_checkm2: bool = True,
    checkm2_db: str | Path | None = None,
    force: bool = False,
    verbose: bool = True,
) -> ManagedToolConfig:
    root = (Path(base_dir).expanduser().resolve() if base_dir else DEFAULT_BASE_DIR)
    env_root = root / "envs"
    prokka_prefix = env_root / "prokka"
    checkm2_prefix = env_root / "checkm2"
    conda_exe = detect_conda_executable(prefer_mamba=False)

    cfg = ManagedToolConfig(
        base_dir=str(root),
        conda_exe=conda_exe,
        prokka_prefix=str(prokka_prefix) if with_prokka else None,
        checkm2_prefix=str(checkm2_prefix) if with_checkm2 else None,
        checkm2_db=None,
    )

    if with_prokka:
        _create_env(prokka_prefix, ["prokka"], verbose=verbose, force=force)

    if with_checkm2:
        _create_env(checkm2_prefix, ["python=3.12", "checkm2"], verbose=verbose, force=force)
        if checkm2_db is not None:
            db_path = normalize_checkm2_db(checkm2_db)
        elif os.environ.get("CHECKM2DB"):
            db_path = normalize_checkm2_db(os.environ["CHECKM2DB"])
        else:
            existing = _find_dmnd(DEFAULT_DB_DIR if root == DEFAULT_BASE_DIR else root / "databases" / "CheckM2_database")
            if existing is not None:
                db_path = existing
            else:
                db_dir = root / "databases" / "CheckM2_database"
                db_path = _download_checkm2_db(checkm2_prefix, db_dir, conda_exe=conda_exe, verbose=verbose)
        cfg.checkm2_db = str(db_path)

    save_config(cfg)
    return cfg

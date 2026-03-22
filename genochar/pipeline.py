from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import List, Sequence

from .managed_envs import load_config
from .utils import infer_strain_name


class PipelineError(RuntimeError):
    pass


_TOOL_TO_PREFIX_KEY = {
    "prokka": "prokka_prefix",
    "checkm2": "checkm2_prefix",
}



def _resolve_tool_invocation(name: str) -> tuple[list[str], dict[str, str] | None]:
    direct = shutil.which(name)
    if direct is not None:
        env = os.environ.copy()
        cfg = load_config()
        if name == "checkm2" and cfg and cfg.checkm2_db:
            env["CHECKM2DB"] = cfg.checkm2_db
        return [name], env

    cfg = load_config()
    if cfg is not None:
        prefix_attr = _TOOL_TO_PREFIX_KEY.get(name)
        prefix = getattr(cfg, prefix_attr, None) if prefix_attr else None
        conda_exe = cfg.conda_exe or os.environ.get("CONDA_EXE") or shutil.which("conda")
        if prefix and conda_exe:
            env = os.environ.copy()
            if name == "checkm2" and cfg.checkm2_db:
                env["CHECKM2DB"] = cfg.checkm2_db
            return [conda_exe, "run", "--no-capture-output", "-p", str(prefix), name], env

    if name == "checkm2":
        raise PipelineError(
            "Required executable not found for CheckM2. Run 'genochar setup' or install checkm2 on PATH."
        )
    if name == "prokka":
        raise PipelineError(
            "Required executable not found for Prokka. Run 'genochar setup' or install prokka on PATH."
        )
    raise PipelineError(f"Required executable not found on PATH: {name}")



def run_command(cmd: Sequence[str], verbose: bool = True, env: dict[str, str] | None = None) -> None:
    if verbose:
        print(">>", " ".join(shlex.quote(c) for c in cmd))
    try:
        subprocess.run(list(cmd), check=True, env=env)
    except subprocess.CalledProcessError as exc:
        raise PipelineError(
            f"Command failed with exit code {exc.returncode}: {' '.join(shlex.quote(c) for c in cmd)}"
        ) from exc



def run_prokka(
    assemblies: Sequence[Path],
    outdir: Path,
    threads: int = 8,
    kingdom: str | None = None,
    extra_args: str | None = None,
    verbose: bool = True,
) -> List[Path]:
    prokka_cmd, env = _resolve_tool_invocation("prokka")
    outdir.mkdir(parents=True, exist_ok=True)

    gffs: List[Path] = []
    for asm in assemblies:
        strain = infer_strain_name(asm)
        sample_out = outdir / strain
        sample_out.mkdir(parents=True, exist_ok=True)

        cmd = [
            *prokka_cmd,
            "--outdir", str(sample_out),
            "--prefix", strain,
            "--force",
            "--cpus", str(threads),
        ]
        if kingdom and kingdom.lower() != "auto":
            cmd += ["--kingdom", kingdom]
        if extra_args:
            cmd += shlex.split(extra_args)
        cmd.append(str(asm))
        run_command(cmd, verbose=verbose, env=env)
        gff = sample_out / f"{strain}.gff"
        gffs.append(gff)

    return gffs



def _resolve_checkm2_database_path() -> str | None:
    cfg = load_config()
    if cfg is not None and cfg.checkm2_db:
        return cfg.checkm2_db
    return os.environ.get("CHECKM2DB")



def run_checkm2(
    assemblies: Sequence[Path],
    outdir: Path,
    threads: int = 8,
    verbose: bool = True,
) -> Path:
    checkm2_cmd, env = _resolve_tool_invocation("checkm2")
    outdir.mkdir(parents=True, exist_ok=True)

    assembly_inputs = [str(Path(asm).resolve()) for asm in assemblies]
    if not assembly_inputs:
        raise PipelineError("No assembly FASTA files were provided to CheckM2.")

    cmd = [
        *checkm2_cmd,
        "predict",
        "--threads", str(threads),
        "--input",
        *assembly_inputs,
        "--output-directory", str(outdir),
        "--force",
    ]

    db_path = _resolve_checkm2_database_path()
    if db_path:
        cmd += ["--database_path", db_path]

    run_command(cmd, verbose=verbose, env=env)

    report = outdir / "quality_report.tsv"
    if not report.exists():
        raise PipelineError(f"CheckM2 finished without producing: {report}")
    return report



def find_existing_gffs(assemblies: Sequence[Path], gff_inputs: Sequence[Path] | None = None) -> List[Path]:
    if gff_inputs:
        return list(gff_inputs)

    found: List[Path] = []
    for asm in assemblies:
        strain = infer_strain_name(asm)
        candidates = [
            asm.with_suffix(".gff"),
            asm.with_suffix(".gff3"),
            asm.parent / f"{strain}.gff",
            asm.parent / f"{strain}.gff3",
            asm.parent / f"{strain}_prokka" / f"{strain}.gff",
            asm.parent / f"{strain}_bakta" / f"{strain}.gff3",
        ]
        hit = next((p.resolve() for p in candidates if p.exists()), None)
        if hit is None:
            raise PipelineError(
                f"Could not locate an existing GFF/GFF3 for assembly: {asm}. "
                "Provide --gff explicitly or choose --annotate prokka/none."
            )
        found.append(hit)
    return found

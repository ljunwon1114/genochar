from __future__ import annotations

import shlex
import shutil
import subprocess
from pathlib import Path
from typing import List, Sequence

from .utils import infer_strain_name


class PipelineError(RuntimeError):
    pass


def require_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise PipelineError(f"Required executable not found on PATH: {name}")


def run_command(cmd: Sequence[str], verbose: bool = True) -> None:
    if verbose:
        print(">>", " ".join(shlex.quote(c) for c in cmd))
    try:
        subprocess.run(list(cmd), check=True)
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
    require_executable("prokka")
    outdir.mkdir(parents=True, exist_ok=True)

    gffs: List[Path] = []
    for asm in assemblies:
        strain = infer_strain_name(asm)
        sample_out = outdir / strain
        sample_out.mkdir(parents=True, exist_ok=True)

        cmd = [
            "prokka",
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
        run_command(cmd, verbose=verbose)
        gff = sample_out / f"{strain}.gff"
        gffs.append(gff)

    return gffs


def run_bakta(
    assemblies: Sequence[Path],
    outdir: Path,
    threads: int = 8,
    extra_args: str | None = None,
    verbose: bool = True,
) -> List[Path]:
    require_executable("bakta")
    outdir.mkdir(parents=True, exist_ok=True)

    gffs: List[Path] = []
    for asm in assemblies:
        strain = infer_strain_name(asm)
        sample_out = outdir / strain
        sample_out.mkdir(parents=True, exist_ok=True)

        cmd = [
            "bakta",
            "--output", str(sample_out),
            "--prefix", strain,
            "--keep-contig-headers",
            "--force",
            "--threads", str(threads),
        ]
        if extra_args:
            cmd += shlex.split(extra_args)
        cmd.append(str(asm))
        run_command(cmd, verbose=verbose)

        candidate_gffs = [
            sample_out / f"{strain}.gff3",
            sample_out / f"{strain}.gff",
        ]
        gff = next((p for p in candidate_gffs if p.exists()), candidate_gffs[0])
        gffs.append(gff)

    return gffs


def _prepare_checkm2_inputs(assemblies: Sequence[Path], input_dir: Path) -> None:
    if input_dir.exists():
        shutil.rmtree(input_dir)
    input_dir.mkdir(parents=True, exist_ok=True)

    for asm in assemblies:
        strain = infer_strain_name(asm)
        staged = input_dir / f"{strain}.fa"
        source = Path(asm).resolve()
        if staged.exists() or staged.is_symlink():
            staged.unlink()
        try:
            staged.symlink_to(source)
        except OSError:
            shutil.copy2(source, staged)


def run_checkm2(
    assemblies: Sequence[Path],
    outdir: Path,
    threads: int = 8,
    verbose: bool = True,
) -> Path:
    require_executable("checkm2")
    outdir.mkdir(parents=True, exist_ok=True)

    input_dir = outdir / "input_bins"
    _prepare_checkm2_inputs(assemblies, input_dir)

    cmd = [
        "checkm2",
        "predict",
        "--threads", str(threads),
        "--input", str(input_dir),
        "--output-directory", str(outdir),
        "--force",
    ]
    run_command(cmd, verbose=verbose)

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
                "Provide --gff explicitly or choose --annotate prokka/bakta/none."
            )
        found.append(hit)
    return found

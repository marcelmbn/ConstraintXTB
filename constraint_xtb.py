from __future__ import annotations

from pathlib import Path
import subprocess as sp
from typing import Iterable, Literal

from pydantic import BaseModel, PositiveInt, field_validator

ConstraintValue = float | Literal["auto"]
ConstraintValueInput = float | int | str


def coerce_constraint_value(value: ConstraintValueInput) -> ConstraintValue:
    """
    Normalise a constraint value to either ``float`` or the literal ``'auto'``.

    Raises:
        ValueError: If a string different from ``'auto'`` is provided.
        TypeError: If the value cannot be converted to a floating point number.
    """
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized != "auto":
            raise ValueError("String constraint values must be 'auto'.")
        return "auto"
    try:
        return float(value)
    except (TypeError, ValueError) as exc:  # pragma: no cover - defensive programming
        raise TypeError("Constraint value must be numeric or 'auto'.") from exc


class ConstraintBase(BaseModel):
    """Common functionality shared by the constraint data models."""

    value: ConstraintValue

    @field_validator("value", mode="before")
    @classmethod
    def _validate_value(cls, value: ConstraintValueInput) -> ConstraintValue:
        return coerce_constraint_value(value)


class DistanceConstraint(ConstraintBase):
    """Representation of a distance constraint between two atoms."""

    atom1: PositiveInt
    atom2: PositiveInt


class AngleConstraint(ConstraintBase):
    """Representation of an angle constraint between three atoms."""

    atom1: PositiveInt
    atom2: PositiveInt
    atom3: PositiveInt


class DihedralConstraint(ConstraintBase):
    """Representation of a dihedral constraint between four atoms."""

    atom1: PositiveInt
    atom2: PositiveInt
    atom3: PositiveInt
    atom4: PositiveInt


class DetailedInput:
    """
    Build xTB detailed (xcontrol) input files supporting distance, angle and dihedral constraints.

    The class collects constraints and serialises them to the `$constrain` block required
    by xTB. Constraints can be provided up-front via the constructor or appended later
    using the respective ``add_*`` methods.
    """

    def __init__(
        self,
        *,
        force_constant: float | None = None,
        distances: Iterable[DistanceConstraint | tuple[int, int, ConstraintValueInput]]
        | None = None,
        angles: Iterable[AngleConstraint | tuple[int, int, int, ConstraintValueInput]]
        | None = None,
        dihedrals: Iterable[
            DihedralConstraint | tuple[int, int, int, int, ConstraintValueInput]
        ]
        | None = None,
    ) -> None:
        self.force_constant: float | None = force_constant
        self.distances: list[DistanceConstraint] = []
        self.angles: list[AngleConstraint] = []
        self.dihedrals: list[DihedralConstraint] = []

        if distances:
            for distance in distances:
                if isinstance(distance, DistanceConstraint):
                    self.distances.append(distance)
                else:
                    self.add_distance_constraint(*distance)
        if angles:
            for angle in angles:
                if isinstance(angle, AngleConstraint):
                    self.angles.append(angle)
                else:
                    self.add_angle_constraint(*angle)
        if dihedrals:
            for dihedral in dihedrals:
                if isinstance(dihedral, DihedralConstraint):
                    self.dihedrals.append(dihedral)
                else:
                    self.add_dihedral_constraint(*dihedral)

    def add_distance_constraint(
        self,
        atom1: int,
        atom2: int,
        value: ConstraintValueInput,
    ) -> DetailedInput:
        """Add a distance constraint to the detailed input.

        Parameters
        ----------
        atom1 : int
            Index (1-based) of the first atom.
        atom2 : int
            Index (1-based) of the second atom.
        value : ConstraintValueInput
            Target distance in ångström or the literal ``'auto'``.

        Returns
        -------
        DetailedInput
            The current instance for fluent chaining.
        """
        normalized_value = coerce_constraint_value(value)
        constraint = DistanceConstraint(
            atom1=atom1, atom2=atom2, value=normalized_value
        )
        self.distances.append(constraint)
        return self

    def add_angle_constraint(
        self,
        atom1: int,
        atom2: int,
        atom3: int,
        value: ConstraintValueInput,
    ) -> DetailedInput:
        """Add an angle constraint to the detailed input.

        Parameters
        ----------
        atom1 : int
            Index of the first atom.
        atom2 : int
            Index of the vertex atom.
        atom3 : int
            Index of the third atom.
        value : ConstraintValueInput
            Target angle in degrees or ``'auto'``.

        Returns
        -------
        DetailedInput
            The current instance for fluent chaining.
        """
        normalized_value = coerce_constraint_value(value)
        constraint = AngleConstraint(
            atom1=atom1, atom2=atom2, atom3=atom3, value=normalized_value
        )
        self.angles.append(constraint)
        return self

    def add_dihedral_constraint(
        self,
        atom1: int,
        atom2: int,
        atom3: int,
        atom4: int,
        value: ConstraintValueInput,
    ) -> DetailedInput:
        """Add a dihedral constraint to the detailed input.

        Parameters
        ----------
        atom1 : int
            Index of the first atom.
        atom2 : int
            Index of the second atom.
        atom3 : int
            Index of the third atom.
        atom4 : int
            Index of the fourth atom.
        value : ConstraintValueInput
            Target dihedral angle in degrees or ``'auto'``.

        Returns
        -------
        DetailedInput
            The current instance for fluent chaining.
        """
        normalized_value = coerce_constraint_value(value)
        constraint = DihedralConstraint(
            atom1=atom1,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            value=normalized_value,
        )
        self.dihedrals.append(constraint)
        return self

    @staticmethod
    def _format_numeric(value: float) -> str:
        """Convert numbers to compact strings while keeping a decimal point if needed."""
        formatted = f"{value:.10f}".rstrip("0").rstrip(".")
        return formatted if formatted else "0"

    def _format_value(self, value: ConstraintValue) -> str:
        """Return the textual representation for a constraint value."""
        if isinstance(value, str):
            return value
        return self._format_numeric(float(value))

    def _has_constraints(self) -> bool:
        """Return True when at least one constraint or force constant is defined."""
        return bool(
            self.distances
            or self.angles
            or self.dihedrals
            or self.force_constant is not None
        )

    def compile(self) -> str:
        """Serialise the stored constraints to an xTB `$constrain` block."""
        if not self._has_constraints():
            return ""

        lines: list[str] = ["$constrain"]

        if self.force_constant is not None:
            lines.append(
                f"   force constant={self._format_numeric(self.force_constant)}"
            )

        for distance in self.distances:
            lines.append(
                f"   distance: {distance.atom1}, {distance.atom2}, {self._format_value(distance.value)}"
            )

        for angle in self.angles:
            lines.append(
                f"   angle: {angle.atom1}, {angle.atom2}, {angle.atom3}, {self._format_value(angle.value)}"
            )

        for dihedral in self.dihedrals:
            lines.append(
                "   dihedral: "
                f"{dihedral.atom1}, {dihedral.atom2}, {dihedral.atom3}, {dihedral.atom4}, "
                f"{self._format_value(dihedral.value)}"
            )

        lines.append("$end")
        return "\n".join(lines) + "\n"

    def write(self, detailed_input_path: Path = Path("detailed_input.inp")) -> Path:
        """Persist the compiled detailed input to disk and return the resulting path."""
        detailed_input_path.write_text(self.compile())
        return detailed_input_path


class XTB:
    """
    Class to run xtb calculations.
    """

    def __init__(self, path: Path, detailedinput: DetailedInput) -> None:
        self.path = path
        self.detailedinput = detailedinput
        # mapping from gfn method to cli argument
        self.dict_gfn_map = {
            "gfn0": "--gfn 0",
            "gfn1": "--gfn 1",
            "gfn2": "--gfn 2",
            "gfnff": "--gfnff",
        }

    def _run(self, temp_path: Path, arguments: list[str]) -> tuple[str, str, int]:
        """
        Run xtb with the given arguments.

        Arguments:
        arguments (list[str]): The arguments to pass to xtb.

        Returns:
        tuple[str, str, int]: The output of the xtb calculation (stdout and stderr)
                              and the return code
        """
        try:
            xtb_out = sp.run(
                [str(self.path)] + arguments,
                cwd=temp_path,
                capture_output=True,
                check=True,
            )
            # get the output of the xtb calculation (of both stdout and stderr)
            xtb_log_out = xtb_out.stdout.decode("utf8", errors="replace")
            xtb_log_err = xtb_out.stderr.decode("utf8", errors="replace")
            return xtb_log_out, xtb_log_err, 0
        except sp.CalledProcessError as e:
            xtb_log_out = e.stdout.decode("utf8", errors="replace")
            xtb_log_err = e.stderr.decode("utf8", errors="replace")
            return xtb_log_out, xtb_log_err, e.returncode

    def run_xtb(
        self,
        input_structure: str,
        gfn_method: str,
        ncores: int,
        *,
        verbose: bool = False,
        workdir: Path | None = None,
    ) -> int:
        """Run an xtb optimization with the given number of cores.

        Parameters
        ----------
        input_structure : str
            Path to the coordinate file that should be optimised.
        gfn_method : str
            Identifier of the GFN method to use.
        ncores : int
            Number of CPU cores.
        verbose : bool, optional
            When ``True`` print the captured stdout/stderr streams.
        workdir : Path | None, optional
            Directory used as working directory for the xtb execution. Defaults to the
            current directory when ``None``.

        Returns
        -------
        int
            Return code of the xtb execution.
        """

        working_directory = workdir or Path("xtb_constraint_opt")
        working_directory.mkdir(parents=True, exist_ok=True)
        detailed_input_file = self.detailedinput.write(
            working_directory / "detailed_input.inp"
        )

        try:
            gfn_flag = self.dict_gfn_map[gfn_method.lower()]
        except KeyError as exc:
            available = ", ".join(self.dict_gfn_map)
            raise ValueError(
                f"Unsupported GFN method '{gfn_method}'. Choose one of: {available}."
            ) from exc

        arguments = [
            input_structure,
            "--opt",
            gfn_flag,
            "-P",
            f"{ncores}",
            "--input",
            str(detailed_input_file.name),
        ]
        xtb_out, xtb_err, return_code = self._run(
            temp_path=working_directory, arguments=arguments
        )

        if verbose:
            if xtb_out:
                print("--- xtb stdout ---")
                print(xtb_out.rstrip())
            if xtb_err:
                print("--- xtb stderr ---")
                print(xtb_err.rstrip())
            if not xtb_out and not xtb_err:
                print("--- xtb produced no output ---")

        return return_code


def run_constraint_opt(
    structure: Path | str,
    *,
    xtb_path: Path | str = Path("xtb"),
    gfn_method: str = "gfn2",
    ncores: int = 4,
    force_constant: float | None = None,
    distances: Iterable[tuple[int, int, ConstraintValueInput]] | None = None,
    angles: Iterable[tuple[int, int, int, ConstraintValueInput]] | None = None,
    dihedrals: Iterable[tuple[int, int, int, int, ConstraintValueInput]] | None = None,
    detailed_input: DetailedInput | None = None,
    verbose: bool = False,
    workdir: Path | None = None,
) -> int:
    """
    Run an xTB optimisation with optional distance, angle and dihedral constraints.

    Parameters
    ----------
    structure : Path | str
        Path to the structure file (e.g. xyz) that should be optimised.
    xtb_path : Path | str, optional
        Path to the xtb executable. Defaults to ``Path("xtb")``.
    gfn_method : str, optional
        Identifier of the GFN method to use (``gfn0``, ``gfn1``, ``gfn2``, ``gfnff``). Defaults to ``"gfn2"``.
    ncores : int, optional
        Number of CPU cores to pass via ``-P``. Defaults to ``4``.
    force_constant : float | None, optional
        Global force constant applied to the constraints. Defaults to ``None``.
    distances : Iterable[tuple[int, int, ConstraintValueInput]] | None, optional
        Iterable of ``(atom_i, atom_j, value)`` tuples defining distance constraints. Defaults to ``None``.
    angles : Iterable[tuple[int, int, int, ConstraintValueInput]] | None, optional
        Iterable of ``(atom_i, atom_j, atom_k, value)`` tuples defining angle constraints. Defaults to ``None``.
    dihedrals : Iterable[tuple[int, int, int, int, ConstraintValueInput]] | None, optional
        Iterable of ``(atom_i, atom_j, atom_k, atom_l, value)`` tuples defining dihedral constraints. Defaults to ``None``.
    detailed_input : DetailedInput | None, optional
        Preconfigured ``DetailedInput`` instance. When provided, additional constraint arguments are merged into this instance.
        Defaults to ``None``.
    verbose : bool, optional
        When ``True``, print the captured xtb stdout/stderr streams. Defaults to ``False``.
    workdir : Path | None, optional
        Directory used as working directory for the xtb execution. Defaults to the current directory.

    Returns
    -------
    int
        Return code of the xtb execution.
    """
    if isinstance(structure, Path):
        structure = str(structure.resolve())

    if detailed_input is None:
        detailedinput = DetailedInput(
            force_constant=force_constant,
            distances=distances,
            angles=angles,
            dihedrals=dihedrals,
        )
    else:
        detailedinput = detailed_input
        if force_constant is not None:
            detailedinput.force_constant = force_constant
        if distances:
            for distance in distances:
                detailedinput.add_distance_constraint(*distance)
        if angles:
            for angle in angles:
                detailedinput.add_angle_constraint(*angle)
        if dihedrals:
            for dihedral in dihedrals:
                detailedinput.add_dihedral_constraint(*dihedral)

    workdir_path = Path(workdir) if workdir is not None else None
    xtb = XTB(path=Path(xtb_path), detailedinput=detailedinput)
    return xtb.run_xtb(
        input_structure=structure,
        gfn_method=gfn_method,
        ncores=ncores,
        verbose=verbose,
        workdir=workdir_path,
    )

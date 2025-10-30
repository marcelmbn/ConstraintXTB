from __future__ import annotations

from pathlib import Path

from constraint_xtb import DetailedInput, run_constraint_opt


WATER_GEOMETRY = """3
water molecule
O      0.000000     0.000000     0.000000
H      0.757160     0.586260     0.000000
H     -0.757160     0.586260     0.000000
"""

DEMO_DIRECTORY = Path("demo_output")


def write_water_xyz(path: Path) -> Path:
    """
    Persist a minimal water geometry to the provided XYZ path.
    """
    path.write_text(WATER_GEOMETRY)
    return path


def build_water_constraints() -> DetailedInput:
    """
    Create a DetailedInput instance constraining the H–O–H angle to 170°.
    """
    detailed_input = DetailedInput()
    detailed_input.add_angle_constraint(2, 1, 3, 170.0)
    detailed_input.force_constant = 0.5
    return detailed_input


def main() -> None:
    """
    Demonstrate running a constrained optimisation on water.
    """
    DEMO_DIRECTORY.mkdir(parents=True, exist_ok=True)

    structure_path = write_water_xyz(DEMO_DIRECTORY / "water.xyz")
    detailed_input = build_water_constraints()

    print(f"Wrote water geometry to {structure_path}")

    try:
        exit_code = run_constraint_opt(
            structure_path,
            detailed_input=detailed_input,
            verbose=True,
            workdir=DEMO_DIRECTORY,
        )
    except FileNotFoundError as exc:
        print("xtb executable not found. Skipping optimisation.\n", exc)
        return

    print(f"xTB finished with exit code {exit_code}")


if __name__ == "__main__":
    main()

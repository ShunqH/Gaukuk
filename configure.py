import argparse
from pathlib import Path


CONFIG_FILE = Path("config.mk")


def write_config(setup: str, eos: str, flux: str, use_openmp: bool):
    """
    Write configuration options to config.mk.
    """
    content = f"""# configurable options
SETUP = {setup}
EOS = {eos}
FLUX = {flux}
USE_OPENMP = {1 if use_openmp else 0}
"""

    CONFIG_FILE.write_text(content)


def main():
    parser = argparse.ArgumentParser(
        description="Configure build options for Gaukuk"
    )

    # options with values
    parser.add_argument(
        "--setup",
        default="kh",
        help="Setup problem name (default: kh)"
    )

    parser.add_argument(
        "--eos",
        default="adiabatic",
        help="Equation of state (default: adiabatic)"
    )

    parser.add_argument(
        "--flux",
        default="hllc",
        help="Flux solver (default: hllc)"
    )

    # switch option
    parser.add_argument(
        "-openmp",
        action="store_true",
        help="Enable OpenMP"
    )

    args = parser.parse_args()

    write_config(
        setup=args.setup,
        eos=args.eos,
        flux=args.flux,
        use_openmp=args.openmp
    )

    print("Configuration written to config.mk")
    print(f"SETUP       = {args.setup}")
    print(f"EOS         = {args.eos}")
    print(f"FLUX        = {args.flux}")
    print(f"USE_OPENMP  = {1 if args.openmp else 0}")


if __name__ == "__main__":
    main()
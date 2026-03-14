"""ANDES CLI entry point."""

from .main import create_cli
from ..core.andes import ANDES_BANDS, SUBSLIT_CHOICES

cli = create_cli('ANDES', ANDES_BANDS, SUBSLIT_CHOICES)


def main():
    cli()

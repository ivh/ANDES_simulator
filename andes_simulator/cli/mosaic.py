"""MOSAIC CLI entry point."""

from .main import create_cli
from ..core.mosaic import MOSAIC_BANDS, SUBSLIT_CHOICES

cli = create_cli('MOSAIC', MOSAIC_BANDS, SUBSLIT_CHOICES)


def main():
    cli()

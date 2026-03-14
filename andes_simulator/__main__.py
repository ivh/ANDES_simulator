"""
Main entry point for simulation framework.

Allows the package to be run as: python -m andes_simulator
Defaults to ANDES CLI.
"""

from .cli.andes import main

if __name__ == '__main__':
    main()

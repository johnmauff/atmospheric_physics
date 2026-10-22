"""Adiabatic suite: Python version of suites/suite_adiabatic.xml.

No physics: only the global energy fixer and diagnostics (for testing a
dynamical core with the physics framework around it).
"""

from xdsl_ccpp.frontend.py_api import ccpp_suite, export_suite

from categories.diagnostics import state_diagnostics, tendency_diagnostics
from categories.energy import energy_fixer, save_total_energy


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("adiabatic", version="1.0")
class AdiabaticSuite:
    """The 'adiabatic' suite."""

    def physics_before_coupler():
        """Physics before the coupler."""
        energy_fixer(diagnostics=True)
        save_total_energy()
        state_diagnostics()

    def physics_after_coupler():
        """Physics after the coupler."""
        tendency_diagnostics()


if __name__ == "__main__":
    export_suite(AdiabaticSuite)

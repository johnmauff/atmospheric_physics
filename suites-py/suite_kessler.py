"""Kessler suite: Python version of suites/suite_kessler.xml.

Moist dynamical-core test with Kessler warm-rain microphysics.
"""

from xdsl_ccpp.frontend.py_api import ccpp_suite, export_suite

from categories.diagnostics import state_diagnostics, tendency_diagnostics
from categories.energy import dycore_energy_consistency
from categories.idealized import idealized_physics
from parameterizations.kessler import kessler


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("kessler", version="1.0")
class KesslerSuite:
    """The 'kessler' suite."""

    def physics_before_coupler():
        """Physics before the coupler."""
        idealized_physics(scheme=kessler, check_energy=True)
        state_diagnostics(parameterizations=[kessler])

    def physics_after_coupler():
        """Physics after the coupler."""
        # Not yet: thermo_water_update, dme_adjust (atmospheric_physics#222).
        dycore_energy_consistency()
        tendency_diagnostics()


if __name__ == "__main__":
    export_suite(KesslerSuite)

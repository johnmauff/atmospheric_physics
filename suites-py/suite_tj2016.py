"""TJ2016 suite: Python version of suites/suite_tj2016.xml.

Moist dynamical-core test with Thatcher-Jablonowski (2016) simple physics:
large-scale precipitation before the coupler, surface fluxes and boundary
layer after it.
"""

from xdsl_ccpp.frontend.py_api import ccpp_suite, export_suite

from categories.diagnostics import state_diagnostics, tendency_diagnostics
from categories.energy import dycore_energy_consistency
from categories.idealized import idealized_physics
from parameterizations.tj2016 import tj2016_precipitation, tj2016_surface_pbl


# In a CCPP suite class, each method is a group recorded by @ccpp_suite,
# not a real method: it takes no 'self', and the class holds only groups.
# pylint: disable=no-method-argument,too-few-public-methods
@ccpp_suite("tj2016", version="1.0")
class Tj2016Suite:
    """The 'tj2016' suite."""

    def physics_before_coupler():
        """Physics before the coupler."""
        idealized_physics(scheme=tj2016_precipitation, check_energy=True)
        state_diagnostics()

    def physics_after_coupler():
        """Physics after the coupler."""
        idealized_physics(scheme=tj2016_surface_pbl)
        # Not yet: thermo_water_update, dme_adjust (atmospheric_physics#222).
        dycore_energy_consistency()
        tendency_diagnostics()


if __name__ == "__main__":
    export_suite(Tj2016Suite)

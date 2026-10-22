"""Category: idealized (simple) physics -- a single parameterization standing
in for the full physics, for dynamical-core tests."""

from xdsl_ccpp.frontend.py_api import ccpp_composite

from categories.energy import energy_check_without_fluxes
from categories.interface import run_setup


@ccpp_composite
def idealized_physics(*, scheme, check_energy=False):
    """Roles:
        scheme: the idealized parameterization (required).
    check_energy: follow it with the energy and water budget check (with
    zero fluxes: simple physics provides none).
    """
    run_setup(scheme)
    scheme.run()
    if check_energy:
        energy_check_without_fluxes()

"""Category: convection."""

from xdsl_ccpp.frontend.py_api import ccpp_composite

from categories.energy import energy_checked
from categories.interface import run_setup


@ccpp_composite
def convection(*, dry_adjustment=None, deep=None, shallow=None):
    """Convective adjustment of the column, in three roles, in this order.

    Roles (each a Parameterization, or None when not used):
        dry_adjustment: removes dry static instability
        deep:           deep convection
        shallow:        shallow convection
    Deep and shallow convection each run inside an energy and water budget
    check (with their own flux correction, if they have one); dry
    adjustment does not.  Each role's setup (e.g. its namelist) runs just
    before it.
    """
    if dry_adjustment is not None:
        run_setup(dry_adjustment)
        dry_adjustment.run()
    for role in (deep, shallow):
        if role is not None:
            run_setup(role)
            energy_checked(process=role.run, flux_correction=role.energy_flux_correction)

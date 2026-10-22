"""Energy categories, and the energy and water budget check other categories use.

``energy_checked`` is the check CAM puts around a process: zero the boundary
fluxes, run the process, then check that the change in energy and water
matches the fluxes.  Categories (convection, gravity wave drag) and some
parameterizations (Rasch-Kristjansson) wrap their steps in it.
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import nothing
from categories.state import apply_tendency_of_air_temperature, update_temperature

check_energy_gmean               = ccpp_scheme_ref("check_energy_gmean")
check_energy_gmean_diagnostics   = ccpp_scheme_ref("check_energy_gmean_diagnostics")
check_energy_zero_fluxes         = ccpp_scheme_ref("check_energy_zero_fluxes")
check_energy_fix                 = ccpp_scheme_ref("check_energy_fix")
check_energy_scaling             = ccpp_scheme_ref("check_energy_scaling")
check_energy_chng                = ccpp_scheme_ref("check_energy_chng")
check_energy_fix_diagnostics     = ccpp_scheme_ref("check_energy_fix_diagnostics")
check_energy_save_teout          = ccpp_scheme_ref("check_energy_save_teout")
dycore_energy_consistency_adjust = ccpp_scheme_ref("dycore_energy_consistency_adjust")


@ccpp_composite
def energy_checked(*, process, flux_correction=None):
    """Run *process* inside an energy and water budget check.

    flux_correction: run after the process, before the check, when the
    process's boundary fluxes need adjusting first.
    """
    check_energy_zero_fluxes()
    process()
    (flux_correction or nothing)()
    check_energy_scaling()
    check_energy_chng()


@ccpp_composite
def energy_fixer(*, diagnostics=False):
    """Add the global mean energy difference (between this step's input and
    the previous step's output) as a uniform heating.

    diagnostics: also write the fixer's history output (after the check,
    which updates the global energies).
    """
    check_energy_gmean()
    check_energy_gmean_diagnostics()

    def fix():
        check_energy_fix()
        update_temperature()

    energy_checked(process=fix)
    if diagnostics:
        check_energy_fix_diagnostics()


@ccpp_composite
def energy_check_without_fluxes():
    """Energy and water budget check for simple physics, after it has run:
    simple physics provides no surface heat fluxes, so the check runs with
    zero fluxes."""
    check_energy_zero_fluxes()
    check_energy_scaling()
    check_energy_chng()


@ccpp_composite
def save_total_energy():
    """Save the current total energy, for the next step's energy fixer."""
    check_energy_save_teout()


@ccpp_composite
def dycore_energy_consistency():
    """Scale temperature for energy consistency with the dycore (MPAS and
    SE)."""
    check_energy_scaling()
    dycore_energy_consistency_adjust()
    apply_tendency_of_air_temperature()


@ccpp_composite
def end_of_step_energy():
    """Save total energy for the next step's fixer, then scale temperature
    for energy consistency with the dycore."""
    save_total_energy()
    dycore_energy_consistency()

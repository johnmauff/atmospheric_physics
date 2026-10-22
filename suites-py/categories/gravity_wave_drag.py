"""Category: gravity wave drag."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.energy import energy_checked
from categories.state import (
    apply_constituent_tendencies,
    apply_heating_rate,
    apply_tendency_of_eastward_wind,
    apply_tendency_of_northward_wind,
    geopotential_temp,
    qneg,
    update_dry_static_energy,
)

gravity_wave_drag_common             = ccpp_scheme_ref("gravity_wave_drag_common")
gravity_wave_drag_prepare_profiles   = ccpp_scheme_ref("gravity_wave_drag_prepare_profiles")
gravity_wave_drag_top_taper          = ccpp_scheme_ref("gravity_wave_drag_top_taper")
convert_dry_constituent_tendencies_to_dry_air_basis = ccpp_scheme_ref(
    "convert_dry_constituent_tendencies_to_dry_air_basis")
gravity_wave_drag_common_diagnostics = ccpp_scheme_ref("gravity_wave_drag_common_diagnostics")


@ccpp_composite
def gravity_wave_drag(*, sources):
    """Drag from unresolved gravity waves, inside an energy and water budget
    check.

    Roles:
        sources: list of Parameterizations, one per gravity wave source,
                 run in the order given; each writes its own diagnostics.
    """
    gravity_wave_drag_common()          # initialise all sources, read the namelist

    def drag():
        gravity_wave_drag_prepare_profiles()
        gravity_wave_drag_top_taper()   # taper drag towards the model top
        for source in sources:
            source.run()
        convert_dry_constituent_tendencies_to_dry_air_basis()
        gravity_wave_drag_common_diagnostics()
        apply_tendency_of_eastward_wind()
        apply_tendency_of_northward_wind()
        apply_constituent_tendencies()
        apply_heating_rate()
        qneg()
        geopotential_temp()
        update_dry_static_energy()

    energy_checked(process=drag)

"""Category: cloud macrophysics and microphysics."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import run_setup

cloud_fraction_options = ccpp_scheme_ref("cloud_fraction_options")


@ccpp_composite
def clouds(*, macro_micro=None, macrophysics=None, aerosol_activation=None,
           microphysics=None):
    """Clouds: condensation, cloud fraction, cloud water and ice, and the
    cloud properties radiation uses.

    Roles (each a Parameterization, or None):
        macro_micro:        a combined macro- and microphysics package whose
                            steps interleave (CAM4's Rasch-Kristjansson); it
                            owns its own energy check.
    or, as separate parameterizations, run in this order:
        macrophysics:       condensation and cloud fraction
        aerosol_activation: droplet activation and ice nucleation from
                            aerosols, for the microphysics
        microphysics:       cloud water, ice and precipitation
    Each role's setup (init only) runs just before it.

    The cloud fraction namelist is read first (init only): every cloud
    scheme uses the cloud fraction settings.
    """
    cloud_fraction_options()
    for role in (macro_micro, macrophysics, aerosol_activation, microphysics):
        if role is not None:
            run_setup(role)
            role.run()

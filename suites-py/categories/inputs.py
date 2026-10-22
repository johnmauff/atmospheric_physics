"""Category: prescribed inputs read from file at the start of a step."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

to_be_ccppized_temporary               = ccpp_scheme_ref("to_be_ccppized_temporary")
prescribe_radiative_gas_concentrations = ccpp_scheme_ref("prescribe_radiative_gas_concentrations")
rebin_seasalt                          = ccpp_scheme_ref("rebin_seasalt")


@ccpp_composite
def prescribed_inputs(*, ozone=None, aerosols=None, volcanic_aerosols=None):
    """Constituents and radiative gases, then the prescribed fields a model
    uses instead of interactive chemistry.

    Roles (each a Parameterization, or None when not prescribed):
        ozone:             prescribed ozone
        aerosols:          prescribed bulk aerosols
        volcanic_aerosols: prescribed volcanic aerosol
    Prescribed aerosols come with 4 sea-salt bins; they are rebinned to the
    2 the rest of the physics uses, after all prescribed fields are read.
    """
    # Helper schemes: all possible constituents from file or namelist.
    to_be_ccppized_temporary()
    prescribe_radiative_gas_concentrations()
    for role in (ozone, aerosols, volcanic_aerosols):
        if role is not None:
            role.run()
    if aerosols is not None:
        rebin_seasalt()

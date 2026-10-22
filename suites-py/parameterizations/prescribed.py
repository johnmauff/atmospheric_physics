"""Prescribed fields read from file (no interactive chemistry)."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

prescribed_ozone_scheme            = ccpp_scheme_ref("prescribed_ozone")
prescribed_aerosols_scheme         = ccpp_scheme_ref("prescribed_aerosols")
prescribed_aerosol_deposition_flux = ccpp_scheme_ref("prescribed_aerosol_deposition_flux")
tropopause_find                    = ccpp_scheme_ref("tropopause_find")
prescribed_volcanic_aerosol        = ccpp_scheme_ref("prescribed_volcanic_aerosol")


@ccpp_composite
def prescribed_bulk_aerosols_with_deposition():
    """Bulk aerosols and their deposition flux to the surface."""
    prescribed_aerosols_scheme()
    prescribed_aerosol_deposition_flux()


@ccpp_composite
def prescribed_volcanic():
    """Volcanic aerosol, placed relative to the tropopause."""
    tropopause_find()
    prescribed_volcanic_aerosol()


prescribed_bulk_aerosols = Parameterization(name="prescribed bulk aerosols",
                                            run=prescribed_aerosols_scheme)
prescribed_ozone = Parameterization(name="prescribed ozone", run=prescribed_ozone_scheme)
prescribed_aerosols_and_deposition = Parameterization(
    name="prescribed bulk aerosols, with deposition flux",
    run=prescribed_bulk_aerosols_with_deposition)
prescribed_volcanic_aerosols = Parameterization(name="prescribed volcanic aerosol",
                                                run=prescribed_volcanic)

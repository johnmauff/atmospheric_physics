"""MUSICA chemistry (MICM solver; TUV-x photolysis)."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

calc_dry_air_ideal_gas_density = ccpp_scheme_ref("calc_dry_air_ideal_gas_density")
musica_ccpp                    = ccpp_scheme_ref("musica_ccpp")


@ccpp_composite
def musica_chemistry():
    calc_dry_air_ideal_gas_density()        # air density MUSICA needs
    musica_ccpp()


musica = Parameterization(name="MUSICA chemistry", run=musica_chemistry)

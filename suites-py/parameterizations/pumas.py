"""PUMAS two-moment cloud microphysics (CAM5, CAM7), run on its own grid."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization
from categories.state import apply_constituent_tendencies, qneg, update_temperature

pumas_pre_main                   = ccpp_scheme_ref("pumas_pre_main")
micro_pumas_ccpp_dimensions_pre  = ccpp_scheme_ref("micro_pumas_ccpp_dimensions_pre")
micro_pumas_ccpp                 = ccpp_scheme_ref("micro_pumas_ccpp")
micro_pumas_ccpp_dimensions_post = ccpp_scheme_ref("micro_pumas_ccpp_dimensions_post")
clamp_number_concentrations      = ccpp_scheme_ref("clamp_number_concentrations")
pumas_post_main                  = ccpp_scheme_ref("pumas_post_main")
micro_pumas_ccpp_optics_limiter  = ccpp_scheme_ref("micro_pumas_ccpp_optics_limiter")


@ccpp_composite
def pumas_microphysics():
    pumas_pre_main()                        # PUMAS namelist (init only)
    micro_pumas_ccpp_dimensions_pre()       # host grid -> PUMAS grid
    micro_pumas_ccpp()
    micro_pumas_ccpp_dimensions_post()      # PUMAS grid -> host grid
    apply_constituent_tendencies()
    clamp_number_concentrations()
    qneg()
    # Size distributions for radiation from the updated constituents; needs
    # the air temperature before the heating is applied.
    pumas_post_main()
    micro_pumas_ccpp_optics_limiter()
    update_temperature()


pumas = Parameterization(name="PUMAS microphysics", run=pumas_microphysics)

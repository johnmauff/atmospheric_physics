"""Kessler warm-rain microphysics as a package (a composite scheme).

Kessler works in potential temperature and dry mixing ratios. The package
converts temperature to potential temperature and back itself; the
conversions between wet and dry mixing ratios are the host's (its
constituents, its convention), so they are hooks.

Potential temperature, the Exner function and the dry air density are the
package's own state: nothing outside it uses them. On a GPU they stay
resident on the device between the package's segments (composite-scheme 2g).
"""

from xdsl_ccpp.frontend.py_api import ccpp_package, ccpp_scheme_ref

calc_exner                     = ccpp_scheme_ref("calc_exner")
temp_to_potential_temp         = ccpp_scheme_ref("temp_to_potential_temp")
calc_dry_air_ideal_gas_density = ccpp_scheme_ref("calc_dry_air_ideal_gas_density")
kessler_scheme                 = ccpp_scheme_ref("kessler")
potential_temp_to_temp         = ccpp_scheme_ref("potential_temp_to_temp")
kessler_update                 = ccpp_scheme_ref("kessler_update")

kessler = ccpp_package("kessler")


@kessler.entry_point
def kessler_microphysics(*, to_dry, to_wet):
    """HOOKS: to_dry converts the water species Kessler uses to dry mixing
    ratios, to_wet converts them back."""
    calc_exner()
    temp_to_potential_temp()
    calc_dry_air_ideal_gas_density()
    to_dry()
    kessler_scheme()
    potential_temp_to_temp()
    to_wet()
    kessler_update()

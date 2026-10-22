"""Category: vertical diffusion (turbulent mixing in the boundary layer and
above).

The framework is the same for every model: prepare inputs, get diffusion
coefficients from the boundary-layer scheme (the role), compute provisional
updates for momentum, dry static energy and tracers, turn them into
tendencies, apply them.
"""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import run_setup
from categories.state import (
    apply_constituent_tendencies,
    apply_heating_rate,
    apply_tendency_of_eastward_wind,
    apply_tendency_of_northward_wind,
    geopotential_temp,
    qneg,
    update_dry_static_energy,
)

vertical_diffusion_options                   = ccpp_scheme_ref("vertical_diffusion_options")
zero_upper_boundary_condition                = ccpp_scheme_ref("zero_upper_boundary_condition")
hb_diff_set_vertical_diffusion_top           = ccpp_scheme_ref("hb_diff_set_vertical_diffusion_top")
tms_beljaars_zero_stub                       = ccpp_scheme_ref("tms_beljaars_zero_stub")
vertical_diffusion_prepare_inputs            = ccpp_scheme_ref("vertical_diffusion_prepare_inputs")
hb_diff_prepare_vertical_diffusion_inputs    = ccpp_scheme_ref("hb_diff_prepare_vertical_diffusion_inputs")
vertical_diffusion_not_use_rairv             = ccpp_scheme_ref("vertical_diffusion_not_use_rairv")
vertical_diffusion_set_temperature_at_toa_default = ccpp_scheme_ref(
    "vertical_diffusion_set_temperature_at_toa_default")
vertical_diffusion_interpolate_to_interfaces = ccpp_scheme_ref("vertical_diffusion_interpolate_to_interfaces")
vertical_diffusion_set_total_surface_stress  = ccpp_scheme_ref("vertical_diffusion_set_total_surface_stress")
vertical_diffusion_sponge_layer              = ccpp_scheme_ref("vertical_diffusion_sponge_layer")
implicit_surface_stress_add_drag_coefficient = ccpp_scheme_ref("implicit_surface_stress_add_drag_coefficient")
vertical_diffusion_wind_damping_rate         = ccpp_scheme_ref("vertical_diffusion_wind_damping_rate")
vertical_diffusion_diffuse_horizontal_momentum = ccpp_scheme_ref(
    "vertical_diffusion_diffuse_horizontal_momentum")
vertical_diffusion_set_dry_static_energy_at_toa_zero = ccpp_scheme_ref(
    "vertical_diffusion_set_dry_static_energy_at_toa_zero")
vertical_diffusion_diffuse_dry_static_energy = ccpp_scheme_ref("vertical_diffusion_diffuse_dry_static_energy")
vertical_diffusion_diffuse_tracers           = ccpp_scheme_ref("vertical_diffusion_diffuse_tracers")
vertical_diffusion_tendencies                = ccpp_scheme_ref("vertical_diffusion_tendencies")
vertical_diffusion_tendencies_diagnostics    = ccpp_scheme_ref("vertical_diffusion_tendencies_diagnostics")


@ccpp_composite
def vertical_diffusion(*, boundary_layer, mountain_stress=None, form_drag_damping=None,
                       undiffused_tracer_fluxes=None):
    """Vertical diffusion.

    Roles (Parameterizations, or None when not used):
        boundary_layer:  the boundary-layer scheme (required): ``setup``
                         reads its namelist, ``run`` computes the diffusion
                         coefficients, ``diagnostics`` (if any) writes them
                         out once the sponge layer has modified them.
        mountain_stress: turbulent mountain stress -- adds to the surface
                         drag coefficient.
        form_drag_damping: orographic form drag (Beljaars) -- adds to the
                         wind damping rate.
        undiffused_tracer_fluxes: applies the surface fluxes of tracers
                         that are not diffused (they would be lost
                         otherwise), after the tracer diffusion.
    Without mountain stress and form drag, their fields are zeroed.
    """
    run_setup(boundary_layer)                   # its namelist, first
    vertical_diffusion_options()
    zero_upper_boundary_condition()             # stub: zero upper boundary condition
    hb_diff_set_vertical_diffusion_top()        # init only: top of vertical diffusion
    if mountain_stress is None and form_drag_damping is None:
        tms_beljaars_zero_stub()

    # Inputs
    vertical_diffusion_prepare_inputs()         # pressure coordinates, potential temperature
    hb_diff_prepare_vertical_diffusion_inputs() # surface fluxes from the coupler
    vertical_diffusion_not_use_rairv()          # not WACCM-X
    vertical_diffusion_set_temperature_at_toa_default()
    vertical_diffusion_interpolate_to_interfaces()  # t, rho, p at interfaces
    vertical_diffusion_set_total_surface_stress()

    # Diffusion coefficients
    boundary_layer.run()
    vertical_diffusion_sponge_layer()           # modifies the momentum coefficient
    if boundary_layer.diagnostics is not None:
        boundary_layer.diagnostics()

    # Surface drag
    implicit_surface_stress_add_drag_coefficient()
    if mountain_stress is not None:
        mountain_stress.run()
    vertical_diffusion_wind_damping_rate()
    if form_drag_damping is not None:
        form_drag_damping.run()

    # Provisional updates, then tendencies
    vertical_diffusion_diffuse_horizontal_momentum()   # dissipates kinetic energy
    vertical_diffusion_set_dry_static_energy_at_toa_zero()
    vertical_diffusion_diffuse_dry_static_energy()
    vertical_diffusion_diffuse_tracers()
    if undiffused_tracer_fluxes is not None:
        undiffused_tracer_fluxes.run()
    vertical_diffusion_tendencies()
    vertical_diffusion_tendencies_diagnostics()  # before the tendencies are applied

    # State update
    apply_tendency_of_northward_wind()
    apply_tendency_of_eastward_wind()
    apply_constituent_tendencies()
    apply_heating_rate()
    qneg()
    geopotential_temp()
    update_dry_static_energy()

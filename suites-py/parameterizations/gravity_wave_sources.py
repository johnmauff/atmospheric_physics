"""Gravity wave sources for the gravity wave drag category."""

from xdsl_ccpp.frontend.py_api import ccpp_composite, ccpp_scheme_ref

from categories.interface import Parameterization

gravity_wave_drag_orographic             = ccpp_scheme_ref("gravity_wave_drag_orographic")
gravity_wave_drag_orographic_diagnostics = ccpp_scheme_ref("gravity_wave_drag_orographic_diagnostics")


@ccpp_composite
def orographic_gravity_waves():
    """Isotropic orographic source (CAM4)."""
    gravity_wave_drag_orographic()
    gravity_wave_drag_orographic_diagnostics()


orographic = Parameterization(name="orographic gravity waves", run=orographic_gravity_waves)

# CAM5's suite runs the orographic source without its history output (likely
# an omission in suite_cam5.xml).
orographic_without_diagnostics = Parameterization(
    name="orographic gravity waves, no diagnostics", run=gravity_wave_drag_orographic)


# --- CAM7 sources ------------------------------------------------------------

compute_total_precipitation_rate              = ccpp_scheme_ref("compute_total_precipitation_rate")
gravity_wave_drag_moving_mountain             = ccpp_scheme_ref("gravity_wave_drag_moving_mountain")
gravity_wave_drag_moving_mountain_diagnostics = ccpp_scheme_ref("gravity_wave_drag_moving_mountain_diagnostics")
gravity_wave_drag_convection_deep             = ccpp_scheme_ref("gravity_wave_drag_convection_deep")
gravity_wave_drag_convection_deep_diagnostics = ccpp_scheme_ref("gravity_wave_drag_convection_deep_diagnostics")
gravity_wave_drag_frontogenesis               = ccpp_scheme_ref("gravity_wave_drag_frontogenesis")
gravity_wave_drag_frontogenesis_diagnostics   = ccpp_scheme_ref("gravity_wave_drag_frontogenesis_diagnostics")
gravity_wave_drag_ridge                       = ccpp_scheme_ref("gravity_wave_drag_ridge")
gravity_wave_drag_ridge_beta                  = ccpp_scheme_ref("gravity_wave_drag_ridge_beta")
gravity_wave_drag_ridge_beta_diagnostics      = ccpp_scheme_ref("gravity_wave_drag_ridge_beta_diagnostics")


@ccpp_composite
def moving_mountain_gravity_waves():
    """Moving mountain source (SE dycore only); needs the total
    precipitation rate."""
    compute_total_precipitation_rate()
    gravity_wave_drag_moving_mountain()
    gravity_wave_drag_moving_mountain_diagnostics()


@ccpp_composite
def deep_convection_gravity_waves():
    """Convective source from deep convection (Beres et al.)."""
    gravity_wave_drag_convection_deep()
    gravity_wave_drag_convection_deep_diagnostics()


@ccpp_composite
def frontal_gravity_waves():
    """Frontally generated source."""
    gravity_wave_drag_frontogenesis()
    gravity_wave_drag_frontogenesis_diagnostics()


@ccpp_composite
def ridge_gravity_waves():
    """Ridge (anisotropic orographic) source."""
    gravity_wave_drag_ridge()
    gravity_wave_drag_ridge_beta()
    gravity_wave_drag_ridge_beta_diagnostics()


moving_mountain = Parameterization(name="moving mountain gravity waves",
                                   run=moving_mountain_gravity_waves)
deep_convection = Parameterization(name="deep convective gravity waves",
                                   run=deep_convection_gravity_waves)
frontogenesis = Parameterization(name="frontal gravity waves", run=frontal_gravity_waves)
ridge = Parameterization(name="ridge gravity waves", run=ridge_gravity_waves)

"""The interface between physics categories and parameterizations.

A category (convection, clouds, radiation, ...) has *roles* -- deep and
shallow convection, say.  A suite fills each role with a ``Parameterization``:
the object a parameterization module exports.  The category decides when
each part of it runs (and wraps it in energy checks, state updates, ...);
the parameterization decides which CCPP schemes that part calls.
"""

from dataclasses import dataclass
from typing import Callable, Optional


def nothing(*_args, **_kwargs):
    """An empty step: records no scheme call."""


@dataclass(frozen=True)
class Parameterization:
    """One parameterization, as a category uses it.

    Attributes:
        name: Short human-readable name (for documentation and errors).
        run: Records the parameterization's scheme calls -- usually a
            ``@ccpp_composite`` function, so it is a named block in the IR.
        setup: Init-only step (e.g. reading the parameterization's namelist),
            run by the category just before ``run``.
        diagnostics: History output the category places after a step of its
            own (e.g. after vertical diffusion's sponge layer).
        energy_flux_correction: Run after ``run`` and before the category's
            energy and water budget check, for a parameterization whose
            boundary fluxes need adjusting first.
    """

    name: str
    run: Callable[..., None]
    setup: Optional[Callable[[], None]] = None
    diagnostics: Optional[Callable[[], None]] = None
    energy_flux_correction: Optional[Callable[[], None]] = None


def run_setup(parameterization: Optional[Parameterization]) -> None:
    """Run *parameterization*'s setup step, if it has one."""
    if parameterization is not None and parameterization.setup is not None:
        parameterization.setup()

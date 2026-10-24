# Python suites (`suites-py`)

Python versions of the suite definition files in `../suites/`, written with
xdsl-ccpp's Python frontend (`xdsl_ccpp.frontend.py_api`). Each one defines
the same physics as its XML suite; what changes is how the suite is put
together.

## Layers

| Layer | Directory | What it holds |
|---|---|---|
| Suites | `suite_<name>.py` | One line per physics **category**, with the **parameterizations** that fill its roles. No scheme names. |
| Categories | `categories/` | One composite per category (convection, clouds, radiation, vertical diffusion, ...). It fixes the order of the category's roles and the energy checks and state updates around them, once for every CAM version. |
| Parameterizations | `parameterizations/` | One module per parameterization (ZM, Hack, PUMAS, RRTMGP, ...), each exporting a `Parameterization` object (`categories/interface.py`): `run`, and optionally `setup`, `diagnostics`, `energy_flux_correction`. |

For example, `suite_cam4.py`:

```python
convection(dry_adjustment=dadadj, deep=zhang_mcfarlane, shallow=hack)
clouds(macro_micro=rasch_kristjansson)
radiation(scheme=rrtmgp, aerosol_optics=bulk_aerosol_optics)
```

Two parameterizations are **packages** (composite schemes with hooks):
`rrtmgp_package.py` (the RRTMGP radiation package, wrapped by CAM's adapter
in `rrtmgp.py`) and `kessler_package.py` (Kessler microphysics, wrapped by
`kessler.py`, which fills its wet/dry conversion hooks).

The Python test suites (`suite_kessler_test.py`, `suite_rrtmgp.py`) live
in `../test/test_suites-py/` and import from this directory.

## Differences from the XML suites

The flattened scheme order is the XML's, except that a few **init-only**
schemes (namelist readers and similar, with only register/init/final
phases) run inside their category instead of at the top of the group:
`zm_conv_options`, `cloud_fraction_options`,
`compute_cloud_fraction_two_moment`, `set_cloud_fraction_top` (cam4: the
first two; cam5: three; cam7: all four). They run once, before time
stepping, and do not read or write anything another init routine they
moved past uses, so results are unchanged -- confirmed bit-for-bit for cam4
(CAM-SIMA QPC4, all history files identical).

## Using them in CAM-SIMA

Set `CCPP_GENERATOR=xdsl_ccpp` and `CCPP_SUITE_DIR=suites-py` in the case
(test suites: `CCPP_SUITE_DIR=test/test_suites-py`). `cam_autogen.py` runs
the suite script with this directory on `PYTHONPATH`.

## Tools (xdsl-ccpp)

Run from `src/physics/ncar_ccpp` with the xdsl-ccpp Python environment.

- **Verify** that every Python suite matches its XML suite (flattened IR;
  with `--case`, also the caps generated with a built case's host
  metadata):

  ```
  python -m xdsl_ccpp.tools.verify_py_suites \
      --py-dir suites-py test/test_suites-py --xml-dir suites test/test_suites \
      --scheme-dir schemes test/test_schemes --import-path suites-py --legacy-mode \
      [--case CASE_DIR ...] [--strict]
  ```

  It prints PASS/FAIL per suite and names any accepted difference (moved
  init-only schemes, package-private names). `--strict` fails on the moves.

- **Convert** a Python suite to an ordinary suite definition file (the flat
  scheme list), e.g. for a generator that only reads XML:

  ```
  python -m xdsl_ccpp.tools.py2sdf suites-py/suite_cam4.py \
      --import-path suites-py --legacy-mode -o suite_cam4.xml
  ```

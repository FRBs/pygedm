# Dev environment

`pygedm` has compiled C/C++ extensions (`ymw16`, `ne21c`) that need a C++ compiler,
`pybind11`, and the `f2c` runtime library/headers. Use the `pixi.toml` env for this —
it's a plain conda environment (no PyPI resolution), which matters on this network:

- `pip install -e .` / `uv pip install -e .` fails here with `ModuleNotFoundError: No
  module named 'pybind11'` (not installed) or, once pybind11 is present,
  `fatal error: 'f2c.h' file not found` (the `ne21c` extension links against `f2c`,
  which isn't a PyPI package — it must come from conda-forge).
- A `pyproject.toml`-integrated pixi setup (`[tool.pixi.project]` alongside `[project]`,
  as used in `ska-ost-low-uv`) auto-treats `[project.dependencies]` as PyPI deps and
  needs to fetch prefix.dev's conda/PyPI name-mapping file. That fetch reliably fails
  on this network (`failed to download pypi name mapping: io error: unexpected end of
  file`) even though direct `curl` to the same URLs succeeds — looks like a proxy/DPI
  issue specific to pixi's request pattern, not a real outage. Don't waste time
  retrying it; use the standalone `pixi.toml` below instead, which never touches that
  endpoint (pure conda deps + a plain `pip install -e .` task, no `[tool.pixi.pypi-dependencies]`).

## Setup

```bash
pixi install -e test
pixi run -e test install   # pip install -e . --no-build-isolation, builds ymw16/ne21c
```

## Running tests

```bash
pixi run -e test pytest tests/ -q
```

If you need a one-off script or REPL, use `pixi run -e test python ...` — do not use
the system `python`/`pip`, they don't have `pybind11`/`f2c`.

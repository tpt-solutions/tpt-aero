# tpt-aero (Python)

Python bindings for [tpt-aero](https://github.com/tpt-solutions/tpt-aero). The
functions take the same TOML case text as the `tpt-aero` command and return
the same text, so a notebook and the command line agree.

## Build

```text
pip install maturin
cd crates/tpt-aero-py
maturin develop --release          # installs into the current environment
```

`maturin build --release` writes a wheel instead. The module is built against
the stable ABI (`abi3-py39`), so one wheel covers Python 3.9 and later.

## Use

```python
import tpt_aero

case = open("uav.toml").read()        # the same file `tpt-aero wing` reads
print(tpt_aero.wing_csv(case))        # alpha sweep, CSV
print(tpt_aero.wing_warnings(case))   # validity envelope, empty when inside it
svg = tpt_aero.wing_svg(case)         # lift-drag polar as SVG text

print(tpt_aero.isa_csv([0.0, 5000.0]))
rotor = tpt_aero.rotor(open("quadrotor.toml").read())   # dict of the report
```

Errors raise `ValueError` with the message the command line prints.

## Pandas recipe

The rotor RPM sweep is CSV, which pandas reads directly:

```python
import io
import pandas as pd
import tpt_aero

csv = tpt_aero.rotor_sweep_csv(open("quadrotor.toml").read(), 4000, 9000, 500)
df = pd.read_csv(io.StringIO(csv))
df["efficiency_proxy"] = df["thrust_n"] / df["power_w"]
print(df.head())
```

A single operating point loads straight from its dict:

```python
point = pd.DataFrame([tpt_aero.rotor(open("quadrotor.toml").read())])
```

The wing sweep also reads with `pd.read_csv(io.StringIO(tpt_aero.wing_csv(case)))`;
the `l_over_d` column is blank at zero drag, which pandas reads as `NaN`.

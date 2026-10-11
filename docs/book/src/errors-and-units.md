# Errors and units

## Errors

The CLI prints `error: <message>` to standard error and exits with status 1. A
usage error (an unknown flag, a missing argument) exits with status 2 and
prints the help.

Messages name the input that failed and, where it helps, the valid range:

- `cannot read <path>: <reason>` when a case file cannot be opened.
- `case file error: …` for TOML that does not parse or has an unknown key or
  an unknown surface kind.
- Value errors with no prefix, naming the field, for example
  `control rudder: need 0 <= inboard < outboard <= span/2 (0.75 m), got 0 to 1.5`.
- `altitude 90000 m: altitude out of ISA range (0..84852 m geopotential)`.
- `solver error: …` when the solver cannot finish, for example a singular
  influence matrix or a trim that does not converge.

Warnings about the validity envelope (aspect ratio, taper, sweep, angle of
attack, Mach number) go to standard error with the prefix `warning:`. They do
not stop the run.

## Units

Inputs and outputs are SI unless a column says otherwise.

| quantity | unit |
|---|---|
| length, span, chord, altitude | m (ISA altitude is geopotential) |
| speed | m/s |
| angle | deg in case files and the CLI, rad in the library |
| density | kg/m³ |
| mass | kg |
| force | N |
| moment | N·m |
| rotor speed | rpm in case files, rad/s inside the solver |
| power | W |
| coefficients | dimensionless: C_L, C_D, C_m, C_l, C_n |

`isa --imperial` switches the table to feet, °R, psf, slug/ft³ and ft/s. The
input altitudes are then feet.

Angles in the `derivatives` output are per radian, except the non-dimensional
rate derivatives, which use `q·c/2V` and `p·b/2V`.

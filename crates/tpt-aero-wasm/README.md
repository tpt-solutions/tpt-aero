# tpt-aero-wasm

WebAssembly bindings for tpt-aero, with a small browser demo in `www/`. The
functions take the same TOML case text as the `tpt-aero` command and return the
same text, so the browser and the command line give the same numbers.

## Build the demo

```text
rustup target add wasm32-unknown-unknown
cargo install wasm-bindgen-cli --version 0.2.129   # match the wasm-bindgen crate version
cargo build -p tpt-aero-wasm --target wasm32-unknown-unknown --release
wasm-bindgen --target web --out-dir crates/tpt-aero-wasm/www/pkg \
  target/wasm32-unknown-unknown/release/tpt_aero_wasm.wasm
python -m http.server --directory crates/tpt-aero-wasm/www 8000
```

Open <http://localhost:8000>. The page loads `pkg/tpt_aero_wasm.js`, which
`wasm-bindgen` writes; the `pkg/` folder is not committed.

## Exported functions

| function | returns |
|---|---|
| `wing_csv(case)` | the alpha sweep as CSV |
| `wing_svg(case)` | the lift-drag polar as an SVG document |
| `wing_warnings(case)` | validity warnings, one per line |
| `trim_text(case, speed, altitude)` | the level trim report |
| `report_markdown(case, speed, altitude)` | the Markdown stability report |
| `rotor_text(case)` | one BEMT operating point |
| `isa_table(altitudes)` | the standard atmosphere, CSV |
| `wing_template()` | a starting UAV wing case |

Errors are thrown to JavaScript as strings with the same messages the command
line prints.

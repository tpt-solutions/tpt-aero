# Installing and building

tpt-aero builds on stable Rust (see `rust-toolchain.toml` in the repository).
It depends on the sibling `tpt-math` workspace by path, so check out both
side by side:

```text
parent/
├── tpt-aero/
└── tpt-math/
```

Then, from `tpt-aero/`:

```sh
cargo test --workspace             # the full test suite
cargo run -p tpt-aero-cli -- --help
cargo run --release -p tpt-aero --example propeller_flight
```

Install the command with `cargo install --path crates/tpt-aero-cli`.

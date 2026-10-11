# Contributing to tpt-aero

Thanks for looking at tpt-aero. It is a pure-Rust aerodynamics library, and
most of the value is in the numbers being right, so a change that alters a
result needs evidence.

## Before you start

- Read [docs/CONVENTIONS.md](docs/CONVENTIONS.md) for frames, signs and units.
- Check [todo.md](todo.md) for open items. Work that is listed
  there is welcome; say in the issue which item you are taking.
- Open an issue first for anything that changes an API or a published number.

## Making a change

1. `cargo test --workspace` must pass, including doctests.
2. `cargo clippy --workspace --all-targets -- -D warnings` must be clean.
3. `cargo fmt --all` for the files you touched. (The repository has one
   pre-existing formatting diff in `crates/tpt-aero-rotor/src/lib.rs`; leave it
   alone unless you are already editing that file.)
4. A bug fix lands with a regression test that fails before the fix.
5. A change to a validated result updates [docs/VALIDATION.md](docs/VALIDATION.md)
   in the same pull request.

## What the tests are for

- Unit tests sit next to the code.
- `tests/` holds integration and validation tests against analytical limits,
  AVL and published data.
- `tpt-aero-verify` holds property tests and the Kani harnesses (`cfg(kani)`).

## Licence

Contributions are dual-licensed MIT OR Apache-2.0, the same as the project.
Do not add a dependency that is Apache-2.0 only (see `deny.toml`).

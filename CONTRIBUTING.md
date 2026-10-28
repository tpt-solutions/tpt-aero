# Contributing to tpt-aero

Thanks for looking at tpt-aero. It is a pure-Rust aerodynamics library, and
most of the value is in the numbers being right.

Contributions are made through issues only. Pull requests are not accepted.
Open an issue to report a bug, ask a question or propose a change, and include
enough evidence for any change that alters a result.

## Before you open an issue

- Read [docs/CONVENTIONS.md](docs/CONVENTIONS.md) for frames, signs and units.
- Check [todo.md](todo.md) for open items and mention the item your issue
  relates to.
- Say in the issue if it changes an API or a published number, and what the
  expected value is and where it comes from.

## What a useful bug report includes

- The version or commit you are using.
- A minimal reproduction, such as a short code snippet or test case.
- The actual result, the expected result and how the expectation was derived.

## What the tests are for

- Unit tests sit next to the code.
- `tests/` holds integration and validation tests against analytical limits,
  AVL and published data.
- `tpt-aero-verify` holds property tests and the Kani harnesses (`cfg(kani)`).

## Licence

Contributions are dual-licensed MIT OR Apache-2.0, the same as the project.
If an issue proposes a new dependency, state its licence; it must be allowed
by `deny.toml`.

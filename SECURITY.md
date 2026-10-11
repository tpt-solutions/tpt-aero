# Security policy

tpt-aero is a numerical library and a command-line tool. It has no network
access and no `unsafe` code (`unsafe_code = "forbid"` in the workspace lints).
The realistic security issues are a crafted case file that crashes the CLI or
exhausts memory, and a dependency advisory.

To report a vulnerability, use GitHub's private advisory form:
<https://github.com/tpt-solutions/tpt-aero/security/advisories/new>. Please do
not open a public issue for it. We will acknowledge the report, agree a fix
and disclosure date with you, and credit you if you wish.

Dependency advisories are checked by `cargo deny`; the CI workflow runs it.

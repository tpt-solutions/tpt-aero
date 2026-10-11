# AVL reference decks

Input decks for the AVL 3.52 comparison in `crates/tpt-aero-stability/tests/avl_reference.rs`.
Results are compared with `docs/VALIDATION.md`.

| Deck | Geometry | Used for |
|---|---|---|
| `wing_rectangular.avl` | Rectangular wing, span 10 m, chord 1.2 m, flat, 20 × 4 panels | C_Lα, C_Di |
| `wing_dihedral.avl` | As above with 5° dihedral | C_Yβ, C_lβ, C_nβ |
| `t_tail_fine.avl` | Dihedral wing (40 strips), T-tail at z = 1.8 m (20), fin 0–1.5 m (24) | Full set, refined grid |
| `t_tail_coarse.avl` | As `t_tail_fine` at 20 / 10 / 6 strips | Convergence check only |

`derivatives.in` is the AVL command script: `oper`, `a a 2` (alpha 2°), `x`
(execute), `st` (stability derivatives), `quit`.

## Running

```
avl t_tail_fine.avl < derivatives.in
```

The derivatives are per radian for α and β, and per unit `pb/2V`, `rb/2V` for the
rates, matching `tpt-aero-stability`.

## Building AVL from source (Windows, MinGW gfortran)

The AVL 3.52 source is at <https://web.mit.edu/drela/Public/web/avl/avl3.52.tgz>.
From the extracted tree:

```
cd plotlib && make mingwDP        # libPlt_gDP.a
cd ../eispack && make -f Makefile.gfortran
cd ../bin && make -f Makefile.mingwDP   # avl.exe
```

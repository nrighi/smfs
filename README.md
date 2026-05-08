# SMFs

`SMFs` is a Wolfram Language / Mathematica package for working with truncated genus-2 Siegel modular forms.

It provides two complementary ways to work:

- **directly in nome variables** through truncated Fourier expansions,
- **symbolically in geometric variables** `(T, Z, U)` through replacement rules.

The package includes:

- truncated genus-2 Eisenstein series via `Ek`,
- truncated expansions of the Igusa cusp forms via `Chi10` and `Chi12`,
- symbolic modular-form heads `E4`, `E6`, `E10`, `E12`, `\[Chi]10`, `\[Chi]12`,
- precomputed replacement rules via `PrecomputedRules[n]`,
- support for derivatives of the symbolic modular forms.

---

## Files

- `SMFs.m` — the package
- `examples.nb` — example notebook illustrating the main workflows

---

## Loading the package

In a notebook located in the same directory as `SMFs.m`, load the package with

```mathematica
SetDirectory[NotebookDirectory[]];
Get["SMFs.m"];
```

To inspect the exported symbols, use

```mathematica
Names["SMFs`*"]
```

---

## Main workflows

### 1. Direct expansion in nome variables

The functions

- `Ek[k, n, q1, r, q2]`
- `Chi10[n, q1, r, q2]`
- `Chi12[n, q1, r, q2]`

return truncated Fourier expansions in the nome variables

- `q1 = Exp[2 Pi I T]`
- `r  = Exp[2 Pi I Z]`
- `q2 = Exp[2 Pi I U]`

Examples:

```mathematica
Chi10[2, q1, r, q2]
Chi12[2, q1, r, q2]
Ek[4, 2, q1, r, q2]
```

---

### 2. Symbolic replacement in `(T, Z, U)`

The symbolic heads

- `E4[T, Z, U]`
- `E6[T, Z, U]`
- `E10[T, Z, U]`
- `E12[T, Z, U]`
- `\[Chi]10[T, Z, U]`
- `\[Chi]12[T, Z, U]`

can be replaced by truncated Fourier expansions using

```mathematica
PrecomputedRules[n]
```

For example:

```mathematica
rules = PrecomputedRules[2];
\[Chi]10[T, Z, U] /. rules
```

or for a composite expression:

```mathematica
expr = E4[T, Z, U]^3/\[Chi]12[T, Z, U];
expr /. PrecomputedRules[2]
```

This is especially useful for numerical evaluation, plotting, and symbolic experimentation.

---

## Derivatives

The precomputed rules also support derivatives with respect to `T`, `Z`, and `U`.

Examples:

```mathematica
Derivative[1, 0, 0][E4][T, Z, U] /. PrecomputedRules[2]
Derivative[0, 1, 0][\[Chi]12][T, Z, U] /. PrecomputedRules[2]
Derivative[1, 1, 0][\[Chi]10][T, Z, U] /. PrecomputedRules[2]
```

---

## Fixed points of the Siegel upper half space

The example notebook also includes a collection of distinguished fixed points of the Siegel upper half space, stored as replacement rules

- `pt1`, `pt2`, ..., `pt6`

and collected into the list

- `allpts`

These can be used to evaluate modular forms and their derivatives numerically at special points.


---

## Exported symbols

### Direct nome expansions

- `Ek`
- `Chi10`
- `Chi12`

### Symbolic modular forms

- `E4`
- `E6`
- `E10`
- `E12`
- `\[Chi]10`
- `\[Chi]12`

### Replacement rules

- `PrecomputedRules`

---

## Notes

- All expansions are **truncated**, so results depend on the chosen order `n`.
- The symbolic variable order is consistently **`(T, Z, U)`**.
- `Chi10` and `Chi12` are the **nome-series generators**.
- `\[Chi]10` and `\[Chi]12` are **symbolic modular-form heads** intended to be replaced using `PrecomputedRules[n]`.

---

## Example notebook

The file `examples.nb` provides a minimal introduction to the package and demonstrates:

- how to load `SMFs.m`,
- how to inspect the exported symbols,
- how to compute direct truncated expansions in nome variables,
- how to replace symbolic modular forms in `(T, Z, U)`,
- how to evaluate derivatives,
- how to test expressions at fixed points of the Siegel upper half space.

---

## Citation

If you use `SMFs` in work that leads to a publication, please cite the companion paper:

```bibtex
@article{Leedom:2026nby,
    author = "Leedom, Jacob M. and Righi, Nicole and Westphal, Alexander",
    title = "{Automorphic Structures of Heterotic Vacua}",
    eprint = "2605.05322",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    reportNumber = "DESY-26-063",
    month = "5",
    year = "2026"
}
```

---

## Author

Nicole Righi

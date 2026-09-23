# Julia model library

This is the repository's only model-implementation layer. It contains
calcium–phosphate homeostasis models and direct PTH/bone-resorption mechanisms
that support the integrated calcium models.

Included models:

- `KrollParams` / `simulate_kroll` — PTH timing and bone-cell response. It
  does not expose a plasma-calcium state; it is retained as a direct
  PTH-to-bone-resorption mechanism.
- `LemaireParams` / `simulate_lemaire` — PTH/vitamin-D-regulated osteoblast and
  osteoclast dynamics used in the integrated-model lineage.
- `RaposoParams` / `simulate_raposo` — 11-state systemic calcium/phosphate
  homeostasis model.
- `PetersonKhuranaParams` / `simulate_peterson_khurana` — 32-state integrated
  calcium, phosphate, endocrine, kidney, bone, and dose-event variant.
- `simulate_peterson_khurana_bmd` — BMD observation extension of the integrated
  model.
- `GranjonParams` / `simulate_granjon` — 22-state calcium/phosphate/FGF23
  model translated from the public CaPO4Sim source artifact.

Run the test suite from the repository root:

```sh
julia --project=julia -e 'using Pkg; Pkg.test()'
```

The library keeps source-faithful behavior separate from physiological repairs
or new mechanisms. Each addition needs traceable equations, parameters, and
automated physical/numerical tests. See `../LITERATURE_REVIEW.md` for the
scope rule and primary sources.

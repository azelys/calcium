# Julia model library

This directory is the reproducible Julia implementation layer for the
calcium-homeostasis and bone-remodeling project.

Current models:

- `KomarovaParams` and `simulate_komarova`: Komarova et al. (2003) power-law
  osteoblast/osteoclast model with the published bone-mass balance.
- `LemaireParams` and `simulate_lemaire`: Lemaire et al. (2004) mass-action
  osteoblast/osteoclast model, with the public CellML extension represented by
  explicit optional osteoblast removal parameters.
- `PivonkaParams` and `simulate_pivonka`: Pivonka et al. (2008) RANKL/OPG/
  TGF-beta signaling model with algebraic PTH, OPG, and RANKL observables.
- `KrollParams` and `simulate_kroll`: Kroll (2000) PTH pulse-schedule model
  for preosteoblasts, osteoblasts, osteoclasts, and IL-6.
- `RaposoParams` and `simulate_raposo`: Raposo et al. (2002) 11-state minimal
  calcium/phosphate homeostasis model, translated from the public Physiome
  CellML artifact.
- `PetersonKhuranaParams` and `simulate_peterson_khurana`: faithful Julia port
  of the local 32-state R model, including its default daily dosing events.
- `simulate_peterson_khurana_bmd`: the 33-state extension from Peterson & Riggs
  (2012), which adds the indirect nonlinear lumbar-spine BMD equation driven by
  the model's BSAP and CTx proxies.
- `GranjonParams` and `simulate_granjon`: 22-state calcium/phosphate/FGF23
  model translated from the public CaPO4Sim compiled core. It includes rapid
  and slow bone mineral pools, intracellular phosphate, Ca/P chemical
  complexes, renal handling, PTH, calcitriol, and FGF23. Time is in minutes;
  the public artifact's disease-state initial conditions are available through
  `granjon_initial_state(:hypopara)`, `:hypoD3`, `:hyperD3`, and `:php1`.
- `NelsonParams` and `simulate_nelson`: seven-state surgical-menopause model
  translated from the public MATLAB repository. It imposes natural or
  surgical estrogen trajectories, tracks pre-/mature osteoblast and
  osteoclast populations, osteocytes, sclerostin, and normalized bone density,
  and can activate the paper's post-surgical apoptosis/differentiation effects.

The Granjon source artifact is the [CaPO4Sim repository](https://github.com/DivadNojnarg/CaPO4Sim),
especially its `compiled_core.c`, `calcium_phosphate_core.R`, and model-engine
initial-state CSV files. The Julia translation preserves the artifact's
as-coded conventions and keeps interventions such as `Bispho`, `Furo`, and
`Cinacal` as parameter overrides.

The Nelson source artifact is the [SurgicalMenopauseBone repository](https://github.com/ashleefv/SurgicalMenopauseBone),
which supplies the MATLAB equations, parameter file, steady-state
initialization, estrogen forcing, and spine-BMD data used by the paper. The
Julia port preserves the source's seven-state ODE and its distinction between
the baseline model and the optional post-surgical effects.

Run the tests from this directory's parent with:

```sh
julia --project=julia -e 'using Pkg; Pkg.test()'
```

The implementation policy is to preserve an as-published/as-coded mode and
keep corrections or extensions explicit. Every future model should add its
equations, provenance, baseline simulation, and physical/numerical tests
together.

The local integrated model's full default run is intentionally not part of the
fast package test suite. A 210-day run with daily PTH/vitamin-D events should
reproduce the R baseline: 10,081 rows, final plasma-calcium amount
21.9080333, final intracellular-bone calcium amount 5.20782635, and PTH range
1.164587–679.338954. With the 2012 BMD extension, the same run gives a
baseline-normalized lumbar-spine BMD of approximately 104.65671 at day 210
under the source example dosing schedule.

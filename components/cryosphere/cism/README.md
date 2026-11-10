# CISM stage 2 — basal traction, effective pressure, flow factor, cap coupling boundary, temperature step, SMB application, and the DIVA velocity-solve pointwise laws

Stage-2 implementation of the CISM ice-sheet component as built in CESM
(`components/cism` = cismwrap_3_00_002, source `source_cism` =
cism_main_3.00.002, cesm3.0-alphabranch checkout `stage1/cesm`, reference run
T1850Gg/f09_g17_gris20) as esm 1.3.0 files: the pointwise laws of the
reference run's basal/thermodynamic/transport/velocity physics that
CISM_LANE.md §2 marks **same law** as EqIce — the pseudo-plastic basal drag
(row 8), the Bueler–van Pelt effective pressure (row 9), and the Paterson–Budd
flow factor (row 4) — plus the CESM-only **cap coupling boundary** (row 15,
missing in EqIce), the **temperature step** (row 5, the prognostic ice
temperature solve on CISM's σ grid, `cism_thermo.esm` +
`cism_thermo_discretization.esm`), the **SMB application / tracer surface-bed
handling / post-transport cleanup** slice of the thickness/tracer solve
(`cism_transport.esm`, rows 6/7/14), and the **DIVA velocity-solve pointwise
laws** (`cism_velo.esm` + `cism_velo_discretization.esm`, row 2).

The `stage1/...` paths and `bugs.md` citations throughout refer to the EqESM
project repository (the CESM checkout instrumented to produce the dump tuples
the tests pin, and the project's centralized bug list).

The shared math that equals an EqIce law is **imported by reference from the
EqIce stage-2 libraries**; because the EqIce project is not yet upstream in
EarthSciModels (no EqIce branch or PR exists as of 2026-11-10), this PR ships
clearly-marked **temporary copies** of exactly the imported laws
(`eqice_basal_templates.esm`, `eqice_ice_rheology_templates.esm`, bodies
byte-identical to `eqice/stage2/`), which the CISM files import by reference.
When the EqIce components land in EarthSciModels, a follow-up PR re-points the
imports at the EqIce files and deletes the copies. The CISM-specific pieces —
the bed-elevation φ ramp, the bwat=0 → N=overburden branch, the min(·, P) cap
and floating zeroing, the T* = clamp(T−T_pmp, −50, 0) + 273.15 construction,
the cap tsfc −tkfrz conversion with the −250 °C sentinel clamp, and the CISM
constants (arrmlh/arrmll/actenh/actenl/gascon, tkfrz, glide_types defaults) —
live in `cism_templates.esm`. The `tests` blocks pin the stage-1
instrumentation dumps (`stage1/dumps/cism/`) — see [Testing](#testing).

All model logic lives in `.esm` files; this README is the human-readable map.
Sources: `stage1/cesm/components/cism/source_cism/libglissade/`
(`glissade_basal_traction.F90`, `glissade_therm.F90`,
`glissade_grid_operators.F90`, `glissade_velo_higher.F90`,
`glissade_mass_balance.F90`, `glissade.F90`, `glissade_transport.F90`,
`glissade_masks.F90`), `source_cism/libglimmer/glimmer_physcon.F90`,
`source_cism/libglimmer/glimmer_utils.F90`,
`source_cism/libglide/glide_types.F90`, `source_cism/libglide/glide_setup.F90`,
`source_cism/libglad/glad_main.F90`, `source_cism/libglad/glad_timestep.F90`,
`source_cism/libglad/glad_output_states.F90`,
`source_cism/libglad/glad_output_fluxes.F90`,
`source_cism/drivers/cpl/nuopc/glc_import_export.F90`, the run namelists
(`stage1/cime/cases/T1850Gg_cism_inst/run/cism.gris.config`), and the CISM→EqIce
mapping of `stage1/runs/CISM_LANE.md` §2.

Reference configuration pinned (from `cism.gris.config` of the T1850Gg run and
the source defaults): `which_ho_babc = 3` (HO_BABC_PSEUDO_PLASTIC_OLD),
`which_ho_effecpress = 4` (HO_EFFECPRESS_BWAT_BVP), `flow_law = 2`
(FLWA_PATERSON_BUDD), `which_ho_flow_enhancement_factor = 0` (E = 1),
`which_ho_approx = HO_APPROX_DIVA = 4`, `which_ho_gradient = 0`, `margin = 1`,
`assemble_taud = 1`, `diva_slope_correction = .true.`, `max_slope = 0.5`,
`thck_gradient_ramp = 50 m`, `thklim = 1 m` (velo) / `thklim = 0` (SMB masks),
`upn = 11` (10 layers), `nquad = 4`, `which_smb = SMB_APPLY_RATE`,
`which_ho_vertical_remap = HO_VERTICAL_REMAP_FIRST_ORDER`,
`which_ho_ice_age = HO_ICE_AGE_COMPUTE`, `tstep_mbal = 8760 h` (1 yr),
`dt = 0.5 yr`, `eus = 0`, `pseudo_plastic_q = 0.5`,
`pseudo_plastic_u0 = 100 m/yr`, the source defaults
`pseudo_plastic_phimin/phimax = 5/40 deg`, `pseudo_plastic_bedmin/bedmax =
−300/700 m` (glide_types.F90:2380-2383, not exposed in the namelist),
`bwat_till_max = 2 m`, `N_0 = 1000 Pa`, `e_0 = 0.69`, `C_c = 0.12`,
`effecpress_delta = 0.02` (glide_types.F90:2297-2330), `rhoi = 917 kg/m³`,
`rhoo = 1026 kg/m³` (SHR_CONST_RHOSW), `rhow = 1000 kg/m³`, `grav = 9.80616
m/s²`, `pmlt = 9.7456e-8 K/Pa`, `scyr = 31536000 s`, and the flow-law
constants `arrmlh = 1.733e3`, `arrmll = 3.613e-13` (Pa⁻³ s⁻¹), `actenh =
139e3`, `actenl = 60e3` J/mol, `gascon = 8.314` J/mol/K
(glimmer_physcon.F90:89-104). The temperature-step constants (CCSMCOUPLED):
rhoi = 917 kg/m³, shci = 2117.27 J/kg/K, lhci = 333700 J/kg, coni = 2.1
W/m/K, pmp_threshold = 1e-3 degC, eps08 = 1e-8, eps11 = 1e-11, dttem =
1.5768e7 s (the 0.5-yr glissade step).

## File map

| File | Contents | Fortran source |
|---|---|---|
| `cism_templates.esm` | Expression-template library (17 match-less templates): `cism_phi_ramp` — the PSEUDO_PLASTIC_OLD bed-elevation friction-angle ramp (φ = phimin below bedmin, phimax above bedmax, linear ramp between); `cism_tempcor` — the Paterson–Budd pressure-corrected temperature clamp(temp − pmptemp, −50, 0) degC; `cism_tsfc_import` — the cap import −250 °C sentinel clamp tsfc = ifelse(tsfc_c < −250, 0, tsfc_c) on the consumer's converted intermediate tsfc_c = tsfc_raw − tkfrz; `cism_internal_melt_rate` — the layer melt-energy density rhoi·shci·max(0, T−Tpmp) (J/m³) of the TEMP_PROGNOSTIC cap; `cism_basal_melt_rate` — the basal flux-balance melt rate bflx_mlt/(lhci·rhoi) with the eps11 guard on the FLUX; `cism_basal_freeze_on` — the freeze-on predicate; `cism_basal_melt_cap` — the post-melt basal temperature cap; `cism_hybrid_edge_gradient` / `cism_vertex_gradient` / `cism_driving_stress` / `cism_theta_slope_x` / `cism_diva_slope_factor` / `cism_beta_eff` / `cism_stagger_masked_mean` — the DIVA surface-geometry, driving-stress, slope-correction and effective-friction stencil laws and the ice-masked vertex stagger; `cism_sigma` — the glide stretched-sigma coordinate sigma(x, n) = (1 − (x+1)^(−n))/(1 − 2^(−n)); `cism_ocean` — the shared ocean predicate ocean = topg < eus ∧ thck ≤ 0; `cism_tracer_mix` — the shared mixing law T′ = (H·T + s·T′)/(H + s). The Clapeyron pmp itself is imported from the temporary EqIce copies (`pressure_melting_temperature` of `overburden_pressure`, with T_melting = 0 degC and beta_Clausius_Clapeyron = pmlt) — no CISM pmp template | `glissade_basal_traction.F90`: 256-277; `glissade_therm.F90`: 2596-2598, 2127-2155, 1930-2012; `glc_import_export.F90`: 326-333; `glissade_velo_higher.F90`: 2174-2193, 2549-2564, 5756-5770, 595-630; `glissade_grid_operators.F90`: 999-1094, 1482+; `glide_setup.F90`: 303-308, 380-387; `glide_mask.inc`: 33-34; `glissade_masks.F90`: 70-96; `glissade_mass_balance.F90`: 1380-1386, 1458-1465 |
| `cism_cap.esm` | Models `GlcCapImport` (the cap import transformations: the named intermediate tsfc_c = tsfc_raw − tkfrz and the −250 °C sentinel clamp via the one-parameter `cism_tsfc_import`, qsmb pass-through) and `GlcCapExport` (the export laws: topo = usrf on the active grid, ice_covered = (thck > 0) ∧ (usrf > 0), and the coupling-interval averages rofi/rofl/hflx written as storage ODEs — the continuous limit of the Fortran's per-dynamics-step accumulation, see the cap section below) | `glc_import_export.F90`: 326-333; `glad_output_states.F90`: 48-158; `glad_output_fluxes.F90`: 77-172; `glad_main.F90`: 971-1027 |
| `cism_basal.esm` | Models `Effecpress` (cell-centered N: overburden P = ρᵢgH, the Bueler–van Pelt N0(δP/N0)^s·10^((e0/Cc)(1−s)) law with s = clamp(bwat/bwat_till_max, 0, 1) imported by reference from `eqice_basal_templates.esm` `mohr_coulomb_N_interior`, the bwat=0 → N = P branch, the min(·, P) cap and the floating-cell zero) and `BasalTraction` (φ ramp → tanφ, speed = √(u²+v²+0.01²), τ_c = tanφ·N_stag, β = τ_c/(u0^q·speed^(1−q)) — the drag law itself imported from `eqice_basal_templates.esm` `pseudo_plastic_drag`) | `glissade_basal_traction.F90`: calc_effecpress 668-1114 (BWAT_BVP 944-985, init 784, stagger 1096-1098, cap 1082-1086), calcbeta 69-664 (PSEUDO_PLASTIC_OLD 250-292, smallnum 136/189) |
| `cism_flowlaw.esm` | Model `FlowFactor` — flwa = E·a·exp(−Q/(RT*)) with T* = tempcor + 273.15, the (a, Q) pair switching at T* = 263.15 K; the softness law imported by reference from `eqice_ice_rheology_templates.esm` `pb_softness`, the CISM constants (arrmlh/arrmll/actenh/actenl/gascon) bound at the call site, the T* construction CISM-specific | `glissade_therm.F90`: 2408-2640 (FLWA_PATERSON_BUDD 2585-2606, pmptemp 2170); `glimmer_physcon.F90`: 89-102 |
| `cism_thermo_discretization.esm` | Local ESD-format rule file: the `cism_sigma_diffusion` lowering rule + the `cism_interface_flux` stencil template + the `cism_layer_thickness` integral measure, with the UPN/lev_nodes grid contract; imports the ESD `integral_lev_whole` rule (NLEV = UPN−1, lev_nodes renamed to lay_nodes) | `glissade_therm.F90`: 1324-1507, 151-163 |
| `cism_thermo.esm` | Model `IceTemperatureColumn` — the vertical heat equation ρᵢcᵢ ∂T/∂t = ∂/∂z(kᵢ ∂T/∂z) + ρᵢcᵢΦ with the ODE state on the UPN−1 = 10 INTERIOR σ nodes (PDE with the `cism_sigma_diffusion` rule; the surface and bed nodes are algebraic in CISM and are carried as temp_sfc/temp_bed, feeding the interior stencil as boundary values), the pmp profile (EqIce imports), and the basal/internal melt laws + 10 tests | `glissade_therm.F90`: glissade_therm_driver, glissade_temperature_matrix_elements 1324-1507, glissade_basal_melting_ground 1930-2012, glissade_basal_melting_internal 2032-2159, glissade_pressure_melting_point 2163-2192 |
| `cism_transport.esm` | Models `GladSmb` — the cap qsmb→acab conversion acab = qsmb·tstep_mbal·hours2seconds/rhoi·windows_per_year with the named glad/glimmer constants and the glad-level zeroing (ice-free ablation, ice-free ocean via the `cism_ocean` predicate); `MassBalance` — the SMB application: the rate law d(thck)/dt = ifelse(thck ≤ 0 ∧ net < 0, 0, net) with net = acab − bmlt, the per-step algebraic maps of the dumped entry state (melt_potential, the applied rates, the layer-1/layer-10 tracer blends via `cism_tracer_mix`, the column-mean law, the tracer zeroing of melted columns), and the ice-age increment d(age)/dt = 1 where thck > 0; `SurfaceUpdate` — calc_lsrf_usrf: lsrf = topg where grounded, lsrf = eus − (rhoi/rhoo)·thck where floating, usrf = lsrf + thck. Declares a `solver` block (stiffness low, abstol 1e-10, reltol 1e-6) for the non-smooth clip — see 'The clip and the solver block' | `glad_main.F90`: 810-824; `glad_timestep.F90`: 171-187; `glissade_mass_balance.F90`: 819-1571; `glissade.F90`: 1203-1213; `glissade_transport.F90`: 139-208; `glimmer_utils.F90`: 344-378 |
| `cism_velo_discretization.esm` | Local ESD-format rule file: the σ-column `NLAYERS = 10` / `lev` / `lev_nodes` grid contract, the `cism_sigma_measure` / `cism_sigma_midpoint` templates, and the `integral_lev_whole` / `integral_lev_nodes_cumulative_to_top` rules imported by reference from EarthSciDiscretizations with `NLEV = NLAYERS` | `glissade_velo_higher.F90`: compute_integrals_diva 5570-5696 |
| `cism_velo.esm` | Models `DivaGeometry` (the hybrid-margin surface gradient with the max_slope clamp, the 4-cell mean staggered thickness, the driving stress taudx = −rhoi·grav·stagthck·dusrf_dx, the DIVA slope angle θ = atan(\|dz_dx\|) and slope factor F = cos²θ + 4 sin²θ, and the bwat = 0 effective-pressure identity effecpress_stag = rhoi·g·stagthck), `DivaIntegrals` (the Goldberg (2011) σ-column integrals: efvs as the quad-point mean, omega_k and omega written as ESD integral operators, the `cism_stagger_masked_mean` vertex staggers, the active-cell guard), and `DivaVelocity` (beta_internal = f_ground·beta, the effective basal friction beta_eff = F·β/(F + β·ω·cosθ), btractx = β_internal·uvel_basal, the basal-velocity recovery uvel_basal = uvel_2d·beta_eff/beta_internal, and the 3-D profile uvel[k] = uvel_basal + btractx·stag_omega_k[k]·cosθ/F), plus the top-level mount of `BasalTraction` and the `variable_map` edges feeding `DivaVelocity` | `glissade_velo_higher.F90`: glissade_velo_higher_solve, load_vector_gravity 4353-4384, compute_integrals_diva 5570-5696, compute_3d_velocity_diva 5700-5821, F90:5091, 5266-5272, 5578-5583, 5756-5770, 5780, 5800-5821; `glissade_grid_operators.F90`: surface_elevation_gradient 999-1079, glissade_stagger 69-157, glissade_slope_angle_staggered 1482+ |
| `eqice_basal_templates.esm` | **TEMPORARY COPY** of the EqIce basal_templates library (5 laws: `overburden_pressure`, `mohr_coulomb_N_interior`, `basal_speed_squared`, `basal_speed`, `pseudo_plastic_drag`), bodies byte-identical to `eqice/stage2/basal_templates.esm`; imported by reference by cism_basal/cism_thermo until EqIce is upstream | PISM v2.3.2 (MohrCoulombPointwise.cc, MohrCoulombYieldStress.cc, basal_resistance.cc) |
| `eqice_ice_rheology_templates.esm` | **TEMPORARY COPY** of the EqIce ice_rheology_templates library (3 laws: `pressure_melting_temperature`, `arrhenius_softness`, `pb_softness`), bodies byte-identical to `eqice/stage2/ice_rheology_templates.esm`; imported by reference by cism_flowlaw/cism_thermo until EqIce is upstream | PISM v2.3.2 (FlowLaw.cc softness_paterson_budd, EnthalpyConverter.cc melting_temperature) |

The cell→vertex stagger of the effective pressure (`effecpress_stag`, the
plain 4-cell average with `stagger_margin = 0`, glissade_grid_operators.F90:125-133
via glissade_basal_traction.F90:1096-1098) is a grid discretization step: it
belongs to the stage-3 assembly with the velocity solve, not to these pointwise
law components. In this lane the staggered N enters `BasalTraction` as the
input `N_stag` (the traced N column of `cism_basal_traction.csv`, identical to
the `effecpress_stag` column of `cism_velo.csv` at the same vertex). The
horizontal transport and the first-order vertical remap
(`libglissade/glissade_transport.F90`) are likewise discretization operations
that need the stage-3 assembly with the EarthSciDiscretizations stencils.

## EqIce reuse map (CISM_LANE.md §2 rows 4, 8, 9)

| CISM law | EqIce import | CISM-specific additions here |
|---|---|---|
| Effective pressure (row 9): N = N0·(δP/N0)^s·10^((e0/Cc)(1−s)) | `eqice_basal_templates.esm` `mohr_coulomb_N_interior` (the identical functional form, Bueler & van Pelt 2015 Eq. 23 / Tulaczyk 2000; free names N0, W_till_max, e0, Cc resolved to the glide_types defaults) | s = clamp(bwat/bwat_till_max, 0, 1) (F90:960; applied to the water thickness before the law, algebraically identical since W_till_max > 0); **bwat = 0 → N = overburden** (F90:784 initialization; EqIce's law has no such branch — PISM always evaluates the formula); min(·, P) cap (F90:974-975); floating-cell zero (F90:982-985); P = ρᵢgH via EqIce `overburden_pressure` |
| Basal drag (row 8): β = τ_c/(u0^q·speed^(1−q)) | `eqice_basal_templates.esm` `pseudo_plastic_drag` (β = τ_c·(ε²+u²+v²)^((q−1)/2)·u0^(−q) with ε = 0.01 m/yr — exactly CISM's formula; the speed regularization is EqIce's `basal_speed`) | τ_c = tanφ·N with φ from the **bed-elevation ramp** (EqIce derives φ from data/inversion), using the cell-center topg (the stagtopg quirk, bugs.md #23); tanφ = tan(φ·π/180); the ramp parameters are glide_types defaults |
| Flow factor (row 4): flwa = E·a·exp(−Q/(RT*)) | `eqice_ice_rheology_templates.esm` `pb_softness` (same law, same constants family) | T* = clamp(T−T_pmp, −50, 0) + 273.15 (CISM's T* construction; PISM takes the pressure-adjusted temperature directly); the branch predicate is T* ≥ 263.15 K — CISM's `tempcor >= −10` (F90:2602) and EqIce's `T_pa < T_critical` agree at the boundary; **CISM's arrmll = 3.613e-13 differs from PISM's A_cold = 3.61e-13 in the third digit**, so the constants are bound at the CISM call site (A_cold = arrmll, Q_cold = actenl, A_warm = arrmlh, Q_warm = actenh, T_critical = 263.15 K, R_ideal_gas = gascon); E multiplies outside the law (F90:2579, 2603, 2605) |
| Pressure melting point | `eqice_ice_rheology_templates.esm` `pressure_melting_temperature` applied to `eqice_basal_templates.esm` `overburden_pressure`, with T_melting = 0 degC and beta_Clausius_Clapeyron = pmlt bound at the call sites | nothing (the standard Clapeyron law) |

Not in this lane (later lanes / different laws): the local-till hydrology
(row 10) and calving (row 12). The φ ramp is the only basal piece without an
EqIce counterpart; the cap coupling (row 15) is CESM-only and is covered by
`cism_cap.esm` (see the cap section below); the temperature step (row 5) is
covered by `cism_thermo.esm` (see the temperature-step section below) and
supplies the pmptemp column law pmptemp = −ρᵢg·pmlt·thck·stagsigma (F90:2170)
internally.

## Pointwise laws — Fortran vs esm transcription notes

- **PSEUDO_PLASTIC_OLD** (glissade_basal_traction.F90:250-292): bed =
  topg − eus with the **cell-center topg** (the code comment at F90:261-263
  admits the stagtopg bug; kept for backward compatibility, recorded in
  bugs.md #23); φ from the ramp; tanφ = tan(φ·π/180); τ_c = tanφ·N_stag;
  speed = √(u²+v²+smallnum²) with smallnum = 1.0d-2 m/yr (F90:136,189);
  β = τ_c/(u0^q·speed^(1−q)). The speed cap of F90:196-198
  (`beta_powerlaw_umax > 0`) and the β floor/f_ground weighting of
  F90:597-633 (`beta_grounded_min = 100 Pa yr/m`, `which_ho_beta_limit = 0`)
  are applied **after** the traced law: not transcribed here because they are
  solve-assembly steps (and no-ops on the traced rows: f_ground = 1 and
  β > 100 at every traced vertex; beta_powerlaw_umax = 0 in the config).
- **BWAT_BVP** (F90:944-985): relative_bwat = max(0, min(bwat/bwat_till_max,
  1)); N = N0·(δ·P/N0)^relative_bwat·10^((e0/Cc)(1−relative_bwat)); N =
  min(N, P); floating cells zeroed. The bwat > 0 branch only fires where
  bwat > 0 — **the reference run never fires it** (see How N is pinned).
  The final 0..P capping of F90:1082-1086 is a no-op on this construction and
  is documented, not transcribed.
- **FLWA_PATERSON_BUDD** (F90:2585-2606): tempcor = min(0, temp − pmptemp),
  tempcor = max(−50, tempcor); flwa = E·a·exp(−Q/(R·(tempcor + 273.15))) with
  the (a, Q) pair switching at tempcor ≥ −10; E = flow_enhancement_factor
  (1.0 in the reference run). The waterfrac enhancement of F90:2608-2610 and
  the damage branch are off in the reference config (waterfrac = 0,
  no damage).

## How N is pinned (the effecpress trace is empty)

`cism_effecpress.csv` has **no rows in this run**: the gland20 initial state
has bwat = 0 everywhere and the 2-year cold-base run produces no basal melt
(cism_thermo.csv: bwat = bmlt_ground = 0 on all 256 rows), so the bwat > 0
branch of glissade_calc_effecpress never fires. The effective pressure is
therefore pinned three ways:

1. **The bwat = 0 branch from the dumps**: with bwat = 0, N (cell) =
   overburden P = ρᵢg·thck (the F90:784 initialization). The Effecpress tests
   `bwat0_overburden_*` assert N = ρᵢg·thck at the velo-entry thck of three
   dump cells (steps 1 and 4).
2. **The staggered value from the dumps**: at every traced vertex,
   effecpress_stag = ρᵢg·stagthck to 2e-16 relative (verified across all 4
   velo-entry rows checked: cells 25,11 / 26,12 / 28,11 step 2 / 19,13 step 4)
   — the full N chain (N = ρᵢgH → stagger_margin=0 4-cell average) is
   consistent with the traced N columns; the stagger itself is the stage-3
   grid operation.
3. **The bwat > 0 branches from the source formula**: the tests
   `bwat_partial_plateau`, `bwat_saturated`, `bwat_clamped` and
   `floating_zero` evaluate the F90:960-985 law at the glide_types defaults
   with synthetic bwat (0.5, 2.0, 3.0 m) and dump-derived overburden; expected
   values are hand-computed from the source formula (documented in each test
   description), in the same spirit as the MOSART Fortran-driver pins for
   branches the dump window never reaches. The trace machinery fires whenever
   bwat > 0 (a wetter stage-1 window would replace these with dump pins).

## Cap coupling boundary (row 15) — `cism_cap.esm`

The mediator↔CISM coupling layer around `glad_gcm` (CISM_LANE.md §2 row 15,
"missing in EqIce — CESM-specific coupling conversions"). The stage-1 dump
`stage1/dumps/cism/cism_cap.csv` traces the boundary in `source_glc/glc_RunMod.F90`
(instrumented fork): **import** rows with tsfc [degC] and qsmb [kg/m²/s]
*exactly as consumed by CISM* (after the cap `import_fields` processing), and
**export** rows with topo [m], hflx [W/m²], rofi/rofl [kg/m²/s = mm H₂O/s] and
ice_covered [−] *exactly as produced* by `glad_set_output_fields`. Two cap
steps (one per coupling year, GLC_NCPL=1; tstep_mbal = 8760 h; each year = 2
CISM dynamics steps of dt = 0.5 yr), 32 cells × 2 phases × 2 steps = 128 rows.

### Import law (`GlcCapImport`)

`glc_import_export.F90:326-333` (`import_fields`):

- **tsfc** = ifelse(tsfc_c < −250, 0, tsfc_c) with the named intermediate
  tsfc_c = tsfc_raw − tkfrz: tsfc_raw the mediator `Sl_tsrf` [K] and tkfrz =
  SHR_CONST_TKFRZ = 273.15 K (`glc_constants.F90:22-23`,
  `shr_const_mod.F90:37`). The −250 °C clamp is the "Jer hack" sentinel
  (F90:328-332): CLM-ocean land points send raw 0 K, which the conversion
  would drop to −273.15 °C; the clamp resets them to 0 °C. The conversion is
  the consumer's named intermediate (written once), and the clamp is the
  one-parameter `cism_tsfc_import` template — the Fortran converts and clamps
  in place on the same array.
- **qsmb** passes through unchanged (`state_getimport`, F90:679-714, applies
  no conversion).
- The mediator→model halo conversion (`vector_to_spatial` /
  `parallel_convert_nonhaloed_to_haloed`) is a grid operation, deferred to
  the stage-3 assembly.
- The model-side continuation — the mbal averaging of qsmb/tsfc over the
  coupling interval and the qsmb → acab conversion
  (acab = qsmb·tstep_mbal·3600/1000·rhow/rhoi m/yr ice,
  `glad_main.F90:800-824`, `glad_timestep.F90:200`, `glide_set_acab`) — is the
  surface-forcing/transport lane's boundary (already pinned by the stage-1
  `smb_chain` validator check and the transport dumps), not the cap's.

### Export law (`GlcCapExport`)

`glad_set_output_fields` (`glad_main.F90:971-1027`) → `set_output_states`
(`glad_output_states.F90:48-158`) + `glad_output_fluxes`
(`glad_output_fluxes.F90:77-172`):

- **topo** = usrf on the active grid (usrf > 0, `is_in_active_grid`
  F90:110-133), 0 elsewhere (F90:91, 97-99).
- **ice_covered** = 1 where thck > min_thck = 0 (`glad_constants.F90:63`,
  F90:152) **and** usrf > 0, else 0 (F90:82-86).
- **rofi** = coupling-interval average of ρᵢ·(calving_thck + removal_thck)/dt
  (F90:104-105) — calving + the other solid-ice removals (icebergs/isthmuses,
  forced retreat, tiny-thickness cleanup, F90:91-96).
- **rofl** = coupling-interval average of ρᵢ·bmlt + ρᵢ·melt_thck/dt (F90:128-129)
  — basal melt + lateral/submarine melt.
- **hflx** = coupling-interval average of −ucondflx (F90:138-139; ucondflx
  positive down, `glissade_therm.F90:669`). The sign-flip comment says
  "positive up", which contradicts the "positive down" labels of
  `glad_main.F90:625` and the CMEPS field dictionary — see `bugs.md` (#32
  entry); unpinned by the dump because hflx = 0 throughout this window.

**Time average as a storage ODE, not `integral`.** The Fortran forms each flux
as a sum over the n_icetstep = 2 dynamics steps of the coupling interval
(reset at the interval start, `glad_timestep.F90:149`; one accumulation per
`glissade_tstep`, F90:290), then divides by the step count
(`calculate_average_output_fluxes`, F90:151-159). With piecewise-constant
integrands over the steps, Σₖ fₖ/count = (1/T_coupling)·∫₀^T_coupling f dt
with T_coupling = count·dt = 8760 h = 3.1536e7 s — so the export is the
coupling-interval average in the integral sense. It is written here as the
storage ODEs D(S_rofi)/dt = ρᵢ·(calving_rate + removal_rate) (and
analogously S_rofl, S_hflx) with the export S/T_coupling: the `integral` op of
esm 1.3.0 is a *spatial* PIDE term that the CLI cannot lower or evaluate
(probed: `unlowered_operator`, esm-spec §4.2/§9.6.8 — a file using it
validates but cannot simulate), so the differential form — the same
continuous-limit pattern the MOSART lane uses for `D(wh)/dt` — is used
instead. The ODE is exact for the Fortran's piecewise-constant integrands.
Model time is in seconds (the dump's thour column × 3600; cap_step 1 at 0 h,
cap_step 2 at 8760 h). The Fortran exports the average only at the interval
end; the .esm evaluates S/T_coupling at any t, and the reset of S at each
interval boundary is a discrete event for the stage-3 integration.

**Inputs the dump does not carry.** The integrands — calving_thck,
removal_thck, bmlt, lateral_melt_thck, ucondflx — are per-dynamics-step CISM
state; `cism_cap.csv` carries only the final interval averages (all zero in
this window: the 2-year cold-base run has no calving, no basal/lateral melt,
and ucondflx = 0). They enter as parameters. usrf/thck are likewise not in
the cap dump and are taken from the transport dumps in the topo tests
(cross-file identity topo == transport-exit usrf verified on all 64 export
rows). `ice_sheet_grid_mask` (the same active-grid gate) is not dumped and is
not modeled.

### Cap tests (22 tests, 48 assertions)

All instantaneous pins are rel 1e-12 `time_span {0,0}`; the average-law pins
integrate one coupling interval [0, 3.1536e7 s] with constant integrands
(linear ODEs — exact to round-off).

| Test group | Tests | What is pinned |
|---|---|---|
| Import conversion (dump) | 5 | tsfc = clamp(tsfc_raw − 273.15) and the qsmb pass-through at cap_step 1/2 rows (cells 25,11 / 28,11 / 19,13 / 27,14 / 28,12, incl. the negative-SMB cell 28,11, the coldest −9.678 °C and warmest −0.004 °C traced tsfc); the raw input is rebuilt as tsfc_raw = tsfc_dump + tkfrz and the intermediate round-trips to tsfc_c = tsfc_dump — the clamp did not fire on any traced row (pinned synthetically) |
| Import clamp (source formula) | 4 | the sentinel firing (raw 0 K — the Jer hack case), the strict-< boundary at exactly −250 °C, just below it, and the warm linear branch |
| Export states (dump) | 6 | topo = usrf and ice_covered = 1 at both cap steps (cells 25,11 / 26,12 / 19,13 / 28,11), incl. the thinnest traced ice (thck = 1.65e-5 m, the min_thck = 0 boundary exercised from above) |
| Export states (source formula) | 2 | the ice-free branch (thck = 0 → ice_covered = 0, topo unchanged) and the inactive-grid branch (usrf = 0 → topo = ice_covered = 0) |
| Flux averages | 5 | the dumped zeros (all 64 export rows have hflx = rofi = rofl = 0; the storage ODEs integrate over the full [0, T_coupling] so the integration is exercised, not just the t = 0 defaults), and the averaging law with constant integrands: solid (rofi = ρᵢ·1.75 m/yr, incl. an intermediate S_rofi pin at t = T/2), liquid, heat, and all three combined |

### Coverage limitations (this data window)

- **The −250 °C clamp never fires in the dumps** (traced tsfc ∈ [−9.68,
  −0.004] °C — all CLM land/ice points); pinned synthetically from the source
  formula.
- **rofi/rofl/hflx are identically zero** in the window (no calving, no
  melt, ucondflx = 0); the dump pins the zero case, the averaging law is
  pinned with synthetic constant integrands. A wetter/dynamic stage-1 window
  would dump-pin nonzero averages (and the hflx sign).
- **ice_covered = 1 on every export row** and usrf > 0 everywhere; the
  ice-free and inactive-grid branches are source-formula pins.
- **hflx sign convention** is unpinned (all zeros) — bugs.md #32.
- The qsmb → acab SMB chain is the transport lane's pin
  (`smb_chain` validator check; the cap dump's qsmb is its input).

## Temperature step (CISM_LANE.md §2 row 5)

The prognostic temperature solve (`whichtemp = TEMP_PROGNOSTIC`) as
`IceTemperatureColumn` (`cism_thermo.esm`): the vertical heat equation
ρᵢcᵢ ∂T/∂t = ∂/∂z(kᵢ ∂T/∂z) + ρᵢcᵢΦ on CISM's σ coordinate, the pointwise laws of
the step, and the local EarthSciDiscretizations-format rule that lowers the
diffusion term (`cism_thermo_discretization.esm`, a candidate ESD PR). Pinned
against the extended `cism_thermo.csv` dump (256 rows × 78 columns). The dump
extension — `dttem`, `dissip_001..010`, `sigma_001..011`, `stagsigma_001..010`
(fork branch `instrument/cism-dumps`, committed at all three submodule levels,
not pushed) — was required: a dissip = 0 reproduction fails the operator
identity at 3–5% of the step-1 temperature change, and a pin needs the exit
rows and the grid geometry, not just the entry rows. Source:
`glissade_therm.F90` (glissade_therm_driver, glissade_temperature_matrix_elements
F90:1324-1507, glissade_basal_melting_ground F90:1930-2012,
glissade_basal_melting_internal F90:2032-2159, glissade_pressure_melting_point
F90:2163-2192).

### Representation: interior-state formulation (algebraic boundary nodes)

CISM's temperature lives on UPN+1 = 12 nodes per column: the surface and bed
faces σ(1) = 0, σ(11) = 1 and the 10 layer midpoints stagsigma(1..10)
(glide_calc_sigma, SIGMA_COMPUTE_GLIDE). The matrix rows 1 and upn+1 are NOT
diffusion rows: the surface row pins the Dirichlet value T(0) = min(0, artm)
(F90:843, 1405-1408) and the bed row is the algebraic basal BC (frozen: flux
balance; thawed: T = bpmp; F90:1473-1476, 1497-1504). The model therefore
carries the ODE state on the UPN−1 = 10 INTERIOR nodes (index set `lev` = the
dumped temp_001..010 columns) with d(temp)/dt = dTdt, and the two boundary
nodes as the algebraic unknowns `temp_sfc` and `temp_bed` — so the boundary
laws FEED BACK into the interior operator: the rule's stencil reads
temp_sfc/temp_bed at the surface/bed faces instead of state rows. The full
12-node geometry (sigma_node = [0, stagsigma, 1], the layer widths dsigma_cv,
thck) enters as parameters from the dump extension.

The diffusion term lowers, by the rule `cism_sigma_diffusion` (priority 10,
match-scoped to lev-shaped u, fired on D(coni·D(temp, lev), lev)), to the exact
dups stencil of glissade_temperature_matrix_elements: L·T[n] =
(G[n+2] − G[n+1])/(thck·dsigma_cv[n+1]) over the CISM layer control volume
[σ(n), σ(n+1)], with the interface flux G[m] =
cism_interface_flux(coni, u[m−1], u[m], σ[m−1], σ[m], thck) =
coni·(u[m] − u[m−1])/(thck·(σ[m] − σ[m−1])) (the node-to-node gradient;
stagsigma(0) ≡ σ(1), stagsigma(upn) ≡ σ(upn), F90:151-163) — at the surface
face m = 1 and the bed face m = UPN+1 the flux reads the algebraic
temp_sfc/temp_bed. The same pointwise cism_interface_flux builds the basal
conductive flux lcondflx_ground = −G[UPN+1] (F90:1070). dTdt = L·T/(ρᵢcᵢ) +
dissip on the interior nodes; there are no boundary dTdt rows.

### Pin strategy: the instantaneous operator identity at the exit state

The CISM step is ONE backward-Euler solve of the tridiagonal system
(diag − dt·L)·T_exit = T_entry + dt·Φ. The harness integrators cannot reproduce
one exact stiff backward-Euler step, so the model does not integrate the
entry→exit map; instead the operator is pinned INSTANTANEOUSLY at the exit
state: for the interior rows the identity (T_exit − T_entry)/dttem = L·T_exit + Φ
holds exactly (the dump round-trips at 1e-13..6e-9 relative), and
`IceTemperatureColumn`'s dTdt is asserted against (T_exit − T_entry)/dttem with
the ICs set to the exit profile and the dumped dissip. Five such pins span the
window (step 1 and step 4 of cell (27,14), step 1 of cells (27,13) and (25,11) —
the thinnest solved column at thck = 1.3 m — and step 4 of cell (19,13));
per-node residuals ≤ 6.0e-9 relative → rel 1e-8 tolerances. The algebraic laws
of the step are pinned at rel 1e-12: temp_sfc = min(0, artm) (all 128 exit
rows), bpmp = −ρᵢg·pmlt·thck, the frozen basal row temp_bed = temp[10] +
(bfricflx + bhydroflx − bheatflx)·dsigbot·thck/coni (the bflx_src named
intermediate; closes to 6.2e-17 W/m² in flux terms), lcondflx_ground =
−cism_interface_flux at the bed face (one thin-column assertion at rel 1e-8 —
there the flux is the ratio of two cancelling differences, temp_bed − temp[10]
= 2.5e-8 K from ~1.87 degC values), and the melt chain bmlt_raw / bmlt_ground /
btemp_ground (all-zero in the window: the closed basal balance |bfricflx +
lcondflx_ground| ≤ 1.5e-15 W/m² is zeroed by the eps11 guard on the FLUX,
F90:1927-1929). The operator pins close the boundary feedback: the stencil at
the node adjacent to each boundary reads temp_sfc/temp_bed, so the pinned
algebraic values enter the asserted dTdt. The thawed-bed row, the freeze-on
branch, and the internal-melt cap are pinned synthetically from the source
formulas (see gaps); the thin-cell surface law is pinned on cell (28,11)
(thck = 0.02 m ≤ thklim_temp = 1 m: the column bypasses the solve — temp =
min(0, artm) on all nodes, bpmp/temp_bed/bflx_mlt/btemp_ground = 0 via the
ice_mask gate, dissip = 0).

### Files

| File | Contents | Fortran source |
|---|---|---|
| `cism_templates.esm` | +4 thermo templates: `cism_internal_melt_rate` (layer melt-energy density rhoi·shci·max(0, T−Tpmp), J/m³ — the Fortran's melt_energy per unit layer thickness), `cism_basal_melt_rate` (flux-balance rate bflx_mlt/(lhci·rhoi) with the eps11 zeroing applied to the FLUX, F90:1927-1929), `cism_basal_freeze_on` (predicate bmlt < 0 ∧ bwat = 0 ∧ btemp ≥ bpmp − eps08, evaluated ONCE as the named intermediate freeze_on), `cism_basal_melt_cap` (btemp_ground = min(btemp, bpmp); freeze-on drop to min(bpmp − 2·pmp_threshold, bpmp − eps08), consuming freeze_on). The Clapeyron pmp is imported from the temporary EqIce copies (no CISM pmp template) | `glissade_therm.F90`: 2127-2155, 1930-1967, 2007-2012, 1979 |
| `cism_thermo_discretization.esm` | Local ESD-format rule file: the `cism_sigma_diffusion` lowering rule (interior-state stencil, algebraic boundary values via the free-name contract), the pointwise `cism_interface_flux` face-flux template (reused by the stencil and by lcondflx_ground), `cism_layer_thickness` (dz = thck·dsigma_cv, the integral measure), and the ESD `integral_lev_whole` import (NLEV = UPN−1, lev_nodes renamed to lay_nodes) | `glissade_therm.F90`: 1324-1507, 151-163 |
| `cism_thermo.esm` | Model `IceTemperatureColumn` (the ODE d(temp)/dt = dTdt on the UPN−1 interior nodes, the algebraic temp_sfc/temp_bed laws feeding the operator, the pmp profile from EqIce, the basal/internal melt laws — melt_internal as the whole-column integral of f_melt over z_sfc..z_top divided by rhoi·lhci·dttem, F90:2127-2155) + 10 tests (5 operator pins, thin-cell law, thawed bed, internal melt, freeze-on, eps11) | `glissade_therm.F90` as above |

### Gaps and deferrals

- **Entry→exit backward-Euler map**: not integrated (see pin strategy). A
  stage-3 coupling that must reproduce the CISM trajectory exactly would
  transcribe the tridiagonal solve itself (an algebraic map, not an ODE) or
  accept first-order-in-dt trajectory differences.
- **Surface re-pin jump**: the surface value is re-pinned inside the solve, so
  the ENTRY rows of steps 1 and 3 carry the previous value (the initial 0 and
  the previous year's artm); the surface law is pinned on the exit rows (all
  128).
- **waterfrac ≡ 0**: TEMP_PROGNOSTIC carries no internal water (the dump's
  waterfrac columns are all 0); the enthalpy scheme (waterfrac > 0) is deferred.
- **Thin-ice bypass**: cells with thck ≤ thklim_temp = 1 m skip the solve
  entirely (temp = min(0, artm) on all nodes; F90:1205). Transcribed as the
  ice_mask gate: the pmptemp/bpmp/temp_bed/bflx_mlt/btemp_ground laws return 0
  where ice_mask = 0 (the dumped column), so the gate itself is dump-pinned on
  cell (28,11); the temp = min(0, artm) column profile is pinned there too.
  21 of 128 traced rows are such cells.
- **Thawed-bed / freeze-on / internal-melt branches**: no traced row exercises
  them (the window is frozen: every entry temp_011 < bpmp − pmp_threshold, worst
  margin −1.425 degC; bmlt = bwat = bheatflx = bhydroflx = 0); pinned
  synthetically from the source formulas.
- **Internal-melt temperature reset**: the melt cap's paired state update
  temp → pmptemp (F90:2148) is a post-solve mutation, documented, not
  transcribed as an equation.
- **Internal melt added to bmlt_ground** (F90:1298): melt_internal is computed
  as the whole-column integral but is NOT added into bmlt_ground in the model —
  the sum is unpinned because melt_internal = 0 throughout the window; the
  integral itself is pinned synthetically (the internal_melt_rate test
  integrates the melt-energy density over the bottom layer).
- **pmptemp profile**: the dump carries only the bed value bpmp; the per-node
  profile law −ρᵢg·pmlt·thck·sigma_node is computed via the EqIce imports
  (pressure_melting_temperature of overburden_pressure), pinned at the bed
  against the dumped bpmp and synthetically for the interior (the
  internal_melt_rate test).
- **Dump extension**: dttem, dissip_001..010, sigma_001..011, stagsigma_001..010
  were added to `cism_thermo.csv`; sigma/stagsigma are derivable from the grid
  construction, but dissip and dttem are not, and the operator pins need them.

## Transport slice (CISM_LANE.md §2 rows 6/7/14) — `cism_transport.esm`

Stage-2 implementation of the pointwise slice of the CISM thickness/tracer
solve (`glissade_thickness_tracer_solve`): the parts of the solve that are
**not** grid operations — the cap qsmb→acab conversion and glad-level zeroing
(row 14 boundary), the surface/basal mass-balance application with its tracer
blends (rows 6/7), and the post-transport surface/bed elevation cleanup. The
horizontal transport and the first-order vertical remap
(`libglissade/glissade_transport.F90`) are discretization operations: they
need the stage-3 assembly with the EarthSciDiscretizations stencils and are
**out of scope here**. What is pinned: the SMB application, the temperature
surface/bed handling, the thickness clip, and the cleanup, against the
stage-1 instrumentation dump `stage1/dumps/cism/cism_transport.csv`
(phase=entry at the top of the solve, before horizontal transport;
phase=exit after the SMB application and tracer cleanup) and the cap dump
`cism_cap.csv` (the qsmb import).

The per-cap-step chain is: the cap accumulates qsmb (kg m⁻² s⁻¹) over the
mass-balance window `tstep_mbal`; glad converts it to a per-window thickness
and zeroes it where ice-free/ablation and in ice-free ocean
(`glad_timestep.F90:171-187`, `glad_main.F90:810-824`); each dynamics step
(0.5 yr) applies acab and bmlt to the ice column
(`glissade_apply_smb` → `mass_balance_driver` →
`add_surface_and_basal_mass_balance`, `glissade_mass_balance.F90:819-1571`),
remaps the tracers vertically, and recomputes lsrf/usrf
(`calc_lsrf_usrf`, `glimmer_utils.F90:344-378`). The stage-1 fork traces the
chain at the top of the solve (entry) and after the SMB application and
cleanup (exit).

- **GladSmb** (`glad_main.F90:810-824`): the Fortran accumulates qsmb in
  kg m⁻² s⁻¹ over `tstep_mbal` hours and converts to m of ice: acab =
  qsmb·tstep_mbal·hours2seconds/rhow·rhow/rhoi·windows_per_year with the
  named glad/glimmer constants hours2seconds = 3600 s/h
  (`glad_constants.F90:44`), rhow = SHR_CONST_RHOFW = 1000 kg/m³, rhoi =
  SHR_CONST_RHOICE = 917 kg/m³ (`glimmer_physcon.F90:41-43`) and scyr =
  31536000 s/yr (`glimmer_physcon.F90:90`) driving the DERIVED
  windows_per_year = scyr/(tstep_mbal·hours2seconds) [1/yr; = 1 in the
  reference run, tstep_mbal = 8760 h = 1 yr]. The water density cancels
  (qsmb·3600·tstep_mbal is kg/m², ÷ rhoi is m), so the esm conversion is
  acab = qsmb·tstep_mbal·hours2seconds/rhoi·windows_per_year (the dump
  column is labeled m/yr ice). The `acab = −thck_temp` cap of
  `glad_timestep.F90:233-235` ("Adjust glad acab for output") is applied
  **after** `glide_set_acab` and only modifies the glad output variable,
  never the applied field — documented, not transcribed.
- **The zeroing** (`glad_timestep.F90:171-187`): acab = 0 where thck ≤ 0
  and acab < 0 (ice-free ablation), and where GLIDE_IS_OCEAN (bed below sea
  level, no ice). Both fire independently; the esm `acab` equation nests
  them via the shared `cism_ocean` predicate template.
- **MassBalance** (`glissade_mass_balance.F90:819-1571`): the Fortran walks
  the 10 layers: accumulation adds s = acab·dt to layer 1 with the blend
  T1′ = (H₁·T1 + s·Tsurf)/(H₁ + s) at the surface temperature; ablation
  removes the top layers without changing their tracers; freeze-on adds
  b = −bmlt·dt to layer 10 with the blend at the bed temperature; basal
  melt removes the bottom layers; fully removed layers get tracer = 0; a
  column that melts out entirely accumulates the excess into
  `melt_potential` (a per-step diagnostic, zeroed at the top of every
  `mass_balance_driver` call). The esm model is the layer-free continuous
  limit with the RATE-LEVEL guard d(thck)/dt = ifelse(thck ≤ 0 ∧ net < 0,
  0, net): the derivative stops exactly when the column is gone, which
  reproduces the Fortran's layer-by-layer removal that clips the column at
  zero (F90:1398-1419) without carrying a deficit into the next step (the
  Fortran re-derives thck from the clipped layers every step; the
  state-level max(0, thck_raw) of the first review round could not express
  the per-step zeroing and was replaced by the guard). The per-step
  bookkeeping — melt_potential = max(0, −(H_pre + net·dt)), the applied
  rates acab_applied/bmlt_applied, and the tracer blends — is written as
  algebraic maps of the dumped entry state (parameter `H_pre` + entry
  tracers), evaluated at the step start, which is exact at the step
  boundaries the tests evaluate. The tracer laws keep the layer-1
  temperature, the layer-10 temperature, and the column mean through the
  shared `cism_tracer_mix` law T′ = (H·T + s·T′)/(H + s); the ablation
  mean law removes the ablated ice at the layer-1 temperature, which is
  exact while the ablation stays within layer 1. The column-mean output
  Tbar_out is derived (the Fortran has no column-mean state) and is
  UNPINNED: no dump case discriminates it (see 'Coverage gaps').
- **The ocean mask**: `glissade_get_masks` with thklim = 0
  (`glissade_masks.F90:70-96`) — ocean = topg < eus ∧ thck ≤ 0 (the same
  `cism_ocean` predicate as the glad zeroing); the SMB and BMB are not
  applied there and melt_potential stays 0 (the `effective_areafrac = 0`
  of ocean cells, F90:893-901). No traced cell is ocean (all 32 cells have
  topg ≥ 921 m): unpinned (see 'Coverage gaps').
- **The ice age** (`glissade.F90:1203-1213`): age += dt where thck > 0,
  age = 0 on ice-free cells; the tracer surface/bed values age_usrf = 0,
  age_lsrf = bottom-layer age (`glissade_transport.F90:204-208`). The esm
  ODE stops the increment on ice-free cells but cannot express the
  instantaneous reset — documented difference (unobservable in the window,
  which never refreezes). Unpinned: there is no age column in the dump
  (see 'Coverage gaps').
- **SurfaceUpdate** (`glimmer_utils.F90:344-378`): lsrf = topg where
  grounded (topg − eus ≥ −(rhoi/rhoo)·thck), lsrf = eus − (rhoi/rhoo)·thck
  where floating, usrf = lsrf + thck.

### The clip and the solver block

The rate guard is discontinuous at the clip, and the runner's default
diffsol BDF arm chokes at the crossing:

```
simulate failed: diffsol error: ODE solver error: Step size is too small at time = 0.000009337527277120702
```

(exactly the analytic crossing time thck0/|net| of the failing cell; the
same error appears at 0.07506885100532178 for the cell-29,13 pin). The
document therefore declares its own numerics in the esm-spec §2.2 `solver`
block — `stiffness: "low"` (the system is not stiff; the declaration steers
other bindings to their explicit-RK arms, which cross the kink cleanly) and
integration tolerances `abstol: 1e-10, reltol: 1e-6` (the runner's default
reltol 1e-4 / abstol 1e-6 leaves the BDF step control stuck at the kink). Under these
tolerances the BDF arm crosses the clip and freezes the state at its first
negative step, landing within ~1e-9 m of 0; the full-melt pins assert
thck = 0 with abs tolerance 1e-8 — still eight orders below the −0.35 m the
unguarded law would give, so the clip stays fully discriminating. The ERK
arm (TSIT45) lands on the clip to round-off and passes the same assertions.

### How each law is pinned

The dump rows: 32 cells × 4 steps × 2 phases, upn = 11 (temp_000 = surface,
temp_001..010 = layers 1..10, temp_011 = bed); steps 1,2 ↔ cap_step 1,
steps 3,4 ↔ cap_step 2; every row carries time = 0.5·step yr. All test
inputs come from the dump entry phases (the pre-step state) and the
cism_cap.csv import rows; all expected values come from the dump exit
phases or from the source formulas (stated in the test descriptions).
22 tests / 29 assertions.

| Law | Pin | Tolerance | Residual |
|---|---|---|---|
| qsmb→acab | cap import rows at cells 25,11 (cap steps 1, 2) and 26,11 (cap step 1) | rel 1e-12 | ≤ 1.1e-16 (exact) |
| glad zeroing (ice-free ablation) | cell 28,11 step 2 (thck = 0, raw acab = −0.7075 → 0; dump-verified); cell 28,11 step 3 (thck = 1.65e-5 > 0 → NOT zeroed; dump-verified) | exact | 0 |
| d(thck)/dt = net | accumulation: cell 26,11 steps 1/2/4 (rel 1e-6 / 1e-9 / 1e-6); ablation: cell 25,11 steps 3/4 (rel 2e-4 / 3e-4); full melt: cells 28,11 (steps 1, 3), 29,13 (step 1) → thck = 0 (abs 1e-8) | see left | advection: +9.6e-7 / +9.6e-10 / −1.7e-6 m; +1.8e-5 / +1.9e-5 m; BDF clip landing ≤ ~1e-9 m |
| full-melt clip + melt_potential | cells 28,11 (s1: 0.33379825832198695 m; s3: 0.8837069386692753 m), 29,13 (s1: 0.22629773133841058 m) — melt_potential is derived (not a dump column), from the Fortran's own conservation check | rel 1e-12 | exact |
| acab_applied | clipped branch only (cell 28,11 s1: acab_applied = −H_pre/dt = −0.039916999638080596) | rel 1e-12 | exact (source-formula: not dump columns) |
| layer-1 accumulation blend | cells 26,11 s1 (rel 1e-5), 27,11 s1 (rel 7e-5; s > H₁), 25,12 s1 (rel 2.5e-4) | see left | −1.7e-6 / −1.2e-5 / −3.6e-5 degC |
| tracer zeroing on full melt | cell 28,11 s1 (T1_out = 0; the dump exit temp_001..010 = 0 of the melted columns) | exact | 0 |
| lsrf/usrf | grounded: cell 25,11 s1; ice-free: cell 28,11 s3 (thck = 0); thin ice: cell 29,13 s4 | rel 1e-12 | exact (all 256 rows: lsrf = topg, usrf = lsrf + thck to round-off) |

The tolerance convention: 1e-12 where the esm law is exact (pure algebra
or linear-in-time integration), 5-10× the measured residual where the dump
exit value carries a physics effect the pointwise model does not contain
(horizontal advection, vertical remap) — each documented in the test
description, following the MOSART scheme-gap pattern. The full-melt abs
1e-8 covers the BDF arm's clip landing (see 'The clip and the solver
block').

### Coverage gaps (deferred or unpinnable)

- **Horizontal transport and the first-order vertical remap** are grid
  operations (`glissade_transport.F90:100-272, 1101-1441`) and belong to
  the stage-3 assembly with EarthSciDiscretizations. Their fingerprints are
  the advection residuals of the thickness and tracer pins above.
- **bmlt = 0 everywhere in the reference window** (bmlt_ground =
  bmlt_float = 0 on all 256 rows: the cold base melts nothing and the
  window has no shelf ice). The basal-melt (> 0) removal branches, the
  freeze-on (< 0) branch, and the layer-10 freeze-on blend T10′ are
  transcribed from the source but have no dump pin (the earlier
  freeze_on_blend_synthetic was dropped — synthetics pin the esm code, not
  the Fortran).
- **No ocean cells in the window** (all 32 traced cells have topg ≥ 921 m):
  the glad-level ocean zeroing branch and the glissade-level
  ocean-no-application branch (net_applied = 0, melt_potential = 0) are
  transcribed from the source but unpinned (the earlier
  zeroing_ocean_synthetic and ocean_no_application_synthetic were dropped).
- **No ice-age dump column**: the age increment is transcribed from the
  source but unpinned (the earlier age_increment_synthetic and
  age_ice_free_stops synthetics were dropped); the Fortran's ice-free reset
  (age = 0) is not representable as an ODE — documented in the equation
  description.
- **Floating lsrf/usrf**: the window has no shelf ice (the
  floating_synthetic was dropped); the floating branch
  lsrf = eus − (rhoi/rhoo)·thck is transcribed from the source and
  unpinned. All 256 dump rows satisfy lsrf = topg (grounded).
- **Tbar_out / the ablation mean law**: dropped. The column-mean output is
  derived (the Fortran has no column-mean state), and the search of all 256
  dump rows found no ablation case where the law's effect exceeds the
  advection residual (best ratios ~2.2, effects 4.7e-4 / 3.8e-5 relative —
  inside the advection-dominated tolerance), so no assertion could
  discriminate the law: the earlier ablation_mean_step3/step4 pins passed
  even with Tbar_out = Tbar and were removed. The law itself remains in the
  model (Tbar_blend → Tbar_out via the composed cism_tracer_mix law),
  documented and unpinned.
- **acab_applied / bmlt_applied / melt_potential / ocean cells** are not
  dump columns (the dump carries the input acab, not the applied values):
  pinned from the source formulas at the full-melt and applied-rate pins.
- **The layer-free tracer model** tracks the layer-1 temperature, the
  layer-10 temperature, and the column mean only; the per-layer identities
  (which layer the ablation bites into, the full 10-layer profile) are
  stage-3 vertical-discretization content (the sim1d.py reproduction of
  the exit profiles lives in `scratch/` as the verification reference).
- **The (c) 'acab = −thck_temp' output cap** of glad_timestep.F90:233-235
  is not transcribed (it adjusts the glad output variable only, after
  glide_set_acab; the applied field is unaffected) — documented.
- **melt_potential accumulation across steps**: the Fortran resets
  melt_potential to 0 at the top of every mass_balance_driver call; the
  esm form max(0, −(H_pre + net·dt)) is the per-step value of the dumped
  entry state. The two agree for any single step starting from an
  ice-covered state, which is what the tests pin.

## The DIVA velocity-solve pointwise laws (CISM_LANE.md §2 row 2) — `cism_velo.esm`

Stage-2 implementation of the **pointwise** laws of CISM's higher-order
velocity solve (reference run with `which_ho_approx = HO_APPROX_DIVA = 4`,
`which_ho_gradient = 0`, `margin = 1`, `assemble_taud = 1`,
`diva_slope_correction = .true.`, `which_ho_babc = 3`, `max_slope = 0.5`,
`thck_gradient_ramp = 50 m`, `thklim = 1 m`, `upn = 11` (10 layers),
`nquad = 4`, `rhoi = 917`, `g = 9.80616`, `dx = dy = 20000 m`). What is
implemented here is everything in `glissade_velo_higher_solve` that is a
*local law of a single column*: the vertex geometry of the driving stress and
the DIVA slope correction, the Goldberg (2011) vertical integrals of the
σ column, the vertex staggers of `omega`/`omega_k`, and the 3-D velocity
recovery from the 2-D solve. The finite-element matrix solve itself and the
Picard iteration are **stage-3 assembly** concerns and are left to the full
`eqesm.esm`; this file feeds their inputs from the dumps. Sources:
`stage1/cesm/components/cism/source_cism/libglissade/`
(`glissade_velo_higher.F90`, `glissade_grid_operators.F90`,
`glissade_basal_traction.F90`, `glissade_calc_effecpress`), the run namelists
(`stage1/cime/cases/T1850Gg_cism_inst/run/cism.gris.config`), the stage-1
instrumentation dumps `stage1/dumps/cism/cism_velo.csv` + `cism_velo_nbr.csv`
+ `cism_velo_stag.csv` + `cism_basal_traction.csv` +
`cism_velo_marine_exit.csv` + `cism_velo_marine_pass.csv` (2026-11-09 rerun,
jobs 12277039/12279575/12279827/12279959, fork commits `db624ec`/`76d170d`
with the uncommitted stage-2 dump additions `cism_dump_velo_stag` /
`cism_dump_marine_pass` / `cism_dump_marine_exit`, branch
`instrument/cism-dumps`, not pushed — see `stage1/runs/CISM_LANE.md` §4d),
and `stage1/runs/cism_instrument.json`.

### The vertex stagger law — `cism_stagger_masked_mean`

`glissade_stagger` (glissade_grid_operators.F90:69-157) interpolates a
cell-centered field to the staggered grid.  The vertex `(i,j)` is the **NE
corner** of cell `(i,j)`, so its four cells are SW = `(i,j)`, SE = `(i+1,j)`,
NW = `(i,j+1)`, NE = `(i+1,j+1)` (F90:141-143).  With `stagger_margin = 1`
(F90:135-148) the vertex value is the **ice-masked mean**
`Σ mask·var / max(Σ mask, 1)` over the four cells, where the mask is the
`ice_plus_land_mask` of glissade_velo_higher.F90:1532-1537 (1 where
`ice_mask = 1 or land_mask = 1`); where all four masks are 0 the vertex takes
the 0 of the F90:123 initialization, which the `max(Σmask, 1)` denominator
reproduces exactly since the masks are 0/1.  This is the one-parameter
`cism_stagger_masked_mean` expression template (8 params: the 4 cell values
and the 4 masks).  `DivaIntegrals` applies it to `omega` (scalar) and per
velocity level to `omega_k` (the `stag_omega_k` array), using the cell's own
computed values as the SW cell and the dumped `omega_se/nw/ne` +
`omega_k_se/nw/ne` per-level values as the other cells.  The per-level
application is the exact transcription of the two stagger calls of
glissade_velo_higher.F90:2549-2552 (depth-integrated `omega`) and
F90:2561-2564 (per-level `omega_k`).

**Not staggered: θ and F.**  CISM computes the DIVA slope angles and
correction factors *directly at the vertices* from the plain vertex gradient
of `z_mean` (`glissade_slope_angle_staggered`, glissade_grid_operators.F90:
1482+, with the separate center/vertex calls of glissade_velo_higher.F90:
2174-2193 "to avoid issues with angle interpolation"), so `DivaGeometry`
remains the vertex law for `theta_x/y` and `F_x/y` and there is no stagger
template applied to them (a deviation from the task's literal wording,
documented here with the F90 citations).

### DivaGeometry — the vertex geometry of the driving stress

The vertex `(i,j)` is the **NE corner** of the traced cell; the stencil is the
SW/SE/NW/NE cells around it.  Inputs are the four cells' `usrf`, `thck`,
`ice_mask`, `land_mask` plus the run constants (`thklim`, `ramp`, `max_slope`,
`dx`, `rhoi`, `g`).  Laws (all in `glissade_grid_operators.F90` /
`glissade_velo_higher.F90`):

- `z_mean_* = usrf − thck/2` per stencil cell, and the **plain vertex
  gradient** `dz_dx = (z_mean_ne + z_mean_se − z_mean_sw − z_mean_nw)/(2·dx)`
  (and `dz_dy` with the n/s pairing, `glissade_slope_angle_staggered`'s
  argument, F90:1510).
- `dusrf_dx/dy` = the **hybrid-margin surface gradient**
  (`surface_elevation_gradient`, HO_GRADIENT_MARGIN_HYBRID, F90:999-1079):
  edge gradients from the interior pairs where both stencil cells are ice,
  else the `thklim`-guarded single-cell difference, then the vertex gradient
  from the four edge means, with the `max_slope` clamp (F90:1234-1252).  The
  `ramp` smoothing (50 m) applies only in the interior of the ice sheet where
  `thck` differences are large; it is off in the pinned rows by construction
  (checked: `dusrf_dx` matches the unramped hybrid formula on every exit
  row).
- `theta_x = atan(|dz_dx|)`, `theta_y = atan(|dz_dy|)` (slope_mask = 1
  everywhere on gris20), and the DIVA slope factor
  `F = cos²θ + 4·sin²θ` (F90:2189-2193).  These four are the vertex laws of
  the slope correction and are consumed by `DivaVelocity` through the
  `variable_map` edges `DivaGeometry.theta_x/y, F_x/y → DivaVelocity`.
- `stagthck = (thck_sw + thck_se + thck_nw + thck_ne)/4` (the 4-cell mean of
  `glissade_stagger`, `ice_plus_land_mask`).
- `taudx = −rhoi·g·stagthck·dusrf_dx`, `taudy = −rhoi·g·stagthck·dusrf_dy`
  (`load_vector_gravity`, HO_ASSEMBLE_TAUD_LOCAL, F90:4353-4384).
- `effecpress_stag = rhoi·g·stagthck` — the `bwat = 0` identity of the
  effective-pressure staggering (`HO_EFFECPRESS_BWAT_BVP` with `bwat = 0`
  everywhere in this 2-year run), pinned here at the same stencil rather than
  as an input.

### DivaIntegrals — the Goldberg (2011) σ-column integrals

`compute_integrals_diva` (F90:5570-5696), all inside the
`if (active_cell(i,j))` guard (F90:5091):

- `efvs` = the **arithmetic** mean of the four quadrature-point effective
  viscosities, `efvs[l] = (efvs_q1[l] + efvs_q2[l] + efvs_q3[l] + efvs_q4[l])/4`
  (F90:5266-5272 — a mean, not a harmonic mean).
- `omega_k[k]` = the interface-to-bed cumulative integral of
  `sig_mid[l]·thck·(1/4)·Σ_p 1/efvs_qp[l]` — written as **one ESD `integral`
  operator** over the σ measure (the integrand `omega_integrand` factors out
  the quad-point mean of reciprocals by integral linearity), lowered by the
  `integral_lev_nodes_cumulative_to_top` rule (the exact suffix scan of
  F90:5631-5637).  `omega_k[11] = 0` at the surface.
- `omega_layer_avg[l] = (omega_k[l] + omega_k[l+1])/2` (the trapezoid mean,
  F90:5656) and `omega` = the whole-column `integral` of `omega_layer_avg`
  over the σ measure (`integral_lev_whole`; the exact trapezoid sum of
  F90:5655-5661 — note **no `thck` factor**: this `omega` is Goldberg's
  ω/H).
- `stag_omega` and `stag_omega_k` = the **ice-masked 4-cell vertex means**
  (`cism_stagger_masked_mean`, `glissade_stagger` with `ice_plus_land_mask`
  and stagger_margin = 1, F90:135-148): `stag_omega` from the four
  cell-centered `omega` values (the cell's own as SW, the dumped
  `omega_se/nw/ne` as the others) and `stag_omega_k` applied **per velocity
  level** to the four cells' `omega_k` (the cell's own computed `omega_k` as
  SW, the dumped per-level `omega_k_se/nw/ne` as the others; the four
  `ice_plus_land_mask_sw/se/nw/ne` enter from the new `cism_velo_stag.csv`
  dump).  On every traced row of the reference run the four masks are 1, so
  the masked mean equals the plain 4-cell mean exactly; the masked branch is
  exercised by the marine-margin rows only if a partially grounded vertex
  exists (none does in this window — see Coverage).
- **Inactive-cell guard**: on ice-free cells the Fortran skips the whole
  block, so `efvs` and `omega_k` are 0 there (the dumped `efvs_*` columns
  confirm 0 while the `efvs_q*` columns carry **stale** previous-step values).
  The esm equivalent is a scalar-condition `ifelse(ice_mask == 1, …, 0)`,
  which is lazy in the spec (§4.8) — the integral branch is never evaluated
  on inactive cells, so the 1/efvs_qp blow-up is avoided; `omega_layer_avg`
  and `omega` need no guard (0 propagates).  The inactive pin is the step-3
  cell (17,14) test, where `thck` collapsed 97.3 → 0.986 m between steps.

### DivaVelocity — the 2-D recovery and the 3-D profile

`compute_3d_velocity_diva` (F90:5700-5821) and the calcbeta chain.  Every
input that is a *pointwise law of another component* now arrives through the
top-level `variable_map` edges (`param_to_var`, the same pattern as
`BasalTraction.beta → beta_raw`):

- `theta_x/y` and `F_x/y` ← `DivaGeometry` (the vertex slope-angle laws
  above).
- `omega` and `stag_omega_k` ← `DivaIntegrals` (the staggered integrals
  above).
- `beta_raw` ← the mounted `BasalTraction` (the pseudo-plastic law φ(bed)
  ramp → tanφ → τ_c = tanφ·N → β = τ_c/(u₀^q·speed^(1−q))) evaluated at the
  **accepted calcbeta pass** `(u, v, bed, N_stag)` of the Picard loop.

The laws:

- `beta_internal = f_ground·beta_raw` (the calcbeta weighting,
  F90:595-630).  The tests feed the accepted pass (the single
  `cism_basal_traction.csv` row whose `beta` equals the dumped
  `beta_internal`, velo_iter 38 of 39 on all 128 exit rows; its `bed`/`N`
  equal the exit row's `topg − eus`/`effecpress_stag` exactly), so the chain
  is **exact**: `beta_raw` reproduces the pass `beta` at rel 0.0 and
  `beta_internal = f_ground·beta_raw` the dumped value at rel 0.0 on all 128
  rows (the pass velocity iterate differs from the exit basal velocity by the
  final Picard correction — feeding the exit velocity instead would leave a
  ~6e-7 residual, the reason the earlier 1e-5 tolerance existed).
- `beta_eff_x = F_x·β/(F_x + β·ω)` and `beta_eff_y` analogously
  (the `cism_beta_eff` template, F90:5756-5770).
- **Basal recovery** (Goldberg 2011 eq. 34, F90:5801):
  `uvel_basal_recovered = uvel_2d·beta_eff_x/beta_internal`,
  `vvel_basal_recovered = vvel_2d·beta_eff_y/beta_internal` — pinned as
  explicit unknowns so the recovery law is tested rather than assumed.
- `btractx = β_internal·uvel_basal = beta_eff_x·uvel_2d`,
  `btracty` analogously (exit form F90:3874, solve form F90:5780); and the
  vertical profile `uvel[k] = uvel_basal + btractx·stag_omega_k[k]·cosθ_x/F_x`
  (F90:5805-5821) via the `slope_corr_x/y = cos(θ)/F` factors.
  `uvel_basal`/`vvel_basal` are **parameters** (the dumped exit
  `uvel_011`/`vvel_011`): the 3-D solve that produces them is the stage-3
  assembly.  Because `stag_omega_k_011 = 0` on every traced row, the pinned
  profiles satisfy `uvel[11] = uvel_basal` automatically — the same check the
  Fortran performs.

### Velo testing

15 tests / 155 assertions (manager review 2026-11-09 removed 20
lev_nodes = 11 pins the extractor still emits: uvel/vvel at the bed equal the
uvel_basal/vvel_basal overrides, and omega_k/stag_omega_k there are
identically 0).  Pinned rows (from the 128 exit rows of the 2026-11-09 dumps;
all tuples at full dump precision):

| test id | row | why this row |
|---|---|---|
| `geom_step1_cell2511` | step 1, cell (25,11) | the NE-corner stencil of the standard traced column (thin ice, 1.37 m) |
| `geom_step1_cell2713` | step 1, cell (27,13) | thick ice (434.5 m) with a negative surface slope in x |
| `geom_step3_cell2414` | step 3, cell (24,14) | the fastest traced column (uvel_2d = 5.02 m/yr) |
| `geom_step3_cell2813` | step 3, cell (28,13) | the steepest pinned slope (θ_x = 0.043) |
| `int_step1_cell2511` / `cell2612` / `cell2713` | step 1 | thin / mid (83.6 m) / thick columns: the σ-column integrals across the thickness range |
| `int_step3_cell2414` | step 3 | fastest column's integrals |
| `int_step3_cell1714_inactive` | step 3, cell (17,14) | **inactive-cell guard**: `ice_mask` dropped to 0 between steps (thck 97.3 → 0.986 m) with stale nonzero `efvs_q*` inputs |
| `velo_step1_cell2511` / `cell2713` / `step3_cell2414` | — | the coupled recovery + profile laws on the standard rows |
| `velo_step1_cell2213` | step 1, cell (22,13) | traced cell **ice-free** (thck = 0) with an active vertex — exercises the full 11-level profile through the β chain on a zero-thickness column |
| `marine_step1_cell2749` | step 1, cell (27,49) | the first marine-margin cell (topg = −6.8 m < eus): fully grounded, deep column with the largest `omega` of the pinned rows (2.36e-4 m/(Pa yr)) — the strongest β_eff correction of all pinned rows (β_eff/β_internal = 1.35e-3) |
| `marine_step3_cell3951` | step 3, cell (39,51) | the deep-bed marine cell (topg = −56.7 m) |

Tolerances: **rel 1e-12 everywhere** except `dusrf_dx/dy` + `taudx/taudy` at
rel 1e-11.  The `dusrf`/`taud` loosening is the dump round-trip of ~1800 m
elevation differences (worst deviation 6.6e-13 over the 128 rows; >15×
margin).  The β chain needs no loosening any more: feeding the accepted
calcbeta pass makes `beta_internal` exact to machine precision (rel 0.0 on
all 128 rows and all 32 marine rows), and every downstream law then holds at
≤ 2.9e-16 (worst profile deviation over all rows).

### Velo coverage limitations (this data window)

- **The 2-D finite-element matrix solve and the Picard iteration are not
  here** — they are the stage-3 assembly.  Consequently `uvel_2d`, `vvel_2d`,
  `uvel_basal`, `vvel_basal` are dump-fed inputs in this file
  (independent-model tests), and `effecpress_stag` is pinned via the bwat = 0
  identity instead of being computed.
- **The vertex stagger of `omega`/`omega_k` IS implemented** (the
  `cism_stagger_masked_mean` law in `DivaIntegrals`), but the four stencil
  cells' `omega`/`omega_k` enter as dumped inputs for the SE/NW/NE cells (the
  new `cism_velo_stag.csv` boundary): computing them from the neighbors'
  `efvs_q`/`thck` is the stage-3 field-level coupling.  The masks
  (`ice_plus_land_mask_sw/se/nw/ne`) are all 1 on every traced row, so the
  masked branch of the stagger law is verified only by construction (the
  `max(Σmask, 1)` guard) in this window.
- **`efvs_q*` provenance**: on inactive cells the dumped `efvs_q*` columns
  are stale (the last values computed while the cell was active) — the tests
  feed them as-is and the guard must zero the outputs, exactly as the
  Fortran's `active_cell` branch does.
- **`f_ground = 1` everywhere**: on all 128 exit rows and on all 32
  marine-margin exit rows of the reference run, so the
  `beta_internal = f_ground·beta_raw` weighting is pinned only in its
  saturated branch; the subgrid floatation weighting itself (F90:595-630) is
  not exercised by this data window.  The `cism_velo_marine_exit.csv` /
  `cism_velo_marine_pass.csv` boundary (the first `dump_marine_cells = 8`
  cells with initial ice over a marine bed, topg < eus) was added precisely
  to look for partially grounded vertices (f_ground < 1): none exist in this
  run, so the marine tests pin the fully grounded marine-bed configuration
  instead, and the f_ground < 1 branch stays documented, not exercised.
- The reference run's velocities are in **solver-internal m/yr** (the
  `scale_output` pass runs after the exit trace); all pinned numbers are in
  those units.

## Testing

All `tests` blocks are `time_span {0,0}` instantaneous pins with rel 1e-12
(the dumps are written `es26.17`; the recomputations land at ≤ 4.5e-16
relative, so 1e-12 is never stretched), except the transport ODE pins (the
solver block and the tolerance convention are documented in the transport
section) and the cap average-law pins that integrate one coupling interval.
Dump rows are cited by (file, cism_step, cell i,j, velo_iter/layer) in every
test description.

Run the suite with

    esm test components/cryosphere/cism

(`esm validate` and `esm units` pass for all eleven files; the two
discretization files and the velocity/thermo files that lower them need
`export ESD_ROOT=/path/to/EarthSciDiscretizations` for the `${ESD_ROOT}`
ref expansion — the pinned rules `integral_lev_whole` and
`integral_lev_nodes_cumulative_to_top` are both on ESD `main`). Suite
counts: **94 tests / 380 assertions** (basal 15/41, cap 22/48, flowlaw
10/19, thermo 10/88, transport 22/29, velo 15/155).

| File | Tests | Assertions | What is pinned |
|---|---|---|---|
| `cism_templates.esm` | — | — | template library only (validated, no models) |
| `eqice_basal_templates.esm` / `eqice_ice_rheology_templates.esm` | — | — | temporary EqIce copies (validated, no models; bodies identical to eqice/stage2) |
| `cism_basal.esm` | 15 | 41 | **Effecpress** (7 tests): the bwat=0 overburden branch at three dump cells/steps, the Bueler–van Pelt mid-range plateau (min(·,P) cap at work), full saturation, the s clamp, and the floating zero — the last four from the source formula; **BasalTraction** (8 tests): φ ramp + tanφ + speed + τ_c + β on dump rows covering both ramp-branch (φ = 34.33°, bed = 538 m) and phimax-branch (φ = 40°) cells, a 3-iteration Picard chain of one vertex (N/φ/τ_c invariant, β tracking the speed), the largest-N vertex (N = 1.79e6 Pa), a step-2 row, and the last iteration of the last step |
| `cism_flowlaw.esm` | 10 | 19 | warm-branch flwa at 6 (cell, layer, step) dump rows (incl. the largest \|pmptemp\| = 0.072 K traced), the cold branch, the branch boundary at T* = 263.15 K exactly, and both tempcor clamp bounds — the last four from the source formula |
| `cism_thermo.esm` | 10 | 88 | **IceTemperatureColumn**: the instantaneous operator identity at the exit state (5 pins, the 10 interior dTdt nodes each), temp_sfc/bpmp/temp_bed/lcondflx_ground/bmlt_ground/btemp_ground at the pinned cells, the thin-ice surface law (the ice_mask gate zeroing bpmp/temp_bed/bflx_mlt/btemp_ground), and the thawed-bed / internal-melt / freeze-on / eps11 branches from the source formula |
| `cism_transport.esm` | 22 | 29 | **GladSmb** (5 tests): the qsmb→acab chain at cap rows, the ice-free-ablation zeroing at thck = 0 and its non-firing at thck = 1.65e-5 m; **MassBalance** (14 tests): d(thck)/dt = net on accumulation/ablation/full-melt dump rows (rel 1e-6..3e-4, the advection residual), the melt_potential conservation check, acab_applied in the clipped branch, the layer-1 blend (rel 1e-5..2.5e-4), the tracer zeroing on full melt, and the ice-age increment; **SurfaceUpdate** (3 tests): lsrf/usrf on grounded, ice-free, and thin-ice rows |
| `cism_velo.esm` | 15 | 155 | **DivaGeometry** (4 tests): dusrf/taud at rel 1e-11, θ/F/effecpress_stag at rel 1e-12, the max_slope clamp from the source formula; **DivaIntegrals** (5 tests): efvs, omega_k, omega, the staggers and the inactive-cell guard; **DivaVelocity** (6 tests): the exact β chain at rel 0.0 and the 3-D profile at ≤ 2.9e-16 |
| `cism_cap.esm` | 22 | 48 | see the cap section's test table |

(`esm units --check` is clean on all eleven files; basal and velo report
the standard "no dimensional mismatch found; run without --check to see why
equations were not checked" message for their fractional-exponent rows,
unchanged from stage2.)

### Coverage limitations (this data window)

- **Cold flow-law branch**: all 1280 traced (cell, layer) rows have tempcor in
  [−5.16, −1.13] degC — the −50 floor, the 0 cap, and the tempcor < −10 cold
  branch never fire. Pinned synthetically from the source formula
  (`branch_cold`, `branch_pmp_exceeded`, `branch_floor`, `branch_boundary`).
- **bwat > 0 effective pressure**: empty trace (see How N is pinned).
- **φ = phimin branch**: no traced bed is below bedmin (traced beds span
  538..1786 m), so the ramp's lower plateau is covered only through the
  template structure (the ramp formula is pinned at both ends by the
  ramp/phimax branch tests).
- **Floating cells**: none in the window (all 256 velo rows floating_mask =
  0); the floating zeroing and the f_ground weighting are source-formula
  pins / documented no-ops.
- **Speed cap and β floor** (F90:196-198, 597-633): no-ops on the traced rows
  (beta_powerlaw_umax = 0; β > 100 Pa yr/m and f_ground = 1 everywhere);
  documented, not transcribed.
- **pmptemp column law** (F90:2170): the flowlaw rows carry pmptemp as input
  (the dump has no stagsigma column); the −ρᵢg·pmlt·thck·stagsigma law is
  part of the temperature-solve lane (the stage-1 validator's bpmp check
  pinned its bed value −ρᵢg·pmlt·thck in `cism_instrument.json`).
- The transport and cap coverage gaps are listed in their sections above.

## Next steps (stage 2.5 / stage 3)

1. Mount `Effecpress` + `BasalTraction` (+ `FlowFactor`) into the CISM
   subassembly with the cell→vertex stagger and the velocity solve; integrate
   with earthsciintegrator against the full-CESM trajectory
   (velo-exit uvel/vvel vs the solve, `btractx/btracty` = β·u·f_ground).
2. Couple `DivaGeometry` + `DivaIntegrals` + `DivaVelocity` (+ `BasalTraction`)
   into the top-level `eqesm.esm`: implement the 2-D DIVA matrix solve
   (`glissade_velo_higher_solve`'s Picard loop) and the `effecpress`
   staggering, then test against full-model simulations with
   `earthsciintegrator` time stepping.  The vertex staggers implemented here
   (`cism_stagger_masked_mean`) move to the field-level assembly by feeding the
   neighbor cells' computed integrals through the grid stencil instead of the
   dump columns.
3. Mount `IceTemperatureColumn` (+ the cap lanes) into the CISM subassembly:
   wire artm/thck/bwat/bfricflx/bhydroflx/bheatflx from the surface/velocity
   components, the dissip heating from the velocity solve, and the
   temperature state into the flow factor; the hydrology lane (row 10)
   supplies bwat. The basal-melt outputs (bmlt_ground, btemp_ground) feed the
   thickness/geometry evolution (not in this window).
4. A wetter stage-1 window (nonzero bwat / basal melt) to dump-pin the
   bwat > 0 effecpress branch, the cold flow-law branch, and the cap's
   nonzero flux averages.
5. When the EqIce project's components land in EarthSciModels, re-point the
   `eqice_*_templates.esm` imports at the EqIce files and delete the
   temporary copies.

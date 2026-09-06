# EarthSciModels repository layout

This repository holds authoritative `.esm` (EarthSciML Serialization Format) files
for MTK-derived Earth-science model components. It is a *data* repo: the files
here are the canonical serialized models, and CI validates that each one's inline
tests pass.

## Top-level layout

```
components/        # All .esm files, grouped by science domain
lib/               # Standard-library .esm subsystems (included from
                   # components via §4.7 reference; e.g. solar.esm)
docs/              # Migration tracker + layout convention
src/               # Julia shim (EarthSciModels.jl)
test/              # Shim tests + fixtures
.github/workflows/ # CI
```

Everything authoritative lives under `components/` and `lib/`. The shim, tests,
and CI exist to load and validate those files.

`lib/` is reserved for reusable, low-dependency subsystems that components
include rather than redefine — currently `lib/solar.esm` (NOAA Spencer-Fourier
solar declination, equation of time, zenith angle). Stdlib files use the same
`.esm` schema as components and are validated by the same inline-test
machinery; they are *not* organized by science domain because they cut across
domains.

## Why one flat `components/` tree (not per-schema-section dirs)

A `.esm` file can carry any mix of schema sections (`models:`,
`reaction_systems:`, `coupling:`, `data_sources:`) — see the ESM spec. Splitting the tree by schema section
(`models/`, `reaction_systems/`, ...) would force coupled components from one
upstream repo to scatter across directories and would make multi-section files
ambiguous to file. Organizing by **science domain** instead keeps related
content from one upstream repo together regardless of which schema sections
its `.esm` files use.

## Per-domain subdirs

Inside `components/`, each subdir corresponds to one **science domain**. For
now this is 1:1 with an upstream `earthsciml/*.jl` repo — domain and source
repo coincide. If a domain later spans multiple repos, or a repo crosses
domains, revisit this convention.

```
components/
  aerosol/                   # Aerosol.jl
  atmospheric_deposition/    # AtmosphericDeposition.jl
  atmospheric_dynamics/      # AtmosphericDynamics.jl
  earthsci_data/             # EarthSciData.jl (DataLoaders, ERA5/GEOSFP/...)
  earthsci_discretizations/  # EarthSciDiscretizations.jl (framework primitives)
  earthsci_ml_base/          # EarthSciMLBase.jl (framework primitives)
  environmental_transport/   # EnvironmentalTransport.jl
  gaschem/                   # GasChem.jl (atmospheric chemistry)
  geodynamics/               # Geodynamics.jl
  urban_canopy/              # UrbanCanopy.jl
  vegetation/                # Vegetation.jl
  wildland_fire/             # WildlandFire.jl
  atmospheric_radiation/     # radiation schemes transcribed from WRF (no upstream .jl repo yet; EqWeFiC)
```

Subdir names are lowercase snake_case forms of the upstream repo name (with
`.jl` dropped, hyphens and capitalization normalized).

## File naming

Within each per-domain subdir, `.esm` files are named in lowercase snake_case
matching the dominant component name (e.g. `superfast.esm`, `cloud_physics.esm`,
`era5.esm`). Multi-component files are named after the bundle they represent.

The migration tracker (`docs/migration-tracker.md`) lists `target_path`s as
suggestions for each component — Phase-3 migrators may rename or regroup as
the picture clarifies.

## One file per paper/chapter

**Rule (from `README.md`):** each `.esm` file holds approximately one paper or
chapter of content — a self-contained set of components that together describe
one scientific mechanism, with inline `tests` and `analyses`. This is *finer*
than
"one file per `src/<name>.jl`" in the upstream repo.

For example, `Aerosol.jl/src/aqueous_equilibria.jl` contains 9 `@component`
functions. These all describe the same aqueous-phase equilibria mechanism from
Seinfeld & Pandis, so they can live in a single `.esm` file
(`components/aerosol/aqueous_equilibria.esm`). In contrast,
`Aerosol.jl/src/cloud_physics.jl` covers 9 distinct mechanisms (water
properties, Kohler theory, droplet growth, ice physics, rain formation, ...)
and should be split into several `.esm` files within `components/aerosol/`.

## What is NOT in this repo

- MTK Julia source (lives in `Aerosol.jl`, `GasChem.jl`, etc. under
  github.com/EarthSciML).
- The ESM parser / schema / MTK-conversion code (lives in
  `EarthSciAST.jl`).
- Runtime simulation drivers (each model's inline `tests` + `analyses` are
  self-contained; end-to-end drivers go in the consumer's code).

## Empty-directory policy

`components/` and its per-domain subdirs are created on demand as real `.esm`
files land. There is no `.gitkeep` convention — empty trees are tracked via
the migration tracker, not via stub files.

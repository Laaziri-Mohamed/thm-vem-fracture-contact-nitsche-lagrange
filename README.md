# VEM Fracture Thermo-Poro‑Mechanics (Fortran)

Research code for **linear elasticity with fractures and contact** (Coulomb friction) coupled with **pressure/temperature (poro‑thermo) flow** on general polyhedral meshes using a **Virtual Element Method (VEM) of degree 1**. Two contact formulations are provided:

- **Nitsche** : Nitsche formulation.
- **Mixed VEM with bubble DOFs + Lagrange multipliers**: constant multipliers per fracture face, optional normal stiffness.

This README explains every file in this repository, how to build the solver, and how to run the example setups found in `main.f`.

> **Note**: The source is primarily commented in French. Variable names follow a consistent pattern: `Nb*` = numbers (counts), `by` = “per” (e.g., `NbNodebyCell`), `ff` = fracture face, `CG` = cell/face center of gravity (barycenter), `Inc` = DOF index, `VecTanFace1/2` = intrinsic tangential unit vectors on a face.

---

## Table of Contents

- [What’s in here?](#whats-in-here)
- [Numerical model at a glance](#numerical-model-at-a-glance)
- [Repository layout: file-by-file](#repository-layout-file-by-file)
- [Key data structures & conventions](#key-data-structures--conventions)
- [Build & dependencies](#build--dependencies)
- [Running the code](#running-the-code)
- [Configuration knobs (physics & numerics)](#configuration-knobs-physics--numerics)
- [Output files & visualization](#output-files--visualization)
- [Adding your own meshes](#adding-your-own-meshes)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Citations & references](#citations--references)

---

## What’s in here?

This code base targets 3D polyhedral grids and fractures represented as (internal) faces. The **mechanics** is solved with VEM (degree 1), optionally enriched with **bubble DOFs on fracture faces**. Contact is either **Nitsche**  or **mixed with constant per‑face multipliers** (Coulomb friction). The **porous flow + heat** part solves for **pressure (P)** and **temperature (T)** with Darcy/Fourier fluxes, including high‑order intersections (≥3 faces). A **fixed‑stress style relaxation** is applied in the matrix for the coupling.

---

## Numerical model at a glance

- **Elasticity**: Small‑strain, **linear isotropic** (Lamé `λ`, `μ`). Material parameters live in `Lois.f`.
- **Fractures**: Internal faces with jump operators; **intrinsic tangents** `VecTanFace1/2` and a **local normal** `VecNormalbyFace` (pointing outward of cell 1).
- **Contact**: 
  - **Coulomb friction** (angle set in `Lois.f` via `f_Frottement()`; default ≈25°). 
  - Optional **normal stiffness** `k_n` (via `f_UnSurkn()` → default `1/50 GPa^{-1}`) to allow controlled penetration.
  - **Active‑set** style update for stick/slip (see detailed comments in `main.f`).
- **VEM**: Degree 1 with options for **bubbles on fracture faces** to stabilize/contact handling in the mixed formulation.
- **Coupled P–T**: Darcy flow + Fourier heat conduction; assembly in `JacSm-PT.f` and solve/advance in `SolvePT.f`.
- **Coupling strategy**: `CouplagePoroMeca.f` updates **porosity** (and dφ/dP) and **fracture aperture** `df`, with a **fixed‑stress relaxation** applied in the rock matrix.

---

## Repository layout (file‑by‑file)

Below each file is summarized with its main entry points and responsibilities.

### `main.f` (Program entry point)

Top‑level driver that:
- Reads/sets **parameters** (time stepping, physics constants, switches).
- Loads a **mesh** (`.msh`) from several preconfigured locations (examples include `../../Vor/vmesh_*.msh`, `../../MAILLAGES/Cubic-Cells-transformed/gcube_*.msh`, etc.).
- Performs **DOF numbering** and geometry preprocessing (normals/tangents, gradient operators, jump operators).
- Chooses contact formulation with the **`iMecaVBulle`** flag:
  - `0` → **Nitsche** (no friction) → calls `ElasticiteNitscheFrac(...)`.
  - `1` → **VEM + bubbles + multipliers (with friction)** → calls `ElasticiteVEMBulleFrac(...)`.
- Advances **pressure/temperature** with `SolvePT(...)` and assembles their Jacobian/SM with `JacSmPT(...)` (indirectly).
- Manages **I/O**: writes `.val` diagnostic files and calls visualization hooks (e.g., Ensight/ParaView logic via `InitVisuParaview`, `VisuSolUP`, `FichierEnseightCase`).

**Notable parameters (see the `parameter (...)` lines near the top):**
- Display: `iaffich = 1` (enable visualization outputs).
- **Contact formulation**: `iMecaVBulle = 0` (Nitsche) or `1` (VEM bubbles + multipliers).
- Matrix **permeabilities**: `Permeabilitem1/m2/m3` (default `1e-15`, `1e-19`, `6e-13 m^2`).
- Initial **matrix porosity**: `Phiminit1/2/3` (defaults `0.1`, `0.01`, `0.1`).
- **Fluid viscosity**: `Viscof = 1e-3 Pa·s`.
- **Thermal conductivity**: `CondThermique = 2 W/(m·K)`.
- **Fracture width floor**: `dfmin = 1e-5 m`.
- **Time scales**: `UneHeure`, `UnJour`, `UnAn`; `TempsFinal`, `Deltat{{1..4}}`, `Deltatinit`, `DeltatMin`.
- Output cadence: `naffichpasdetemps`.

**Selected outputs created by `main.f`:**
- `dup_fixedpoint.val`, `dUtan.val`, `dunodes.val`, `SautntNitsche.val`
- `phimoydfmoy.val`, `pmmoy_pfmoy.val`, `dfpf.val`, `dt.val`
- `totalnewtonmeca.val`, `tempsaffich.val`

---

### `NitscheUtils.f` (Nitsche contact & geometrical utilities)

Utilities used across Nitsche and VEM formulations:

- `EvalXSautbyFaceFrac(...)` — **Pointwise jump** evaluation \((v^+ - v^-)(X)\) on a given fracture face using affine reconstructions (constant **tangential gradient per face** + nodal values). Outputs the local shape of the jump as a linear form of DOFs.
- `GradTangentbyFace(...)` — Builds **intrinsic tangential gradients** (2D subspace on the face) used for affine reconstructions and traction evaluation.
- `EvalXPThetaBetabyFaceFrac(...)` — Helper to evaluate **trial quantities** on fracture faces (e.g., jump components and duals), used by contact updates.
- `QuadFace(...)` — **Face quadrature** (weights/points) consistent with VEM/face operators.

These routines take the usual mesh topology arrays (`NbNodebyFace`, `NumNodebyFace`, `NumCellbyFace`, etc.), geometric locations (`XNode`, `XFaceCG`), and the **global/face DOF maps** (`NumIncGlobalbyFaceFrac`, `SautbyFaceFrac`).

---

### `NewtonMecaNitsche.f` (Mechanics: Nitsche )

- `ElasticiteNitscheFrac(...)` — Assembles and solves the **mechanical system with Nitsche terms** on fracture faces. No Coulomb friction here. Uses:
  - Mesh connectivity (`NbNodebyCell`, `NumNodebyFace`, `NumCellbyFace`, …),
  - Orientation (`VecNormalbyFace`, `VecTanFace1/2`),
  - VEM operators (`GradCellVEM`, tangential quantities, `SautbyFaceFrac`),
  - Thermo‑poro coupling through fields like `SolP`, `SolT` when present.
- `ResiduJacContact(...)` — Builds **contact residual/Jacobian** contributions specific to the Nitsche formulation (e.g., stabilization `β`, normal traction, jump consistency).

Use `iMecaVBulle=0` in `main.f` to activate this path.

---

### `NewtonMecaVBulle.f` (Mechanics: mixed VEM with bubbles & friction)

Implements the **mixed formulation** with **constant Lagrange multipliers per fracture face** and **optional bubble DOFs** on faces. Supports **Coulomb friction** and **normal stiffness**.

- `ElasticiteVEMBulleFrac(Temps, ...)` — Top‑level assembly/solve for mechanics in the mixed formulation. Takes the usual mesh & DOF maps, normals/tangents, VEM operators, and the current **P/T** fields for coupling.
- `ResiduContact(beta, rF, ...)` — Builds the **contact residual** (normal + tangential) using an **active‑set** approach. `rF` holds face‑wise friction data (e.g., coefficient).
- `JACContact(beta, rF, ...)` — Assembles the **contact Jacobian** consistent with `ResiduContact` to enable Newton convergence in stick/slip modes.

Use `iMecaVBulle=1` in `main.f` to activate this path.

---

### `Meca.f` (VEM building blocks for mechanics)

A toolbox for geometry, DOF numbering, and basic linear‑algebra utilities used by both mechanical formulations:

- `Calcul_LabelCellbyNode(...)` — Labels, per node, which cells are adjacent; used for local neighborhood logic.
- `NumIncMeca(...)` — **Mechanical DOF numbering** for cells/faces (bulk and fracture DOFs, face multipliers when present).
- `ComputegradvKnonplanar(...)`, `ComputegradvK(...)` — Compute **cellwise gradients** for VEM reconstructions; the `_nonplanar` variant handles non‑planar faces safely.
- `ComputepiK(...)` — VEM **projection operator** for degree‑1 fields on each cell.
- `StressFracMoy(...)` — Averages **face traction** / stress contributions on fracture faces.
- `ComputeSautbyFaceFrac(...)` — Builds the **jump operator** on fracture faces (normal & tangential components).
- `SigmanbyFaceFrac(...)` — Computes **normal traction** on fracture faces.
- Small helpers: `prodscal(...)` (scalar product), `prodtensor(...)` (tensor operation), `ProdMatVec(...)`, `prodscalMatrix(...)`.

---

### `SolvePT.f` (Advance pressure/temperature)

`SolvePT(...)` advances **pressure (P)** and **temperature (T)** with the current time step `Deltat`, assembling transmissibilities for:
- Bulk cells (`TransCell`, `TransCellFourier`),
- Fracture faces (`TransFaceFrac`, `TransFaceFracFourier`),
- **Matrix–fracture coupling** (`TransCellbyFaceFrac`, `TransCellbyFaceFracFourier`).

It uses **boundary condition incidence flags** (`IndIncDir`, `IndIncNeu`, and `...T` for temperature), face/cell DOF maps, and **aperture**/porosity arrays (`dfbyFaceFrac`, `PorobyCell`, and their `nm1` histories). Accumulation terms `AccP`, `AccT` and their previous states are also handled here.

---

### `JacSm-PT.f` (Jacobian & right‑hand side for P/T)

`JacSmPT(...)` assembles the **Jacobian matrix** and **second member** (negative residual) for P and T. It contains special treatment for **intersections with 3+ faces** (for thermal convection, it recenters Fourier + convective flux between face and edge).

Inputs include mesh topology around **interfaces** (`NbInterfacebyCell`, `NumInterfacebyCell`, etc.), fracture topology around **edges** (`NbFaceFracbyArete`, `NumFaceFracbyArete`), transmissibilities (as in `SolvePT`), porosity/aperture (current and `nm1`), and **relaxation factors** `relaxff` (fracture) & `relaxmm` (matrix).

---

### `CouplagePoroMeca.f` (Coupling: porosity & aperture)

`Compute_Porosity_Aperture(...)` updates:
- **Matrix porosity** per cell `PorobyCell` and its derivative **dφ/dP** `DerPorobyCell`,
- **Fracture aperture** per fracture face `dfbyFaceFrac`.

It includes a **fixed‑stress relaxation** *in the matrix only* (not in fractures) to stabilize P–U coupling.

---

### `Lois.f` (Material & physical laws)

Central location for physical parameters and simple constitutive relations:

- **Elasticity**: `f_Lame_mu()`, `f_Lame_lambda()` with defaults `E0=10 GPa`, `ν=0.25`.
- **Biot & storage**: `f_biot()`, `f_C0()` (and related `f_alpha0`, `f_alphaphi`).
- **Fracture/contact**: `f_UnSurkn()` (default `1 / 50 GPa`), `f_cohesion()` (if used), `f_Frottement()` (defaults to tan(25°)).
- **Gravity**: `f_gravite()` (≈10 m/s²).
- **Fluid & thermal props**: density/thermal capacity/enthalpy and derivatives (`f_rho`, `f_Cf`, `f_hf`, `f_ef`, and `dT`/`dp` variants), permeability & conductivity (`f_K0`, `f_Kf`, `f_UnSurKf`, `f_alphaf`), top heat flux `f_htop`, density of solid/fluid mix `f_rhosf`, etc.

> Tweak these functions to change your physical scenario without touching the solvers.

---

### `Makefile`

A minimal Makefile targeting **macOS (Accelerate framework)** with **SuperLU 4.3** and **gfortran/gcc (MacPorts)**:

- Variables:
  - `SuperLUroot`, `SUPERLULIB` — path to your SuperLU 4.3 static library.
  - `LIBS` — links SuperLU + `-framework Accelerate` (BLAS/LAPACK on macOS) + `-lm`.
  - `CC`/`CFLAGS` and `FORTRAN`/`FFLAGS` (note `-fallow-argument-mismatch` for legacy Fortran).
- Pattern rules for building `.o` from `.c` and `.f`.
- A `clean` target.

> If you build on Linux, replace `-framework Accelerate` with `-llapack -lblas` (or `-lopenblas`) and adjust compiler names/paths. You can also modernize to CMake if preferred.

---

## Key data structures & conventions

- **Dimensions & indexing** come from an `include 'include_parameter'` file:
  - Max sizes (e.g., `NbNodebyCellMax`, `NbFacebyCellMax`), space dim `NbDim` (3), etc.
- **Topology arrays** (typical shapes shown):
  - `NbNodebyCell(NbCell)`, `NumNodebyCell(NbCell, NbNodebyCellMax)`
  - `NbNodebyFace(NbFace)`, `NumNodebyFace(NbFace, NbNodebyFaceMax)`
  - `NbFacebyCell(NbCell)`, `NumFacebyCell(NbCell, NbFacebyCellMax)`
  - `NumCellbyFace(NbFace, 2)` — the two cells adjacent to each face
  - Fracture mappings: `NumFaceVersFaceFrac`, `NumFaceFracVersFace`
- **Geometry**:
  - Node coordinates `XNode(NbNode,3)`, centers `XCellCG`, `XFaceCG`
  - Areas/volumes: `SurfaceFace(NbFace)`, `VolCell(NbCell)`
  - Orthonormal frames on faces: `VecNormalbyFace`, `VecTanFace1/2`
    - **Normal** is oriented **outward of cell #1** of `NumCellbyFace`.
    - Tangents are **intrinsic**, so *for visualization* you may need to **flip** tangential components of jumps/multipliers using `OrientationbyFaceFrac` (see comments in `main.f`). 
- **DOF maps**:
  - `NumIncGlobalbyCell`, `NumIncGlobalbyFaceFrac` (and counts) to assemble global vectors/matrices.
  - For mixed contact, **constant multipliers per fracture face** are included.
- **Operators**:
  - `GradCellVEM` — constant cell gradient operator for VEM k=1
  - `GradTangbyFace` — constant tangential gradient per face
  - `SautbyFaceFrac` — jump operator on fracture faces (normal & tangential parts)

---

## Build & dependencies

### Requirements
- **Fortran compiler**: `gfortran` (tested with ≥ 11).
- **BLAS/LAPACK**:
  - macOS: Accelerate (`-framework Accelerate`).
  - Linux: `-llapack -lblas` or `-lopenblas`.
- **SuperLU 4.3** (static `libsuperlu_4.3.a`), headers in `${SuperLUroot}/SRC`.

### Quick build (Linux example)
```bash
# Edit paths as needed
export SuperLUroot=/path/to/SuperLU_4.3
export SUPERLULIB=$SuperLUroot/lib/libsuperlu_4.3.a

# Compile
gfortran -O2 -c main.f Meca.f NitscheUtils.f NewtonMecaNitsche.f \
               NewtonMecaVBulle.f SolvePT.f JacSm-PT.f CouplagePoroMeca.f Lois.f

# Link (with OpenBLAS)
gfortran -O2 -o frac *.o $SUPERLULIB -lopenblas -llapack -lm
```

### Quick build (macOS example using Accelerate)
```bash
export SuperLUroot=$HOME/WORK/PROGRAMMATION/SuperLU/SuperLU_4.3
export SUPERLULIB=$SuperLUroot/lib/libsuperlu_4.3.a

gfortran -O2 -fallow-argument-mismatch -c main.f Meca.f NitscheUtils.f \
    NewtonMecaNitsche.f NewtonMecaVBulle.f SolvePT.f JacSm-PT.f CouplagePoroMeca.f Lois.f

gfortran -O2 -o frac *.o $SUPERLULIB -framework Accelerate -lm
```

> If you get link errors: double‑check BLAS/LAPACK order and that `include_parameter` is present in your include path.

---

## Running the code

1. **Place a mesh** (Gmsh `.msh`) at one of the paths used in `main.f`, or edit the `open(unit=num,file=...)` lines to point to your mesh.
2. **Choose a contact formulation** in `main.f`:
   - `parameter ( iMecaVBulle = 0 )` → Nitsche (friction)
   - `parameter ( iMecaVBulle = 1 )` → Mixed LMs with friction & bubbles
3. (Optional) adjust **physics/time parameters** (permeabilities, porosities, final time `TempsFinal`, time steps `Deltat*`, etc.).
4. Run:
   ```bash
   ./frac
   ```
5. Inspect outputs (see next section).

---

## Configuration knobs (physics & numerics)

Tune these in `main.f` and `Lois.f`:

- **Elasticity**: set `E0`, `pois0` in `Lois.f` (`f_Lame_*`).
- **Contact**: friction via `f_Frottement()` (default `tan(25°)`), normal stiffness via `f_UnSurkn()` (default \(1/50\ \text{GPa}^{-1}\)).
- **Gravity**: `f_gravite()` (default ≈10).
- **Fluid**: viscosity `Viscof`, matrix permeabilities `Permeabilitem*`, fracture permeability `f_UnSurKf()`/`f_Kf()`.
- **Thermal**: `CondThermique`, heat capacity/enthalpy functions in `Lois.f`, top heat flux `f_htop()`.
- **Coupling**: fixed‑stress relaxation factors `relaxff`, `relaxmm` (used in `JacSm-PT.f`/`SolvePT.f` and `CouplagePoroMeca.f`).

---

## Output files & visualization

Common scalar outputs in the working directory (created by `main.f`):
- `dup_fixedpoint.val`, `dUtan.val`, `dunodes.val` — displacement increments (global/tangential/nodal).
- `SautntNitsche.val` — normal/tangential jump (Nitsche mode).
- `phimoydfmoy.val` — averaged porosity/aperture diagnostics.
- `pmmoy_pfmoy.val` — average matrix vs. fracture pressure.
- `dfpf.val` — fracture aperture vs. fracture pressure.
- `dt.val` — time step evolution.
- `totalnewtonmeca.val` — Newton iterations for mechanics.
- `tempsaffich.val` — wallclock timing around visualization steps.

**Visualization**: the code contains calls like `InitVisuParaview`, `VisuSolUP`, and `FichierEnseightCase` (Ensight). These generate files that you can load in **ParaView**. Check/adjust paths in `main.f` if needed.

---

## Adding your own meshes

- Meshes are opened via Fortran `open(... file='*.msh' ...)` calls in `main.f`. Provide your **Gmsh v2**‐style `.msh` files or adapt the reader.
- Make sure fracture faces are **tagged** appropriately so they populate `IndFaceFrac` and related mappings.
- Update the block around `imesh` in `main.f` to point to your file(s).

---


## License

Choose a license that matches your intent (e.g., **MIT**, **BSD‑3‑Clause**, **GPL‑3.0**). A common choice for research code is **MIT**:

```text
MIT License — Copyright (c) 2025-09-30
```

Add a `LICENSE` file at the repository root.

---

## Citations & references

If you use this code in publications, please cite relevant work on:

- **Virtual Element Methods (VEM) for linear elasticity**  
- **Contact mechanics with Coulomb friction** (active‑set / semismooth Newton)
- **Fixed‑stress splitting** for poromechanics

Laaziri, Mohamed, and Roland Masson. "VEM fully discrete Nitsche's discretization of Coulomb frictional contact-mechanics for mixed-dimensional poromechanical models." (2025).

Laaziri, Mohamed, and Roland Masson. "VEM-Nitsche fully discrete polytopal scheme for frictionless contact-mechanics." SIAM Journal on Numerical Analysis 63.1 (2025): 81-102.

---

## Appendix: quick map of major routines

- **Mechanics (Nitsche+ LM + bubbles + friction)**:  `ElasticiteNitscheFrac`, `ResiduJacContact`,`ElasticiteVEMBulleFrac`, `ResiduContact`, `JACContact`
- **P/T**: `SolvePT`, `JacSmPT`
- **Coupling**: `Compute_Porosity_Aperture`
- **VEM & geometry**: `ComputegradvK(_nonplanar)`, `ComputepiK`, `GradTangentbyFace`, `ComputeSautbyFaceFrac`
- **Face utilities**: `EvalXSautbyFaceFrac`, `EvalXPThetaBetabyFaceFrac`, `QuadFace`
- **Materials**: `f_Lame_mu`, `f_Lame_lambda`, `f_biot`, `f_UnSurkn`, `f_Frottement`, `f_K0/Kf`, `f_rho`, `f_Cf`, `f_hf`, `f_ef`, etc.

---

*Last updated: 2025-09-30*

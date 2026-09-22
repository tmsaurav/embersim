# embersim

A framework for coupled large-eddy simulation (LES) and Lagrangian particle
tracking, developed to simulate ember storms at the wildland-urban interface
(WUI). It couples the spectral-element solver [Nek5000](https://github.com/Nek5000/nek5000)
with the particle library [ppiclF](https://github.com/dpzwick/ppiclF), and adds
automated meshing, case templates, particle models and post-processing tools.

The framework was introduced in the PhD thesis cited below. The particle models
are **phenomenological**: their coefficients are calibrated, not derived from
experimental data (see "Particle models").

## Repository layout

```
embersim/
  case/          Case template: LES + particles; collects raw time-averaged
                 statistics during the run (KTH statistics module)
    ppiclf/      Particle setup: PPICLF_USER.h, ppiclf_user.f (forces, collision)
  case_stat/     Template for post-processing the raw statistics files (pstat3D)
                 into time-averaged fields, optionally on a structured grid
  presim/
    meshing/     Gmsh scripts: UniformBox.py, SimpleCubes.py, RefinedCubes.py
  postsim/
    modules/     flowpost.py, particlepost.py
    examples/    Jupyter notebooks (vtu -> npz, binned fields, hit maps)
  patches/       One-line hook applied to ppiclF at build time
  scripts/       env.sh, setup.sh, build_ppiclf.sh
  external/      Pinned dependencies (git submodules): Nek5000, KTH_Toolbox, ppiclF
```

Nothing inside `external/` is edited. Case-specific files live in `case/`, and the
single change to ppiclF is the patch in `patches/`.

## Requirements

- A Fortran/C toolchain with MPI. Built and tested with GNU Fortran/GCC 8.5 through
  Intel MPI 2021.17.2 on NCI Gadi. On other systems, edit `scripts/env.sh`.
- Git, and internet access for the setup and first build (see below).
- Python 3 with `gmsh`, `numpy`, `pyvista` and `matplotlib` for the meshing and
  post-processing tools.

## Getting started

Steps that download code (`setup.sh`, `build_ppiclf.sh` and the first compile fetch
gslib and other libraries) need internet access. On HPC systems, run them on a login
node and run simulations on compute nodes.

```bash
git clone https://github.com/tmsaurav/embersim.git
cd embersim
scripts/setup.sh              # fetches pinned dependencies, builds Nek5000 tools
source scripts/env.sh         # loads modules, sets paths (do this in every new shell)
```

`setup.sh` builds `genmap` and `gmsh2nek` into `external/Nek5000/bin`. If those
commands are not found, add that folder to your `PATH` or call them by full path.

### Workflow

1. **Mesh.** Edit the parameters at the top of a script in `presim/meshing/` and run it
   with Python. It writes a Gmsh `.msh` file with tagged boundaries.
2. **Convert.** Run `gmsh2nek` on the mesh (3 dimensions; give the periodic boundary
   IDs and translation vector). Rename the output to `case.re2` and run `genmap` with
   the name `case` to create `case.ma2`.
3. **Set up the case.** Copy `case/` to a working directory, add the mesh files, and
   edit `SIZE` (set `lelg`, `lelt` and related values for your mesh), `case.par` and
   `case.usr`. Domain-size numbers in `case.usr` (canopy limit, particle seeding
   region, particle walls and periodic ranges) are set for the thesis geometry and
   must be changed for a different domain.
4. **Build.** From the case directory:
   ```bash
   bash <embersim>/scripts/build_ppiclf.sh   # builds ppiclF with your ppiclf/ files + hook
   ./compile_script --all                    # builds the solver
   ```
5. **Run.**
   ```bash
   printf 'case\n%s/\n' "$(pwd)" > SESSION.NAME
   mpiexec -np <ranks> ./nek5000 > case.log.<ranks>
   ```
6. **Post-process.** Use `case_stat` for raw statistics files, and the notebooks in
   `postsim/examples/` for particle data.

`PPICLF_USER.h` sets particle array sizes for a case, so ppiclF is built per case
inside the case directory (`ppiclf_build/`), keeping the submodule clean.

## Particle models

Drag, lift, gravity and ground collision are defined in `case/ppiclf/ppiclf_user.f`.
**The code is the reference definition of the models.** Details of the thesis
(Chapter 5) describe the same approach, but the coefficients and exact expressions
in the code take precedence where they differ.

ppiclF deletes particles that leave the domain after a time step. A particle that
crosses the ground within one step would therefore vanish, so ground collision is
applied immediately after the Runge-Kutta update. `patches/ppiclf-post-rk3.patch`
adds a single call to `ppiclf_user_PostRK3` at that point in ppiclF, and the collision
routine itself sits in `ppiclf_user.f`. Without it, particles are lost through the
ground.

## Relationship to the thesis

The thesis describes an earlier layout in which Nek5000, ppiclF and the KTH Toolbox
were copied into the repository. This version replaces those copies with pinned
dependencies. The thesis-era state is preserved under the git tag `thesis-2025`
(archived on Zenodo: https://doi.org/10.5281/zenodo.22885431).

| Thesis path | This repository |
|---|---|
| `embersim/presim/meshing/` | `presim/meshing/` |
| `embersim/postsim/` | `postsim/` |
| `embersim/KTH_Framework/case`, `case_stat` | `case/`, `case_stat/` |
| `embersim/KTH_Framework/Nek5000` | `external/Nek5000` (pinned) |
| `embersim/KTH_Framework/Toolbox` | `external/KTH_Toolbox` (pinned) |
| `embersim/KTH_Framework/ppiclF` | `external/ppiclF` (pinned) plus `patches/` and `case/ppiclf/` |
| Setup steps of Section 3.2.3 | `scripts/setup.sh`, `scripts/build_ppiclf.sh`, `compile_script` |

Simulation results in the thesis were produced with the thesis-era version. This
version has been built and tested with the toolchain above, and its equivalence to
the thesis-era build has not been checked numerically.

## Not included

Example case inputs for the thesis simulations, and the precursor flow fields used to
restart them, are not part of this repository.

## Dependencies and licence

embersim is released under the MIT licence (see `LICENSE`). Third-party code keeps
its own licences: see `THIRD_PARTY.md` for the pinned versions, licences and
citations. If you use embersim, please also cite Nek5000, ppiclF and the KTH Toolbox.

## Citation

See `CITATION.cff`. Until a software release DOI is available for the maintained
version (v1.0.0 onwards), cite the thesis:

Saurav, T. M. (2025). *A Generalisable LES-Lagrangian Particle Framework for Ember
Storm Simulation at the Wildland-Urban Interface*. PhD thesis, UNSW Canberra.
https://doi.org/10.26190/unsworks/31916

The thesis-era snapshot of this repository (tag `thesis-2025`) is separately
archived and citable via Zenodo: https://doi.org/10.5281/zenodo.22885431

*Acknowledgement: The framework, models and simulations in this repository are the
author's original work, developed during his PhD. Claude (Anthropic) was used only
recently (Sep 2026), under the author's direction, to help restructure the repository into a
cleaner software package and to draft its documentation. The author reviewed and
takes responsibility for the result.*

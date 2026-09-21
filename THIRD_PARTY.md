# Third-party software

embersim does not contain copies of the codes below. They are fetched as git
submodules pinned to the exact commits listed, or downloaded during the build. Each
keeps its own licence. Please check the licence in each upstream repository before
redistributing.

| Component | Upstream | Pinned version | Licence |
|---|---|---|---|
| Nek5000 | https://github.com/Nek5000/nek5000 | commit `9b3d922` (`v19.0-37`) | see upstream `LICENSE` |
| KTH Toolbox | https://github.com/adampep/KTH_Toolbox | commit `b2b7a97` (the commit preceding upstream `6159834`) | see upstream `LICENSE` |
| ppiclF | https://github.com/dpzwick/ppiclF | commit `72483db` (`v1.0.0-43`) | MIT (see upstream `LICENSE`) |
| gslib | https://github.com/gslib/gslib | downloaded during the build by the Nek5000 and ppiclF install scripts | see upstream |
| Gmsh (Python API) | https://gmsh.info | installed with `pip` for the meshing scripts | GPL, see upstream |

Nek5000 also pulls further libraries into `external/Nek5000/3rd_party` during its own
build. Those are governed by their own licences.

## Changes to third-party code

embersim applies one change to ppiclF, at build time: `patches/ppiclf-post-rk3.patch`
adds a single call, `call ppiclf_user_PostRK3`, after the Runge-Kutta update in
`ppiclf_solve_IntegrateRK3`. The routine it calls is defined in this repository, in
`case/ppiclf/ppiclf_user.f`. No other upstream file is modified.

## Please cite

- Fischer, P. F., Lottes, J. W. and Kerkemeier, S. G. (2022). Nek5000: fast
  high-order scalable CFD. http://nek5000.mcs.anl.gov
- Massaro, D., Peplinski, A., Stanly, R., Mirzareza, S., Lupi, V., Mukha, T. and
  Schlatter, P. (2024). A comprehensive framework to enhance numerical simulations in
  the spectral-element code Nek5000. *Computer Physics Communications*, 302, 109249.
- Zwick, D. (2019). ppiclF: A Parallel Particle-In-Cell Library in Fortran. *Journal
  of Open Source Software*, 4(37), 1400.
- Geuzaine, C. and Remacle, J.-F. (2009). Gmsh: A 3-D finite element mesh generator
  with built-in pre- and post-processing facilities. *International Journal for
  Numerical Methods in Engineering*, 79(11), 1309-1331.

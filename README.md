# WUI Flow and Ember Storm

This repository contains the minimal setup for performing a flow simulation at a WUI. The folder named **embersim** contains the Nek5000 files necessary for running a case. Please note that it does not include a .box file, as the user is expected to supply the mesh using the provided Python script called `wuibuilder`. Current features and relevant dependencies are listed below:

## embersim

1. Forest canopy is implemented in `userf`, with the parameters set in `.par`.
2. A recycling boundary condition is set up using an old [example](https://github.com/Nek5000-deprecated/NekExamples-deprecated/blob/master/turbJet/jet.usr).
3. Flow is stabilised using an outflow condition from the turbChannel example, following the work of [Dong et al.](https://doi.org/10.1016/j.jcp.2013.12.042).
4. Currently all particle tracking capabilities are commented out. The user will need ppiclF installed and compiled to try this functionality.

## wuibuilder

This script generates a WUI with a canopy region and an APZ. This requires Gmsh SDK, which is available on Gmsh website. The generated mesh is in `.msh ` format, which can be converted into a Nek-supported `.re2 ` file using built-in `gmsh2nek` tool.

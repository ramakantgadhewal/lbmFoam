## lbmFoam

  DISCLAMER: This repository is under development.

# content

  - ./src/ contains source code for solvers lbm1Foam and lbm2Foam

  - ./test/ contains demo test cases

# usage

Basically:
  1. Compile source code using wmake in the source directory: this will create a new application in $FOAM_USER_APPBIN
      The latter can be created through 'mkdir -p $FOAM_USER_APPBIN'

  2. Test the application in the corresponding directory in the ./test/ folder

Some test case also contain a Makefile which allows a rapid execution:
  Inside the test case type

    make mesh | to create mesh
    make simulation | to run the solver
    make vtk | to eliminate time directories and obtain the results in VTK format

    make clean | to clean up everything
    make cleanMesh | to remove constant/polyMesh
    make cleanSim | to remove simulation logfile and eventually codeStream files
    make cleanVtk | to remove output files

    make test (or make) | to execute all in one shot

  Log file of each execution step will be saved in the log/ directory

Tested on OpenFOAM v9

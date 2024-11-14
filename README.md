# How to install an Actuation Disk Model
1. create a `src` dir into your `$WM_PROJECT_USER_DIR` and copy the source dir of the model you'd like to install to it. *You might need to do this in the OpenFOAM environment to be able to call the directories shortcuts*. For example:
```shell
mkdir -p $WM_PROJECT_USER_DIR/src
cd $WM_PROJECT_USER_DIR/src
cp -r ~/repos/actuationDiskRings/actuationDiskRingsV21_Source/actuationDiskRingsV21_Source .
```
2. compile the source with:
```shell
wmake libso .
```
3. check the installation by copying the `testCase` dir into your running dir and running it. 
    - test cases are set-up for 6 core computation. To modify this you should modify the variables `coresTotal` in `testCase/problem1Files/includeSettings1`.
    - To run it:
    ```shell
    bash RunFull.sh
    ```

# Models available:
- **actuationDiskRingsV21_Source**: combines the Uniform, Analytical, Elliptic, van der Laan Numerical and Navarro Diaz Numerical models into 1 model.
- **actuationDiskRingsV2_Source**: the Airfoil model.
- (*pending*) **actuationDiskRings_uniform**: the Uniform model.
- (*pending*) **actuationDiskRings_elliptic**: the Elliptic model.
- (*pending*) **actuationDiskRings_analytical**: the Analytical model.
- (*pending*) **actuationDiskRings_num_vanDerLaan**: the van der Laan Numerical model.
- (*pending*) **actuationDiskRings_num_navarroDiaz**: the Navarro Diaz Numerical model.

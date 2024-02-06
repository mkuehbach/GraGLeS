Compile structure_generator and twod_obc_solver using the following commands.

```
# assuming we are in the root directory of the project
export MYPROJECTHOME=$PWD
cd $MYPROJECTHOME/code/structure_generator/ && mkdir -p build && cd build
# assuming we are using the GNU compiler in a specific version (installed via classical APT from the GNU compiler apt repo)
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 .. 2>&1 | tee CMake.StructureGenerator.STDOUTERR.txt
make -j16 2>&1 | tee Make.StructureGenerator.STDOUTERR.txt
```

```
# assuming we are in the root directory of the project
export MYPROJECTHOME=$PWD
cd $MYPROJECTHOME/code/twod_obc_solver/ && mkdir -p build && cd build
# assuming we are using the GNU compiler in a specific version (installed via classical APT from the GNU compiler apt repo)
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 .. 2>&1 | tee CMake.TwodObcSolver.STDOUTERR.txt
make -j16 2>&1 | tee Make.TwodObcSolver.STDOUTERR.txt
```

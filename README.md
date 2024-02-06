How to install these software tools on a modern Linux (tested with Ubuntu):

```
# clone this repository
git clone https://github.com/mkuehbach/GraGLeS.git --branch mcl-modernization
cd GraGLeS
git submodule sync --recursive
git submodule update --init --recursive --jobs=4

# set-up and assure that you have a working system for compiling C/C++ applications:
# that is e.g. the GNU/C/C++ compiler, cmake, and (auto)make
# then follow the instructions in
# code/thirdparty/mandatory/fftw
# code/thirdparty/mandatory/hdf5
# code/thirdparty/mandatory/imkl

# compile structure_generator from within code/structure_generator/build via following the instructions in its ../README.md
# compile twod_obc_solver within code/twod_obc_solver/build via following the instructions in its ../README.md
# that is create a build directory and cd into them respectively, thereafter use cmake to create a Makefile
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 ..
# now build (using n many threads here n is 16)
make -j16
```

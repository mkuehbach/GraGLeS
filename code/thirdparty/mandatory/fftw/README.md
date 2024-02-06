The FFTW library is used as the fall-back solution to compile GraGLeS against
when using the Intel Math Kernel Library (IMKL) library is not desired.

The level-set solver routines compile either for single or double floating point
precision. The way how this has been implemented in the code is via linking against
the single and double precision versions of the FFTW library. Currently, this is
realized by compiling the FFTW library on the target system first using single
precision and thereafter using double precision. The two resulting static libraries
are linked via their thirdparty/mandatory/fftw/f32 and thirdparty/mandatory/fftw/f64
locations respectively.

This local compilation can be achieved with the following commands:

```
# assuming we are in the root directory of the project
export MYPROJECTHOME=$PWD
export MYFFTW_VERSION=3.3.10
cd $MYPROJECTHOME/code/thirdparty/mandatory/fftw
wget www.fftw.org/fftw-$MYFFTW_VERSION.tar.gz
mkdir -p f32 && cp fftw-$MYFFTW_VERSION.tar.gz f32
mkdir -p f64 && cp fftw-$MYFFTW_VERSION.tar.gz f64
rm fftw-$MYFFTW_VERSION.tar.gz

cd $MYPROJECTHOME/code/thirdparty/mandatory/fftw/f32
tar -xvf fftw-$MYFFTW_VERSION.tar.gz && rm fftw-$MYFFTW_VERSION.tar.gz cd fftw-$MYFFTW_VERSION
./configure --prefix=$PWD --enable-float --enable-avx512 2>&1 | tee FFTW.F32.Configure.STDOUTERR.txt
make -j16 2>&1 | tee FFTW.F32.Make.STDOUTERR.txt
make check 2>&1 | tee FFTW.F32.MakeCheck.STDOUTERR.txt
make install 2>&1 | tee FFTW.F32.Install.STDOUTERR.txt

cd $MYPROJECTHOME/code/thirdparty/mandatory/fftw/f64
tar -xvf fftw-$MYFFTW_VERSION.tar.gz && rm fftw-$MYFFTW_VERSION.tar.gz && cd fftw-$MYFFTW_VERSION
./configure --prefix=$PWD --enable-avx512 2>&1 | tee FFTW.F64.Configure.STDOUTERR.txt
make -j16 2>&1 | tee FFTW.F64.Make.STDOUTERR.txt
make check 2>&1 | tee FFTW.F64.MakeCheck.STDOUTERR.txt
make install 2>&1 | tee FFTW.F64.Install.STDOUTERR.txt
```
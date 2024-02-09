The HDF5 library is compiled from its source in the following way:

```
# assuming we are in the root directory of the project
export MYPROJECTHOME=$PWD
export MYHDF_VERSION_MAJOR="1"
export MYHDF_VERSION_MINOR="14"
export MYHDF_VERSION_RELEASE="2"
export MYHDF_VERSION="${MYHDF_VERSION_MAJOR}.${MYHDF_VERSION_MINOR}.${MYHDF_VERSION_RELEASE}"
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$MYHDF_VERSION_MAJOR.$MYHDF_VERSION_MINOR/hdf5-$MYHDF_VERSION/src/CMake-hdf5-$MYHDF_VERSION.tar.gz
# unfortunately oftentimes these release tar.gz archives are not proper tar.gz archives but double-packed
# in this case a two-stepped unpacking procedure is required like so
mkdir CMake-hdf5-$MYHDF_VERSION
gzip -d CMake-hdf5-$MYHDF_VERSION.tar.gz && tar -xvf CMake-hdf5-$MYHDF_VERSION.tar
cd CMake-hdf5-$MYHDF_VERSION
./build-unix.sh 2>&1 | tee HDF5.Build.STDOUTERR.txt
cd build && ./HDF5-$MYHDF_VERSION-Linux.sh --include-subdir --skip-license 2>&1 | tee ../HDF5.Install.STDOUTERR.txt
```
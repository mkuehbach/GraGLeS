How to configure, compile, and use jemalloc:

```
# assuming build-essential, m4, autoconf, and make surplus a compiler installed
# run the following code from within code/thirdparty/optional/jemalloc/jemalloc, i.e. in the cloned jemalloc submodule repo
mkdir -p local
# jemalloc will be installed in that local directory
autoconf
./configure --prefix=$PWD/local 2>&1 | tee JEMALLOC.Configure.STDOUTERR.txt
make -j16 2>&1 | tee JEMALLOC.Make.STDOUTERR.txt
make install 2>&1 | tee JEMALLOC.MakeInstall.STDOUTERR.txt
```


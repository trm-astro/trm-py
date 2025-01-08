brew install pcre
wget http://www.iausofa.org/2023_1011_C/sofa_c-20231011.tar.gz
wget ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot5.2.tar.gz

tar -xvzf sofa_c-20231011.tar.gz

brew install automake
brew install libtool

#This makes the libsofa_c.a we will need later
cd sofa/VERSION/c/src
make
# on mac
sudo mkdir /usr/local/lib/sofa
cp *.o *.a *.h /usr/local/lib/sofa/
# I think the sofa bit would be better placed in the usr/src directory
# on scrtpc
make test

# Make libcpgplot.a
# On mac, do this
http://mingus.as.arizona.edu/~bjw/software/pgplot_fix.html
# Otherwise follow the instructions in the README for PGPLOT
cd 
brew install pkg-config
pkg-config --libs --cflags /opt/homebrew/Cellar/pcre/*/lib/pkgconfig/libpcrecpp.pc 

# on mac, do this
export CPPFLAGS="-I/opt/homebrew/Cellar/pcre/8.45/include -I/usr/local/lib/pgplot -I/usr/local/lib/sofa"
export LDFLAGS="-L/opt/homebrew/Cellar/pcre/8.45/lib -L/usr/local/lib/pgplot -L/usr/local/lib/sofa"

./bootstrap
./configure

make

make check

make install

# On SCRTP, do this

module load GCCcore/11.3.0 GCC/11.3.0 Autotools/20220317 libtool/2.4.7 PCRE/8.45 PGPLOT/5.2.2 SOFA_C/19

# wget http://www.iausofa.org/2023_1011_C/sofa_c-20231011.tar.gz
# tar -xvzf sofa_c-20231011.tar.gz
# cd sofa/VERSION/c/src
# make
# make test

# on scrtpc, do this
export CPPFLAGS="-I$HOME/include/ -fPIC"
export LDFLAGS="-L$HOME/lib/ -fPIC"

./bootstrap

./configure -prefix="$HOME"

make

make check

make install


# new CMakeLists.txt version

# on mac
 conan install . --build=missing
 cd build
 conan profile detect --force
 # Note: PLPLOT_BUILD_TYPE= 0 or 4 are strongly preferred when plplot is installed to the system (specify path) or brew default. 2 fot non-default location. 
 cmake .. \
 -DCMAKE_TOOLCHAIN_FILE=/Users/pipgrylls/Code/trm-py-bundle/trm-py/src/cpp-subs/build/Release/generators/conan_toolchain.cmake \
 -DCMAKE_BUILD_TYPE=Release \
 -DPLPLOT_BUILD_TYPE=4
 #-DPLPLOT_USE_PATH=/usr/local/lib/pgplot #if you want to use a system install, non-standard path, or specify the path
 cmake --build .
 cd ..
 sudo cmake --install build
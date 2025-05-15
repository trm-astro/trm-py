#! /bin/bash

# define the cibw before build script for macosx

# This script is run before the build process starts

yum makecache
yum install -y plplot-devel plplot-libs

pip install conan
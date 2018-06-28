#!/bin/sh

set -e

if [ "x$PREFIX" = "x" ]; then
  PREFIX=`pwd`
fi

BUILD_DIR=build
BASEDIR=`pwd`/`dirname $0`

rm -rf "$BASEDIR/$BUILD_DIR"
mkdir -p "$BASEDIR/$BUILD_DIR"
set -e
cd "$BASEDIR/$BUILD_DIR"
cmake -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX="$PREFIX" $* "$BASEDIR"
make -j 8
make install
cd $PREFIX
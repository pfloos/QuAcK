#!/usr/bin/env python3
import os
import sys

if "QUACK_ROOT" not in os.environ:
   os.chdir("..")
   print("")
   print("Please set the QUACK_ROOT environment variable, for example:")
   print("")
   print("$ export QUACK_ROOT={0}".format(os.getcwd()))
   print("")
   sys.exit(1)

QUACK_ROOT=os.environ["QUACK_ROOT"]

compile_gfortran_mac = """
FC = gfortran
AR = libtool
FFLAGS = -I$IDIR -J$IDIR -Wall -Wno-unused -Wno-unused-dummy-argument -O3
CC = gcc
CXX = g++
LAPACK=-lblas -llapack
STDCXX=-lc++
"""

compile_gfortran_mac_debug = """
FC = gfortran -I$IDIR -J$IDIR
AR = libtool
FFLAGS = -I$IDIR -J$IDIR -Wall -g -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
CC = gcc
CXX = g++
LAPACK=-lblas -llapack
STDCXX=-lc++
"""


# TODO Change compiler here
# --------------------------

compiler = compile_gfortran_mac
#compiler = compile_gfortran_mac_debug

header = """#
# This file was automatically generated. Do not modify this file.
# To change compiling options, make the modifications in 
#  {0}/src/make_ninja.py
#

QUACK_ROOT={0}
IDIR=$QUACK_ROOT/include
LDIR=$QUACK_ROOT/lib
BDIR=$QUACK_ROOT/bin
SDIR=$QUACK_ROOT/src

""".format(QUACK_ROOT)

rule_fortran = """
rule fc
  command = $FC $FFLAGS -c $in -o $out && (mv -f *.mod $IDIR &> /dev/null || :) 

"""

rule_build_lib = """
rule build_lib
  command = $AR -static $in -o $out 
  description = Linking $out

"""

rule_build_exe = """

LIBS = $LDIR/libxcf90.a $LDIR/libxc.a $LDIR/libslatec.a $LDIR/libnumgrid.a $LAPACK $STDCXX

rule build_exe
  command = $FC $in $LIBS -o $out
  pool = console
  description = Linking $out

rule build_lib
  command = cd $dir ; ninja $out
  pool = console
  description = Compiling $out

"""

rule_slatec = """
rule make_slatec
  command = cd $QUACK_ROOT/slatec/src ; make static ; cp static/libslatec.a $LDIR
  pool = console
  description = Building Slatec

build $LDIR/libslatec.a: make_slatec
"""

rule_git_clone = """
rule git_clone
  command = cd $QUACK_ROOT ; git clone $url
  pool = console
  description = Cloning $in

"""

build_numgrid = """
rule make_numgrid
  command = cd $QUACK_ROOT/utils ; LDIR="$LDIR" SDIR="$SDIR" CC="$CC" CXX="$CXX" FC="$FC" ./install_numgrid.sh
  description = Building numgrid
  pool = console

build $LDIR/libnumgrid.a: make_numgrid 
"""
  

build_in_lib_dir = "\n".join([
	header,
	compiler,
	rule_fortran,
	rule_build_lib,
])

  
build_in_exe_dir = "\n".join([
	header,
	compiler,
	rule_fortran,
	rule_build_exe,
])

build_main = "\n".join([
	header,
        compiler,
        rule_git_clone,
        build_numgrid,
        rule_slatec,
])

exe_dirs = [ "QuAcK", "eDFT" ]
lib_dirs = list(filter(lambda x: os.path.isdir(x) and \
                                 x not in exe_dirs, os.listdir(".")))

def create_ninja_in_libdir(directory):
    def write_rule(f, source_file, replace):
        obj_file = os.path.join("obj", source_file.replace(replace, ".o"))
        f.write("build {0}: fc {1}\n".format(obj_file,source_file))
        return obj_file

    with open(os.path.join(directory, "build.ninja"),"w") as f:
        f.write(build_in_lib_dir)
        objects = []
        for filename in os.listdir(directory):
            for suffix in [".f", ".f90"]:
                if filename.endswith(suffix):
                   obj_file = write_rule(f, filename, suffix)
                   objects.append(obj_file)
        objects = " ".join(objects)
        f.write("build $LDIR/{0}.a: build_lib {1}\n".format(directory,objects))
        f.write("default $LDIR/{0}.a\n".format(directory))


def create_ninja_in_exedir(directory):
    def write_rule(f, source_file, replace):
        obj_file = os.path.join("obj", source_file.replace(replace, ".o"))
        f.write("build {0}: fc {1}\n".format(obj_file,source_file))
        return obj_file

    with open(os.path.join(directory, "build.ninja"),"w") as f:
        f.write(build_in_exe_dir)
        objects = []
        for filename in os.listdir(directory):
            for suffix in [".f", ".f90"]:
                if filename.endswith(suffix):
                   obj_file = write_rule(f, filename, suffix)
                   objects.append(obj_file)
        objects = " ".join(objects)
        for libname in lib_dirs:
           f.write("build $LDIR/{0}.a: build_lib\n  dir = $SDIR/{0}\n".format(libname))
        libs = " ".join([ "$LDIR/{0}.a".format(x) for x in lib_dirs])
        f.write("build $BDIR/{0}: build_exe {1} {2}\n".format(directory,libs,objects))
        f.write("default $BDIR/{0}\n".format(directory))


def create_main_ninja():

    libs = " ".join([ "$LDIR/{0}.a".format(x) for x in lib_dirs])
    with open("build.ninja","w") as f:
        f.write(build_main)
        f.write("""
rule build_exe
  command = cd $SDIR/$dir ; ninja $out
  pool = console

rule build_lib
  command = cd $dir ; ninja $out
  pool = console
  description = Compiling $out

""")
        for exe_dir in exe_dirs:
            f.write("build $BDIR/{0}: build_exe {1} $LDIR/libnumgrid.a $LDIR/libslatec.a\n".format(exe_dir,libs))
            f.write("  dir = {0} \n".format(exe_dir) )
        for libname in lib_dirs:
           f.write("build $LDIR/{0}.a: build_lib\n  dir = $SDIR/{0}\n".format(libname))
        f.write("build all: phony $BDIR/QuAcK $BDIR/eDFT\n")
        f.write("default all\n")

def create_makefile(directory):
   with open(os.path.join(directory, "Makefile"),"w") as f:
     f.write("""default:
	ninja
	make -C ..
""")

def main():
    for lib_dir in lib_dirs:
       create_ninja_in_libdir(lib_dir)
       create_makefile(lib_dir)

    for exe_dir in exe_dirs:
       create_ninja_in_exedir(exe_dir)
       create_makefile(lib_dir)

    create_main_ninja()

if __name__ == '__main__':
    main()

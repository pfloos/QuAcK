#!/usr/bin/env python3
import os
import sys

DEBUG=False
try:
  DEBUG = sys.argv[1] == "debug"
except:
  pass
 	

if "QUACK_ROOT" not in os.environ:
   os.chdir("..")
   print("")
   print("Please set the QUACK_ROOT environment variable, for example:")
   print("")
   print("$ export QUACK_ROOT={0}".format(os.getcwd()))
   print("")
   sys.exit(1)

QUACK_ROOT=os.environ["QUACK_ROOT"]

if not DEBUG:
    compile_gfortran_mac = """
FC = gfortran
AR = libtool -static -o
FFLAGS = -I$IDIR -J$IDIR -fbacktrace -g -Wall -Wno-unused-variable -Wno-unused -Wno-unused-dummy-argument -O3
CC = gcc
CXX = g++
LAPACK=-lblas -llapack
STDCXX=-lc++
FIX_ORDER_OF_LIBS=
"""

    compile_gfortran_linux = """
FC = gfortran
AR = ar crs
FFLAGS = -I$IDIR -J$IDIR -fbacktrace -g -Wall -Wno-unused -Wno-unused-dummy-argument -O3
CC = gcc
CXX = g++
LAPACK=-lblas -llapack
STDCXX=-lstdc++
FIX_ORDER_OF_LIBS=-Wl,--start-group 
"""
    
    compile_ifort_linux = """
FC = ifort -mkl=parallel -qopenmp
AR = ar crs
FFLAGS = -I$IDIR -g -Ofast -traceback
CC = icc
CXX = icpc
LAPACK=
STDCXX=-lstdc++
FIX_ORDER_OF_LIBS=-Wl,--start-group 
"""
else:
    compile_gfortran_mac = """
FC = gfortran
AR = libtool -static -o
FFLAGS = -I$IDIR -J$IDIR -fbacktrace -Wall -Wno-unused-variable -g -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
CC = gcc
CXX = g++
LAPACK=-lblas -llapack
STDCXX=-lc++
FIX_ORDER_OF_LIBS=
"""

    compile_gfortran_linux = """
FC = gfortran
AR = ar crs
FFLAGS = -I$IDIR -J$IDIR -fbacktrace -Wall -g -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
CC = gcc
CXX = g++
LAPACK=-lblas -llapack
STDCXX=-lstdc++
FIX_ORDER_OF_LIBS=-Wl,--start-group 
"""

compile_olympe = """
FC = ifort -mkl=parallel -qopenmp
AR = ar crs
FFLAGS = -I$IDIR -Ofast -traceback -xCORE-AVX512
CC = icc
CXX = icpc
LAPACK=
STDCXX=-lstdc++
FIX_ORDER_OF_LIBS=-Wl,--start-group 
"""

if sys.platform in ["linux", "linux2"]:
#   compiler = compile_gfortran_linux
   compiler = compile_ifort_linux 
#    compiler = compile_olympe
elif sys.platform == "darwin":
  compiler = compile_gfortran_mac
else:
  print("Unknown platform. Only Linux and Darwin are supported.")
  sys.exit(-1)

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

LIBXC_VERSION=5.0.0

""".format(QUACK_ROOT)

rule_fortran = """
rule fc
  command = $FC $FFLAGS -c $in -o $out

"""

rule_build_lib = """
rule build_lib
  command = $AR $out $in 
  description = Linking $out

"""
LIBS=""
rule_build_exe = """
LIBS = {0} $LAPACK $STDCXX

rule build_exe
  command = $FC $FIX_ORDER_OF_LIBS $in $LIBS -o $out
  pool = console
  description = Linking $out

rule build_lib
  command = cd $dir ; ninja $out
  pool = console
  description = Compiling $out

""".format(LIBS)

rule_git_clone = """
rule git_clone
  command = cd $QUACK_ROOT ; git clone $url
  pool = console
  description = Cloning $in

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
])

exe_dirs = [ "QuAcK"]
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
        libs = " ".join([ "$LDIR/{0}.a".format(x) for x in lib_dirs]) + " "+LIBS
        f.write("build $BDIR/{0}: build_exe {1} {2}\n".format(directory,libs,objects))
        f.write("default $BDIR/{0}\n".format(directory))


def create_main_ninja():

    libs = " ".join([ "$LDIR/{0}.a".format(x) for x in lib_dirs]) + " "+LIBS
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
            sources = [ "$SDIR/{0}/{1}".format(exe_dir,x) for x in  os.listdir(exe_dir) ]
            sources = filter(lambda x: x.endswith(".f") or x.endswith(".f90"), sources)
            sources = " ".join(sources)
            f.write("build $BDIR/{0}: build_exe {1} {2}\n".format(exe_dir,libs,sources))
            f.write("  dir = {0} \n".format(exe_dir) )

        for libname in lib_dirs:
            sources = [ "$SDIR/{0}/{1}".format(libname,x) for x in  os.listdir(libname) ]
            sources = filter(lambda x: x.endswith(".f") or x.endswith(".f90"), sources)
            sources = " ".join(sources)
            f.write("build $LDIR/{0}.a: build_lib {1} \n  dir = $SDIR/{0}\n".format(libname, sources))
        f.write("build all: phony $BDIR/QuAcK\n")
        f.write("default all\n")

def create_makefile(directory):
   with open(os.path.join(directory, "Makefile"),"w") as f:
     f.write("""default:
	ninja
	make -C ..

clean:
	ninja -t clean

debug:
	ninja -t clean
	make -C .. debug
""")

def main():
    for lib_dir in lib_dirs:
       create_ninja_in_libdir(lib_dir)
       create_makefile(lib_dir)

    for exe_dir in exe_dirs:
       create_ninja_in_exedir(exe_dir)
       create_makefile(exe_dir)

    create_main_ninja()

if __name__ == '__main__':
    main()

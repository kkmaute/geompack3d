#!/bin/sh

#cd bin
#gfortran -o f90split f90split.f90

#cd ..
#mkdir src

cd src

#rm -f *.f90

rm -f *.o 

#../bin/f90split ../geompack3.f90
#../bin/f90split ../i4lib.f90

files=`ls *f90`

for item in $files;do
  echo "... compiling $item"
  ifort $FFLAGS -c $item >& log
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $item
    cat log
    exit
  fi
  rm log
done

echo "... linking"
ar qc libgeompack3.icc.a *.o

rm -f *.o

mv libgeompack3.icc.a ../lib/.

cd ../examples

echo "... testing"
ifort -o geopack.icc.exe geompack3_prb.f90 ../lib/libgeompack3.icc.a >& log
if [ $? -ne 0 ]; then
    echo "Errors compiling geompack3_prb.f90"
    cat log
    exit
fi
rm log

geopack.icc.exe >& geopack.log

cd ..

rm -r -f ~/tpls/geompack3d
mkdir ~/tpls/geompack3d
mkdir ~/tpls/geompack3d/intel
mkdir ~/tpls/geompack3d/intel/lib
mv lib/libgeompack3.icc.a ~/tpls/geompack3d/lib/libgeompack.a


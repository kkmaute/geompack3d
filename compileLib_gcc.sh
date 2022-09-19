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
  gfortran -fPIC -O3 -c $item >& log
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $item
    cat log
    exit
  fi
  rm log
done

echo "... linking"
ar qc libgeompack3.a *.o

rm -f *.o

mv libgeompack3.a ../lib/libgeompack3.gcc.a

cd ../examples

echo "... testing"
gfortran -o geopack.exe geompack3_prb.f90 ../lib/libgeompack3.gcc.a >& log
if [ $? -ne 0 ]; then
    echo "Errors compiling geompack3_prb.f90"
    cat log
    exit
fi
rm log

geopack.exe >& geopack.log

cd ..

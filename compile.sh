gfortran -fPIC -O3 -c quadpack_double.f90
f2py -m reflectivity --opt=-O3 --f90flags='-fopenmp' -lgomp \
      -I. quadpack_double.o\
      -c Math.f90 Models.f90 Nanoparticle_Models.f90

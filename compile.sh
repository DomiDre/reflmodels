f2py -m reflectivity --opt=-O3 --f90flags='-fopenmp' -lgomp \
      -I. quadpack_double.o\
      -c MathFunctions.f90 ReflModels.f90 ReflNanoparticle.f90

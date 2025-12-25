mkdir -p ./builds
mkdir -p ./builds

gfortran -O3 -march=native -funroll-loops \
    ./include/mt19937.f90 ./include/mt19937_par.f90 ./include/fhash.f90 \
    ./include/net_loader.f90 ./include/epidemic.f90 \
    ./main.f90 \
    -o ./builds/main.out
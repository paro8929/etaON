g++ -I. -c linear_int.cpp
g++ -fpermissive -fopenmp -std=c++11 -I. shearv3.cpp linear_int.o -o shear

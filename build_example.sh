# Select one of the following 2 lines to generate the build files
cmake -DUSE_OPENMP=ON -DUSE_OPENACC=OFF -B build # for OpenMP
cmake -DUSE_OPENMP=OFF -DUSE_OPENACC=ON -B build # for OpenACC
# Build the code
cmake --build build/
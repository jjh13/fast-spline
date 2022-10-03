#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/wait.h>
#include <chrono>
#include <vector>
#include <tuple>

#define DIMENSION 3

extern "C" {
#if DIMENSION == 2
extern float __reconstruct__(float x, float y, float (*lkup)(int,int));
#elif DIMENSION == 3 && !defined(TEST_TRILINEAR)
extern float __reconstruct__(float x, float y, float z, float (*lkup)(int,int,int));
#elif DIMENSION == 4
extern float __reconstruct__(float x, float y, float z, float w, float (*lkup)(int,int,int,int));
#endif
}


float lattice_accessor_3d(int x, int y, int z) {
    if(x == y && y == z && x == 10)
        return 1.0;
    return 0.0;
}


int main(int argc, char const *argv[]) {
    //
    printf("%f\n", __reconstruct__(10., 9.5, 9.5, lattice_accessor_3d));
}

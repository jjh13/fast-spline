
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

typedef float float4 __attribute__ ((vector_size (16)));

#define DIMENSION 3

extern "C" {
#if DIMENSION == 2
extern float4 __reconstruct__(float4 x, float4 y, float (*lkup)(int,int));
#elif DIMENSION == 3 && !defined(TEST_TRILINEAR)
extern float4 __reconstruct__(float4 x, float4 y, float4 z, float (*lkup)(int,int,int));
#elif DIMENSION == 4
extern float4 __reconstruct__(float4 x, float4 y, float4 z, float4 w, float (*lkup)(int,int,int,int));
#endif
}


float lattice_accessor_3d(int x, int y, int z) {
    static int cnt = 0;
    if(x == y && y == z && x == 10)
        return 1.0;
    return 0.0;
}


int main(int argc, char const *argv[]) {
    //
    float4 x = {10.0, 9.5, 9.25, 9.0};
    float4 y = {10.0, 9.25, 9.25, 9.0};
    float4 z = {10.0, 9.5, 9.25, 9.0};
    float4 res =  __reconstruct__(x, y, z, lattice_accessor_3d);

    printf("%f\n", res[0]);
    printf("%f\n", res[1]);
    printf("%f\n", res[2]);
    printf("%f\n", res[3]);
}

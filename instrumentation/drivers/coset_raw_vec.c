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

struct _coset_info {
    uint32_t *bounds;
    float *buffer;
};

struct _lattice_info {
    uint32_t coset_count;
    uint32_t dimension;
    _coset_info *cosets;
};


typedef float float4 __attribute__ ((vector_size (16)));

extern "C" {
#if DIMENSION == 2
extern float4 __reconstruct__(float4 x, float4 y, _lattice_info *);
#elif DIMENSION == 3 && !defined(TEST_TRILINEAR)
extern float4 __reconstruct__(float4 x, float4 y, float4 z, _lattice_info *);
#elif DIMENSION == 4
extern float4 __reconstruct__(float4 x, float4 y, float4 z, float4 w, _lattice_info *);
#endif
}


float lattice_accessor_3d(int x, int y, int z) {
    if(x == y && y == z && x == 10)
        return 1.0;
    return 0.0;
}


int main(int argc, char const *argv[]) {
    float *coset0 = (float *)calloc(100*100*100, 4);
    float *coset1 = (float *)calloc(100*100*100, 4);

    coset0[(10*100 + 10)*100 + 10] = 1.0;

    // Build simple coset structure
    uint32_t bounds[3] = {100,100,100};
    _coset_info ci[2] = {{bounds, coset0}, {bounds, coset1}};
    _lattice_info li{1, 3, &ci[0]};

    float4 x = {10.0, 9.5, 9.25, 9.0};
    float4 y = {10.0, 9.25, 9.25, 9.0};
    float4 z = {10.0, 9.5, 9.25, 9.0};
    float4 res =  __reconstruct__(x, y, z, &li);

    printf("%f\n", res[0]);
    printf("%f\n", res[1]);
    printf("%f\n", res[2]);
    printf("%f\n", res[3]);
}

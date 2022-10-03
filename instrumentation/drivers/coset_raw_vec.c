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



extern "C" {
#if DIMENSION == 2
extern float __reconstruct__(float x, float y, _lattice_info *);
#elif DIMENSION == 3 && !defined(TEST_TRILINEAR)
extern float __reconstruct__(float x, float y, float z, _lattice_info *);
#elif DIMENSION == 4
extern float __reconstruct__(float x, float y, float z, float w, _lattice_info *);
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

    printf("%f\n", __reconstruct__(10., 9.5, 9.5, &li));
}

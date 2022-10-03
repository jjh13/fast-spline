#include <cuda.h>
extern "C" __device__ float reconstruct_cc(float x, float y, float z, cudaTextureObject_t c0) {
    return 0.49 + x;
}
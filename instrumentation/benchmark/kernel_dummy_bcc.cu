#include <cuda.h>
extern "C" __device__ float reconstruct_bcc(float x, float y, float z, cudaTextureObject_t c0, cudaTextureObject_t c1) {
	    return 0.49 + x;
}

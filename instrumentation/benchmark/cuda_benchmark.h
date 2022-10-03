#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "helper_cuda.h"

#define DIMENSION 3
#define EXP_SEED 420
typedef float float_type;
typedef float VolumeType;

#if defined(LATTICE_cc) || defined(LATTICE_CC)
    #define CARTESIAN_CUBIC
#elif defined(LATTICE_bcc) || defined(LATTICE_BCC)
    #define BODY_CENTERED_CUBIC
#elif defined(LATTICE_fcc) || defined(LATTICE_FCC)
    #define FACE_CENTERED_CUBIC
#endif

#if defined(CARTESIAN_CUBIC)
extern "C" __device__ float reconstruct_cc(float x, float y, float z, cudaTextureObject_t c0);

__device__ float reconstruct_bcc(float x,
                                 float y,
                                 float z,
                                 cudaTextureObject_t c0,
                                 cudaTextureObject_t c1) { return 0.0; }

__device__ float reconstruct_fcc(float x,
                                 float y,
                                 float z,
                                 cudaTextureObject_t c0,
                                 cudaTextureObject_t c1,
                                 cudaTextureObject_t c2,
                                 cudaTextureObject_t c3) { return 0.0; }
#elif defined(BODY_CENTERED_CUBIC)

__device__ float reconstruct_cc(float x, float y, float z, cudaTextureObject_t c0) { return 0.0; }
extern "C"  __device__ float reconstruct_bcc(float x,
                                 float y,
                                 float z,
                                 cudaTextureObject_t c0,
                                 cudaTextureObject_t c1);

__device__ float reconstruct_fcc(float x,
                                 float y,
                                 float z,
                                 cudaTextureObject_t c0,
                                 cudaTextureObject_t c1,
                                 cudaTextureObject_t c2,
                                 cudaTextureObject_t c3) { return 0.0; }

#elif defined(FACE_CENTERED_CUBIC)
__device__ float reconstruct_cc(float x, float y, float z, cudaTextureObject_t c0) { return 0.0; }
__device__ float reconstruct_bcc(float x,
                                 float y,
                                 float z,
                                 cudaTextureObject_t c0,
                                 cudaTextureObject_t c1) { return 0.0; }
extern "C" __device__ float reconstruct_fcc(float x,
                                 float y,
                                 float z,
                                 cudaTextureObject_t c0,
                                 cudaTextureObject_t c1,
                                 cudaTextureObject_t c2,
                                 cudaTextureObject_t c3);

#endif


#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <chrono>

#include "cuda_benchmark.h"

#define N_ITERS 1024*8
#define PROBLEM_SIZE 2048*64
#define BLOCK_SIZE 256
#define RECON_BOUND 128

#if defined(TEST_LINEAR_FETCH)
#define FETCH_LINEAR true
#else
#define FETCH_LINEAR false
#endif


 __global__  void tex_read(
                float *noopt,
                int offset,
                cudaTextureObject_t c0,
                cudaTextureObject_t c1,
                cudaTextureObject_t c2,
                cudaTextureObject_t c3)
{
    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    float rand_sum = 0.0;
    curandState state;
    curand_init(EXP_SEED+offset, tid, 0, &state);

    for(int iter=0; iter < N_ITERS; iter++){
        float x = curand_uniform(&state);
        float y = curand_uniform(&state);
        float z = curand_uniform(&state);

#if defined(CARTESIAN_CUBIC)
        rand_sum += reconstruct_cc(x*RECON_BOUND, y*RECON_BOUND,z*RECON_BOUND,  c0);
#elif defined(BODY_CENTERED_CUBIC)
        rand_sum += reconstruct_bcc(x*RECON_BOUND, y*RECON_BOUND,z*RECON_BOUND,  c0, c1);
#elif defined(FACE_CENTERED_CUBIC)
        rand_sum += reconstruct_fcc(x*RECON_BOUND, y*RECON_BOUND,z*RECON_BOUND,  c0, c1, c2, c3);
#endif //202 BCC, 128 CC, 214 D4S, 128 C4, 161 FCC


    }
    *noopt += rand_sum;
}


cudaArray *d_cosetArray0 = 0, *d_cosetArray1 = 0, *d_cosetArray2 = 0, *d_cosetArray3 = 0;
cudaTextureObject_t texCoset0, texCoset1, texCoset2, texCoset3;
cudaExtent sizeCoset0, sizeCoset1, sizeCoset2, sizeCoset3;

void setupVolumetricTexture(cudaTextureObject_t &texObject, cudaArray *d_volumeArray, bool linearFetch) {
    cudaResourceDesc texRes;
    cudaTextureDesc texDescr;

    /* Destroy current texture if it exists */
    if (texObject)  checkCudaErrors(cudaDestroyTextureObject(texObject));

    /* Setup the resource descriptor */
    memset(&texRes, 0, sizeof(cudaResourceDesc));
    texRes.resType = cudaResourceTypeArray;
    texRes.res.array.array =  d_volumeArray;


    /* We want un-normalized accesses */
    memset(&texDescr, 0, sizeof(cudaTextureDesc));
    texDescr.normalizedCoords = false;
    texDescr.filterMode = linearFetch ? cudaFilterModeLinear : cudaFilterModePoint;
    texDescr.addressMode[0] = cudaAddressModeWrap;
    texDescr.addressMode[1] = cudaAddressModeWrap;
    texDescr.addressMode[2] = cudaAddressModeWrap;
    texDescr.readMode = cudaReadModeElementType;
//    texDescr.readMode = cudaReadModeNormalizedFloat;

    /* Create the texture object */
    checkCudaErrors(cudaCreateTextureObject(&texObject, &texRes, &texDescr, NULL));
}


double time_experiment(
            float *noopt,
            int offset,
            const unsigned int &block_size,
            const cudaTextureObject_t &c0,
            const cudaTextureObject_t &c1,
            const cudaTextureObject_t &c2,
            const cudaTextureObject_t &c3) {
    auto const t0 = std::chrono::steady_clock::now();
    tex_read<<<PROBLEM_SIZE/block_size, block_size>>>(noopt, offset, c0, c1, c2, c3);
    cudaDeviceSynchronize();
    double const time0 = std::chrono::duration_cast<std::chrono::duration<double>>(
            std::chrono::steady_clock::now() - t0)
            .count();

    return ((double)(PROBLEM_SIZE) / time0)*double(N_ITERS);
}



VolumeType *allocate_random_vol(unsigned int x, unsigned int y, unsigned int z){
    VolumeType *ret = (VolumeType *)calloc(x*y, z*sizeof(VolumeType));
    int idx = 0;
    for(int i = 0; i < x; i++)
        for(int j = 0; j < y; j++)
            for(int k = 0; k < z; k++, idx++) {
                ret[idx] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            }
    return ret;
}

cudaArray *move_to_device(float *host, const cudaExtent &e) {
    cudaArray *coset_array;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<VolumeType>();
    checkCudaErrors(cudaMalloc3DArray(&coset_array, &channelDesc, e));

    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr(host, e.width * sizeof(VolumeType),
                          e.width, e.height);
    copyParams.dstArray = coset_array;
    copyParams.extent = e;
    copyParams.kind = cudaMemcpyHostToDevice;
    checkCudaErrors(cudaMemcpy3D(&copyParams));
    return coset_array;
}

//
void allocate_cc(unsigned int x, unsigned int y, unsigned int z) {
    sizeCoset0 = make_cudaExtent(x+1, y+1, z+1);
    VolumeType *h_coset0 = allocate_random_vol(x+1, y+1, z+1);
    d_cosetArray0 = move_to_device(h_coset0, sizeCoset0);
    setupVolumetricTexture(texCoset0, d_cosetArray0, FETCH_LINEAR);
    free(h_coset0);
}

void allocate_bcc(unsigned int x, unsigned int y, unsigned int z) {
	const unsigned int extend_x = x % 2;
	const unsigned int extend_y = y % 2;
	const unsigned int extend_z = z % 2;

    sizeCoset0 = make_cudaExtent((x/2) + 1, (y/2) + 1, (z/2) + 1);
    sizeCoset1 = make_cudaExtent((x/2) + extend_x, (y/2) + extend_y, (z/2) + extend_z);

    VolumeType *h_coset0 = allocate_random_vol((x/2) + 1, (y/2) + 1, (z/2) + 1);
    VolumeType *h_coset1 = allocate_random_vol((x/2) + extend_x, (y/2) + extend_y, (z/2) + extend_z);

    d_cosetArray0 = move_to_device(h_coset0, sizeCoset0);
    d_cosetArray1 = move_to_device(h_coset1, sizeCoset1);

    setupVolumetricTexture(texCoset0, d_cosetArray0, FETCH_LINEAR);
    setupVolumetricTexture(texCoset1, d_cosetArray1, FETCH_LINEAR);

    free(h_coset0);
    free(h_coset1);
}

void allocate_fcc(unsigned int x, unsigned int y, unsigned int z) {
	const unsigned int extend_x = x % 2;
	const unsigned int extend_y = y % 2;
	const unsigned int extend_z = z % 2;

    sizeCoset0 = make_cudaExtent((x/2) + 1, (y/2) + 1, (z/2) + 1);
    sizeCoset1 = make_cudaExtent((x/2) + extend_x, (y/2) + extend_y, (z/2));
    sizeCoset2 = make_cudaExtent((x/2) + extend_x, (y/2), (z/2) + extend_z);
    sizeCoset3 = make_cudaExtent((x/2) + extend_x, (y/2) + extend_y, (z/2));

    VolumeType *h_coset0 = allocate_random_vol((x/2) + 1, (y/2) + 1, (z/2) + 1);
    VolumeType *h_coset1 = allocate_random_vol((x/2) + extend_x, (y/2) + extend_y, (z/2));
    VolumeType *h_coset2 = allocate_random_vol((x/2) + extend_x, (y/2), (z/2) + extend_z);
    VolumeType *h_coset3 = allocate_random_vol((x/2) + extend_x, (y/2) + extend_y, (z/2));

    d_cosetArray0 = move_to_device(h_coset0, sizeCoset0);
    d_cosetArray1 = move_to_device(h_coset1, sizeCoset1);
    d_cosetArray2 = move_to_device(h_coset2, sizeCoset2);
    d_cosetArray3 = move_to_device(h_coset3, sizeCoset3);

    setupVolumetricTexture(texCoset0, d_cosetArray0, FETCH_LINEAR);
    setupVolumetricTexture(texCoset1, d_cosetArray1, FETCH_LINEAR);
    setupVolumetricTexture(texCoset2, d_cosetArray2, FETCH_LINEAR);
    setupVolumetricTexture(texCoset3, d_cosetArray3, FETCH_LINEAR);

    free(h_coset0);
    free(h_coset1);
    free(h_coset2);
    free(h_coset3);
}




int main() {
    float *noopt = 0;
//    setupVolumetricTexture(texCoset0, d_cosetArray0, false);
    cudaMalloc(&noopt, 1*sizeof(float));

#if defined(CARTESIAN_CUBIC)
    allocate_cc(RECON_BOUND, RECON_BOUND, RECON_BOUND);
#elif defined(BODY_CENTERED_CUBIC)
    allocate_bcc(RECON_BOUND, RECON_BOUND, RECON_BOUND);
#elif defined(FACE_CENTERED_CUBIC)
    allocate_fcc(RECON_BOUND, RECON_BOUND, RECON_BOUND);
#endif //202 BCC, 128 CC, 214 D4S, 128 C4, 161 FCC


    double mean_time = time_experiment(noopt, 0, 64, texCoset0, texCoset1, texCoset2, texCoset3);
    printf("%f,", mean_time*1e-6);
}
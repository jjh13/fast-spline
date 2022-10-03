#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "helper_cuda.h"

#define DIMENSION 3
typedef float float_type;
typedef float VolumeType;

#if defined(LATTICE_cc) || defined(LATTICE_CC)
    #define CARTESIAN_CUBIC
#elif defined(LATTICE_bcc) || defined(LATTICE_BCC)
    #define BODY_CENTERED_CUBIC
#elif defined(LATTICE_fcc) || defined(LATTICE_FCC)
    #define FACE_CENTERED_CUBIC
#endif

__device__ __global__ void tex_read(cudaTextureObject_t c0)
{
//    float help =  tex3D<float>(c0, .1,2.,3.);
    float help = 69.;
    printf("%f", help);
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

VolumeType *floatVol = 0;

int main() {
    sizeCoset0 = make_cudaExtent(32, 32, 32);
    sizeCoset1 = make_cudaExtent(32, 32, 32);
    sizeCoset2 = make_cudaExtent(32, 32, 32);
    sizeCoset3 = make_cudaExtent(32, 32, 32);

    floatVol = (VolumeType*) malloc(32*32*32*sizeof(VolumeType));
    for(int x = 0; x < 32*32*32; x++) floatVol[x] = 0.69;

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<VolumeType>();
    checkCudaErrors(cudaMalloc3DArray(&d_cosetArray0, &channelDesc, sizeCoset0));


    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr = make_cudaPitchedPtr(floatVol, sizeCoset0.width * sizeof(VolumeType),
                          sizeCoset0.width, sizeCoset0.height);
    copyParams.dstArray = d_cosetArray0;
    copyParams.extent = sizeCoset0;
    copyParams.kind = cudaMemcpyHostToDevice;
    checkCudaErrors(cudaMemcpy3D(&copyParams));

    setupVolumetricTexture(texCoset0, d_cosetArray0, false);
    printf("help world...");
    dim3 block(1, 1);
    dim3 grid(1, 1);

    tex_read<<<grid, block>>>( texCoset0);
}
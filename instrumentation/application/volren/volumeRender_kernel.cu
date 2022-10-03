/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

// Simple 3D volume renderer

#include <curand.h>
#include <curand_kernel.h>

#include "options.h"

#ifndef _VOLUMERENDER_KERNEL_CU_
#define _VOLUMERENDER_KERNEL_CU_

#include <helper_cuda.h>
#include <helper_math.h>

typedef unsigned int  uint;
typedef unsigned char uchar;


cudaArray *d_volumeArray  = 0;
#if defined(BCC_LATTICE) || defined(FCC_LATTICE)
cudaArray *d_volumeArray1 = 0;
#endif
#ifdef FCC_LATTICE 
cudaArray *d_volumeArray2 = 0;
cudaArray *d_volumeArray3 = 0;
#endif

cudaArray *d_transferFuncArray;

typedef unsigned char VolumeType;
//typedef unsigned short VolumeType;

extern "C" __device__ float reconstruct_cc(float x, float y, float z, cudaTextureObject_t c0);
extern "C" __device__ float reconstruct_bcc(float x, float y, float z, cudaTextureObject_t c0, cudaTextureObject_t c1);
extern "C" __device__ float reconstruct_fcc(float x, float y, float z, cudaTextureObject_t c0, cudaTextureObject_t c1, cudaTextureObject_t c2, cudaTextureObject_t c3);

cudaTextureObject_t	texObject; // For 3D texture

#if defined(BCC_LATTICE) || defined(FCC_LATTICE)
cudaTextureObject_t texObject1; // For 3D texture#endif
#endif
#ifdef FCC_LATTICE 
cudaTextureObject_t texObject2; // For 3D texture
cudaTextureObject_t texObject3; // For 3D texture
#endif

#if defined(TEST_ONE_CELL)

#ifdef CC_LATTICE
__device__ float tex_scale = 1;
#endif

#ifdef FCC_LATTICE
__device__ float tex_scale = 0.25;
#endif

#else
__device__ float tex_scale = 1;

#endif

cudaTextureObject_t transferTex; // For 1D transfer function texture

typedef struct
{
    float4 m[3];
} float3x4;

__constant__ float3x4 c_invViewMatrix;  // inverse view matrix

struct Ray
{
    float3 o;   // origin
    float3 d;   // direction
};

// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm

__device__
int intersectBox(Ray r, float3 boxmin, float3 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float3 invR = make_float3(1.0f) / r.d;
    float3 tbot = invR * (boxmin - r.o);
    float3 ttop = invR * (boxmax - r.o);

    // re-order intersections to find smallest and largest on each axis
    float3 tmin = fminf(ttop, tbot);
    float3 tmax = fmaxf(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
    float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

    *tnear = largest_tmin;
    *tfar = smallest_tmax;

    return smallest_tmax > largest_tmin;
}

// transform vector by matrix (no translation)
__device__
float3 mul(const float3x4 &M, const float3 &v)
{
    float3 r;
    r.x = dot(v, make_float3(M.m[0]));
    r.y = dot(v, make_float3(M.m[1]));
    r.z = dot(v, make_float3(M.m[2]));
    return r;
}

// transform vector by matrix with translation
__device__
float4 mul(const float3x4 &M, const float4 &v)
{
    float4 r;
    r.x = dot(v, M.m[0]);
    r.y = dot(v, M.m[1]);
    r.z = dot(v, M.m[2]);
    r.w = 1.0f;
    return r;
}

__device__ uint rgbaFloatToInt(float4 rgba)
{
    rgba.x = __saturatef(rgba.x);   // clamp to [0.0, 1.0]
    rgba.y = __saturatef(rgba.y);
    rgba.z = __saturatef(rgba.z);
    rgba.w = __saturatef(rgba.w);
    return (uint(rgba.w*255)<<24) | (uint(rgba.z*255)<<16) | (uint(rgba.y*255)<<8) | uint(rgba.x*255);
}

__global__ void
d_render(uint *d_output, uint imageW, uint imageH,
         float density, float brightness,
         float transferOffset, float transferScale, cudaTextureObject_t	tex,

         #if defined(BCC_LATTICE) || defined(FCC_LATTICE)
            cudaTextureObject_t tex1,
         #endif

        #ifdef FCC_LATTICE
            cudaTextureObject_t tex2,
            cudaTextureObject_t tex3,
        #endif

         cudaTextureObject_t	transferTex)
{
    const int maxSteps = 500;
    const float tstep = 0.01f;
    const float opacityThreshold = 0.95f;
    const float3 boxMin = make_float3(-1.0f, -1.0f, -1.0f);
    const float3 boxMax = make_float3(1.0f, 1.0f, 1.0f);

    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;

    if ((x >= imageW) || (y >= imageH)) return;

    float u = (x / (float) imageW)*2.0f-1.0f;
    float v = (y / (float) imageH)*2.0f-1.0f;

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
    eyeRay.d = normalize(make_float3(u, v, -2.0f));
    eyeRay.d = mul(c_invViewMatrix, eyeRay.d);

    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);

    if (!hit) return;

    if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back, accumulating color
    float4 sum = make_float4(0.0f);
    float t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;

    for (int i=0; i<maxSteps; i++)
    {
        // read from 3D texture
        // remap position to [0, 1] coordinates
        #ifdef CC_LATTICE
        float sample = reconstruct_cc(
            32.*(pos.x*0.45 + 0.5), 32.*(pos.y*0.45 + 0.5), 32.*(pos.z*0.45 + 0.5), tex); //  tex3D<float>(tex, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);
         #endif

        #ifdef BCC_LATTICE
        float sample = 4*reconstruct_bcc(
        50.*(pos.x*0.45 + 0.5), 50.*(pos.y*0.45 + 0.5), 50.*(pos.z*0.45 + 0.5), tex, tex1); //  tex3D<float>(tex, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);

        #endif

        #ifdef FCC_LATTICE
        float sample = 2*reconstruct_fcc(
            40*(pos.x*0.45 + 0.5), 40*(pos.y*0.45 + 0.5), 40*(pos.z*0.45 + 0.5), tex, tex1, tex2, tex3); //  tex3D<float>(tex, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);
  
        #endif    

        //sample *= 64.0f;    // scale for 10-bit data

        // lookup in transfer function texture
        float4 col = tex1D<float4>(transferTex, (sample-transferOffset)*transferScale);
        col.w *= density;

        // "under" operator for back-to-front blending
        //sum = lerp(sum, col, col.w);

        // pre-multiply alpha
        col.x *= col.w;
        col.y *= col.w;
        col.z *= col.w;
        // "over" operator for front-to-back blending
        sum = sum + col*(1.0f - sum.w);

        // exit early if opaque
        if (sum.w > opacityThreshold)
            break;

        t += tstep;

        if (t > tfar) break;

        pos += step;
    }

    sum *= brightness;

    // write output color
    d_output[y*imageW + x] = rgbaFloatToInt(sum);
}

extern "C"
void setTextureFilterMode(bool bLinearFilter)
{
    if (texObject){
        checkCudaErrors(cudaDestroyTextureObject(texObject));
    }
    #if defined(BCC_LATTICE) || defined(FCC_LATTICE)
    if (texObject1){
        checkCudaErrors(cudaDestroyTextureObject(texObject1));
    }
    #endif
    #ifdef FCC_LATTICE
    if (texObject2){
        checkCudaErrors(cudaDestroyTextureObject(texObject2));
    }
    if (texObject3){
        checkCudaErrors(cudaDestroyTextureObject(texObject3));
    }
    #endif
    cudaResourceDesc            texRes;
    memset(&texRes,0,sizeof(cudaResourceDesc));

    texRes.resType            = cudaResourceTypeArray;
    texRes.res.array.array    = d_volumeArray;

    cudaTextureDesc             texDescr;
    memset(&texDescr,0,sizeof(cudaTextureDesc));

    texDescr.normalizedCoords = false;
    texDescr.filterMode       = bLinearFilter ? cudaFilterModeLinear : cudaFilterModePoint;

    texDescr.addressMode[0] = cudaAddressModeWrap;
    texDescr.addressMode[1] = cudaAddressModeWrap;
    texDescr.addressMode[2] = cudaAddressModeWrap;

    texDescr.readMode = cudaReadModeNormalizedFloat;

    checkCudaErrors(cudaCreateTextureObject(&texObject, &texRes, &texDescr, NULL));

    #if defined(BCC_LATTICE) || defined(FCC_LATTICE)
    texRes.res.array.array    = d_volumeArray1;
    checkCudaErrors(cudaCreateTextureObject(&texObject1, &texRes, &texDescr, NULL));
    #endif

    #ifdef FCC_LATTICE

    texRes.res.array.array    = d_volumeArray2;
    checkCudaErrors(cudaCreateTextureObject(&texObject2, &texRes, &texDescr, NULL));

    texRes.res.array.array    = d_volumeArray3;
    checkCudaErrors(cudaCreateTextureObject(&texObject3, &texRes, &texDescr, NULL));
    #endif
}

extern "C"
void initCuda(void *h_volume, void *h_volume1, void *h_volume2, void *h_volume3, cudaExtent volumeSize)
{
    // create 3D array
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<VolumeType>();
    checkCudaErrors(cudaMalloc3DArray(&d_volumeArray, &channelDesc, volumeSize));
    #if defined(BCC_LATTICE) || defined(FCC_LATTICE)
    checkCudaErrors(cudaMalloc3DArray(&d_volumeArray1, &channelDesc, volumeSize));
    #endif
    #ifdef FCC_LATTICE
    checkCudaErrors(cudaMalloc3DArray(&d_volumeArray2, &channelDesc, volumeSize));
    checkCudaErrors(cudaMalloc3DArray(&d_volumeArray3, &channelDesc, volumeSize));
    #endif

    // copy data to 3D array
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    checkCudaErrors(cudaMemcpy3D(&copyParams));

    #if defined(BCC_LATTICE) || defined(FCC_LATTICE)
    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume1, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray1;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    checkCudaErrors(cudaMemcpy3D(&copyParams));
    #endif
    #ifdef FCC_LATTICE
    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume2, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray2;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    checkCudaErrors(cudaMemcpy3D(&copyParams));

    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume3, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray3;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    checkCudaErrors(cudaMemcpy3D(&copyParams));
    #endif


    cudaResourceDesc            texRes;
    memset(&texRes, 0, sizeof(cudaResourceDesc));

    texRes.resType            = cudaResourceTypeArray;
    texRes.res.array.array    = d_volumeArray;

    cudaTextureDesc             texDescr;
    memset(&texDescr, 0, sizeof(cudaTextureDesc));

    texDescr.normalizedCoords = false; // access with normalized texture coordinates
    texDescr.filterMode       = LIN_FETCH ? cudaFilterModeLinear : cudaFilterModePoint;


    texDescr.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texDescr.addressMode[1] = cudaAddressModeClamp;
    texDescr.addressMode[2] = cudaAddressModeClamp;

    texDescr.readMode = cudaReadModeNormalizedFloat;

    checkCudaErrors(cudaCreateTextureObject(&texObject, &texRes, &texDescr, NULL));

    #if defined(BCC_LATTICE) || defined(FCC_LATTICE)

    memset(&texRes, 0, sizeof(cudaResourceDesc));
    texRes.resType            = cudaResourceTypeArray;
    texRes.res.array.array    = d_volumeArray1;
    checkCudaErrors(cudaCreateTextureObject(&texObject1, &texRes, &texDescr, NULL));
    #endif

    #ifdef FCC_LATTICE

    texRes.res.array.array    = d_volumeArray2;
    checkCudaErrors(cudaCreateTextureObject(&texObject2, &texRes, &texDescr, NULL));

    texRes.res.array.array    = d_volumeArray3;
    checkCudaErrors(cudaCreateTextureObject(&texObject3, &texRes, &texDescr, NULL));
    #endif

    // create transfer function texture
    float4 transferFunc[] =
    {
        {  0.0, 0.0, 0.0, 0.0, },
        {  1.0, 0.0, 0.0, 1.0, },
        {  1.0, 0.5, 0.0, 1.0, },
        {  1.0, 1.0, 0.0, 1.0, },
        {  0.0, 1.0, 0.0, 1.0, },
        {  0.0, 1.0, 1.0, 1.0, },
        {  0.0, 0.0, 1.0, 1.0, },
        {  1.0, 0.0, 1.0, 1.0, },
        {  0.0, 0.0, 0.0, 0.0, },
    };

    cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();
    cudaArray *d_transferFuncArray;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArray, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArray, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    memset(&texRes,0,sizeof(cudaResourceDesc));

    texRes.resType            = cudaResourceTypeArray;
    texRes.res.array.array    = d_transferFuncArray;

    memset(&texDescr,0,sizeof(cudaTextureDesc));

    texDescr.normalizedCoords = true; // access with normalized texture coordinates
    texDescr.filterMode       = cudaFilterModeLinear;

    texDescr.addressMode[0] = cudaAddressModeClamp; // wrap texture coordinates

    texDescr.readMode = cudaReadModeElementType;

    checkCudaErrors(cudaCreateTextureObject(&transferTex, &texRes, &texDescr, NULL));
}

extern "C"
void freeCudaBuffers()
{
    checkCudaErrors(cudaDestroyTextureObject(texObject));
    checkCudaErrors(cudaDestroyTextureObject(transferTex));
    checkCudaErrors(cudaFreeArray(d_volumeArray));
    checkCudaErrors(cudaFreeArray(d_transferFuncArray));
}


extern "C"
void render_kernel(dim3 gridSize, dim3 blockSize, uint *d_output, uint imageW, uint imageH,
                   float density, float brightness, float transferOffset, float transferScale)
{
    d_render<<<gridSize, blockSize>>>(d_output, imageW, imageH, density,
                                      brightness, transferOffset, transferScale,
                                      texObject, 
                                      #if defined(FCC_LATTICE) || defined(BCC_LATTICE)
                                      texObject1,
                                      #endif

                                      #ifdef FCC_LATTICE
                                      texObject2,
                                      texObject3,

                                      #endif
                                      transferTex);
}

// 
__device__ double d_noopt = 0;

__global__ void d_test(
            cudaTextureObject_t tex0,
         #if defined(BCC_LATTICE) || defined(FCC_LATTICE)
            cudaTextureObject_t tex1,
         #endif

        #ifdef FCC_LATTICE
            cudaTextureObject_t tex2,
            cudaTextureObject_t tex3,
        #endif

         unsigned int *dummy) {

    double x,y,z;

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    curandState state;
    curand_init(100, idx, 0, &state);
\
    double local = 0;
    for(int i = 0; i < 10000; i++) {
        x =  curand_uniform(&state);
        y =  curand_uniform(&state);
        z =  curand_uniform(&state);

        #ifdef CC_LATTICE
        float sample = reconstruct_cc(
            32.*(x*0.45 + 0.5), 
            32.*(y*0.45 + 0.5), 
            32.*(z*0.45 + 0.5), 
            tex0); //  tex3D<float>(tex, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);
         #endif

        #ifdef BCC_LATTICE
        float sample = reconstruct_bcc(
        50.*(x*0.45 + 0.5), 
        50.*(y*0.45 + 0.5), 
        50.*(z*0.45 + 0.5), 
        tex0, 
        tex1); //  tex3D<float>(tex, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);

        #endif

        #ifdef FCC_LATTICE
        float sample = reconstruct_fcc(
            40*(x*0.45 + 0.5), 
            40*(y*0.45 + 0.5), 
            40*(z*0.45 + 0.5), 
            tex0, tex1, tex2, tex3); //  tex3D<float>(tex, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);
  
        #endif    
        local += sample;
    }
    d_noopt = local;
    // atomicAdd(dummy, 1);
}

extern "C"
void random_test(dim3 gridSize, dim3 blockSize){

    float time;
    cudaEvent_t start, stop;

    // d_test<<<gridSize, blockSize>>>(
    //                                   texObject, 
    //                                   #if defined(FCC_LATTICE) || defined(BCC_LATTICE)
    //                                   texObject1,
    //                                   #endif

    //                                   #ifdef FCC_LATTICE
    //                                   texObject2,
    //                                   texObject3,

    //                                   #endif
    //                                   0.0);


    // Create events and record the start time
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);


    unsigned int *d_data; cudaMalloc(&d_data, sizeof(unsigned int));    
    cudaMemset(d_data, 0, sizeof(unsigned int));

    // 
    d_test<<<gridSize, blockSize>>>(
                                      texObject, 
                                      #if defined(FCC_LATTICE) || defined(BCC_LATTICE)
                                      texObject1,
                                      #endif

                                      #ifdef FCC_LATTICE
                                      texObject2,
                                      texObject3,

                                      #endif
                                      d_data);

    unsigned int host_count = 0; 

    cudaMemcpy(&host_count,d_data,sizeof(unsigned int),cudaMemcpyDeviceToHost);
    // printf("%d\n", host_count);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);


    printf("mean:%e, var:%e \n", double(time)/double(1024.*32.*10000.*32.*1000.),0);
}

extern "C"
void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix)
{
    checkCudaErrors(cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeofMatrix));
}


#endif // #ifndef _VOLUMERENDER_KERNEL_CU_

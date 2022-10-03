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

/*
    Volume rendering sample

    This sample loads a 3D volume from disk and displays it using
    ray marching and 3D textures.

    Note - this is intended to be an example of using 3D textures
    in CUDA, not an optimized volume renderer.

    Changes
    sgg 22/3/2010
    - updated to use texture for display instead of glDrawPixels.
    - changed to render from front-to-back rather than back-to-front.
*/

// OpenGL Graphics includes
#include <helper_gl.h>
#if defined (__APPLE__) || defined(MACOSX)
  #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #include <GLUT/glut.h>
  #ifndef glutCloseFunc
  #define glutCloseFunc glutWMCloseFunc
  #endif
#else
#include <GL/freeglut.h>
#endif

#define IN_TEST 0

// CUDA Runtime, Interop, and includes
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <cuda_profiler_api.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>

// CUDA utilities
#include <helper_cuda.h>

// Helper functions
#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_timer.h>

typedef unsigned int uint;
typedef unsigned char uchar;

#define MAX_EPSILON_ERROR 5.00f
#define THRESHOLD         0.30f
#include "options.h"

// #define CC_LATTICE
// #define TEST_ONE_CELL

// Define the files that are to be save and the reference images for validation
const char *sOriginal[] =
{
    "volume.ppm",
    NULL
};

const char *sReference[] =
{
    "ref_volume.ppm",
    NULL
};

const char *sSDKsample = "CUDA 3D Volume Render";

const char *volumeFilename = "Bucky.raw";


#ifdef CC_LATTICE
cudaExtent volumeSize = make_cudaExtent(32, 32, 32);
#endif 

#ifdef BCC_LATTICE
cudaExtent volumeSize = make_cudaExtent(25, 25, 25);
#endif

#ifdef FCC_LATTICE
cudaExtent volumeSize = make_cudaExtent(20, 20, 20);
#endif

typedef unsigned char VolumeType;

//char *volumeFilename = "mrt16_angio.raw";
//cudaExtent volumeSize = make_cudaExtent(416, 512, 112);
//typedef unsigned short VolumeType;

uint width = 512, height = 512;
dim3 blockSize(16, 16);
dim3 gridSize;

float3 viewRotation;
float3 viewTranslation = make_float3(0.0, 0.0, -4.0f);
float invViewMatrix[12];

float density = 0.05f;
float brightness = 1.0f;
float transferOffset = 0.0f;
float transferScale = 1.0f;
bool linearFiltering = true;
bool exit_program = false;

GLuint pbo = 0;     // OpenGL pixel buffer object
GLuint tex = 0;     // OpenGL texture object
struct cudaGraphicsResource *cuda_pbo_resource; // CUDA Graphics Resource (to transfer PBO)

StopWatchInterface *timer = 0;
StopWatchInterface *timerInternal = 0;

// Auto-Verification Code
const int frameCheckNumber = 2;
int fpsCount = 0;        // FPS count for averaging
int fpsLimit = 1;        // FPS limit for sampling
int g_Index = 0;
unsigned int frameCount = 0;

int *pArgc;
char **pArgv;

#ifndef MAX
#define MAX(a,b) ((a > b) ? a : b)
#endif

extern "C" void setTextureFilterMode(bool bLinearFilter);
extern "C" void initCuda(void *h_volume, void *h_volume1, void *h_volume2,void *h_volume3,cudaExtent volumeSize);
extern "C" void freeCudaBuffers();
extern "C" void render_kernel(dim3 gridSize, dim3 blockSize, uint *d_output, uint imageW, uint imageH,
                              float density, float brightness, float transferOffset, float transferScale);

extern "C" void random_test(dim3 gridSize, dim3 blockSize);
extern "C" void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix);

void initPixelBuffer();

void sampleMarshnerLobb_CC(unsigned char *buffer, uint32_t x, uint32_t y, uint32_t z) {
    float dx = 1./float(x);
    float dy = 1./float(y);
    float dz = 1./float(z);

    const float f_m = 6.;
    const float alpha = .25;

    for(uint32_t i = 0; i < x; i++) {
        for(uint32_t j = 0; j < y; j++) {
            for(uint32_t k = 0; k < z; k++) {
                float xx = dx*float(i), yy = dy*float(j), zz=dz*float(k);

                xx-=0.5;
                yy-=0.5;
                zz-=0.5;

                float r = sqrt(xx*xx + yy*yy);
                float rho_r = cos(2*M_PI*f_m*cos(M_PI*r*0.5));

                float result = (1.-sin(M_PI*zz*0.5)) + alpha*(1+rho_r);
                result = result /(2*(1+alpha));
                //result = xx;

                int index = i + x*(j + y*k);

                #ifdef TEST_ONE_CELL
                buffer[index] = 0;
                #else
                buffer[index] = result*254.;
                #endif
            }
        }
    }

    #ifdef TEST_ONE_CELL
    int index = 6 + x*(6 + y*6);
    buffer[index] = 255.0;
    #endif
}


void sampleMarshnerLobb_BCC(unsigned char *buffer0, unsigned char *buffer1, uint32_t x, uint32_t y, uint32_t z) {
    float dx = 1./float(x);
    float dy = 1./float(y);
    float dz = 1./float(z);

    const float f_m = 6.;
    const float alpha = .25;

    int x0 = ceil(volumeSize.width),  y0 = ceil(volumeSize.height),  z0=ceil(volumeSize.depth);
    int x1 = floor(volumeSize.width), y1 = floor(volumeSize.height), z1=floor(volumeSize.depth);

    for(uint32_t i = 0; i < x; i+=2) {
        for(uint32_t j = 0; j < y; j+=2) {
            for(uint32_t k = 0; k < z; k+=2) {
                float xx = dx*float(i), yy = dy*float(j), zz=dz*float(k);

                xx-=0.5;
                yy-=0.5;
                zz-=0.5;

                float r = sqrt(xx*xx + yy*yy);
                float rho_r = cos(2*M_PI*f_m*cos(M_PI*r*0.5));

                float result = (1.-sin(M_PI*zz*0.5)) + alpha*(1+rho_r);
                result = result /(2*(1+alpha));

                int index = (i>>1) + x0*((j>>1) + y0*(k>>1));
                buffer0[index] = result*254.;
            }
        }
    }

    for(uint32_t i = 1; i < x; i+=2) {
        for(uint32_t j = 1; j < y; j+=2) {
            for(uint32_t k = 1; k < z; k+=2) {
                float xx = dx*float(i), yy = dy*float(j), zz=dz*float(k);

                xx -= 0.5;
                yy -= 0.5;
                zz -= 0.5;

                float r = sqrt(xx*xx + yy*yy);
                float rho_r = cos(2*M_PI*f_m*cos(M_PI*r*0.5));

                float result = (1.-sin(M_PI*zz*0.5)) + alpha*(1+rho_r);
                result = result /(2*(1+alpha));

                int index = (i>>1) + x0*((j>>1) + y0*(k>>1));
                buffer1[index] = result*254.;
            }
        }
    }
    int index = ( 1) + x0*((20) + y0*(1));
    // buffer0[index] = 254.;
//  index = ( 2) + x0*((1) + y0*(1));
//  buffer0[0] = 255.;
//  index = ( 1) + x0*((2) + y0*(1));
//  buffer0[index] = 255.;
}

void sampleMarshnerLobb_FCC(unsigned char *buffer0, unsigned char *buffer1, unsigned char *buffer2, unsigned char *buffer3, uint32_t x, uint32_t y, uint32_t z) {
    float dx = 1./float(x);
    float dy = 1./float(y);
    float dz = 1./float(z);

    const float f_m = 6.;
    const float alpha = .25;

    int x0 = ceil(volumeSize.width),  y0 = ceil(volumeSize.height),  z0=ceil(volumeSize.depth);
    int x1 = floor(volumeSize.width), y1 = floor(volumeSize.height), z1=ceil(volumeSize.depth);
    int x2 = floor(volumeSize.width), y2 = ceil(volumeSize.height), z2=floor(volumeSize.depth);
    int x3 = ceil(volumeSize.width), y3 = floor(volumeSize.height), z3=floor(volumeSize.depth);


    for(uint32_t i = 0; i < x; i+=2) {
        for(uint32_t j = 0; j < y; j+=2) {
            for(uint32_t k = 0; k < z; k+=2) {
                float xx = dx*float(i), yy = dy*float(j), zz=dz*float(k);

                xx-=0.5;
                yy-=0.5;
                zz-=0.5;

                float r = sqrt(xx*xx + yy*yy);
                float rho_r = cos(2*M_PI*f_m*cos(M_PI*r*0.5));

                float result = (1.-sin(M_PI*zz*0.5)) + alpha*(1+rho_r);
                result = result /(2*(1+alpha));

                int index = (i>>1) + x0*((j>>1) + y0*(k>>1));
                buffer0[index] = result*254.;
            }
        }
    }

    for(uint32_t i = 1; i < x; i+=2) {
        for(uint32_t j = 1; j < y; j+=2) {
            for(uint32_t k = 0; k < z; k+=2) {
                float xx = dx*float(i), yy = dy*float(j), zz=dz*float(k);

                xx-=0.5;
                yy-=0.5;
                zz-=0.5;

                float r = sqrt(xx*xx + yy*yy);
                float rho_r = cos(2*M_PI*f_m*cos(M_PI*r*0.5));

                float result = (1.-sin(M_PI*zz*0.5)) + alpha*(1+rho_r);
                result = result /(2*(1+alpha));

                int index = (i>>1) + x0*((j>>1) + y0*(k>>1));
                buffer1[index] = result*254.;
            }
        }
    }

    for(uint32_t i = 1; i < x; i+=2) {
        for(uint32_t j = 0; j < y; j+=2) {
            for(uint32_t k = 1; k < z; k+=2) {
                float xx = dx*float(i), yy = dy*float(j), zz=dz*float(k);

                xx-=0.5;
                yy-=0.5;
                zz-=0.5;

                float r = sqrt(xx*xx + yy*yy);
                float rho_r = cos(2*M_PI*f_m*cos(M_PI*r*0.5));

                float result = (1.-sin(M_PI*zz*0.5)) + alpha*(1+rho_r);
                result = result /(2*(1+alpha));

                int index = (i>>1) + x0*((j>>1) + y0*(k>>1));
                buffer2[index] = result*254.;
            }
        }
    }


    for(uint32_t i = 0; i < x; i+=2) {
        for(uint32_t j = 1; j < y; j+=2) {
            for(uint32_t k = 1; k < z; k+=2) {
                float xx = dx*float(i), yy = dy*float(j), zz=dz*float(k);

                xx-=0.5;
                yy-=0.5;
                zz-=0.5;

                float r = sqrt(xx*xx + yy*yy);
                float rho_r = cos(2*M_PI*f_m*cos(M_PI*r*0.5));

                float result = (1.-sin(M_PI*zz*0.5)) + alpha*(1+rho_r);
                result = result /(2*(1+alpha));

                int index = (i>>1) + x0*((j>>1) + y0*(k>>1));
                buffer3[index] = result*254.;
            }
        }
    }
}

void computeFPS()
{
    frameCount++;
    fpsCount++;

    if (fpsCount == fpsLimit)
    {
        char fps[256];
        float ifps = 1.f / (sdkGetAverageTimerValue(&timer) / 1000.f);
        sprintf(fps, "Volume Render: %3.1f fps", ifps);

        glutSetWindowTitle(fps);
        fpsCount = 0;

        fpsLimit = (int)MAX(1.f, ifps);
        sdkResetTimer(&timer);
    }
}

static  int frame_index = 0;

// render image using CUDA
void render()
{
    copyInvViewMatrix(invViewMatrix, sizeof(float4)*3);

    // map PBO to get CUDA device pointer
    uint *d_output;
    // map PBO to get CUDA device pointer
    checkCudaErrors(cudaGraphicsMapResources(1, &cuda_pbo_resource, 0));
    size_t num_bytes;
    checkCudaErrors(cudaGraphicsResourceGetMappedPointer((void **)&d_output, &num_bytes,
                                                         cuda_pbo_resource));
    //printf("CUDA mapped PBO: May access %ld bytes\n", num_bytes);

    // clear image
    checkCudaErrors(cudaMemset(d_output, 0, width*height*4));

    // call CUDA kernel, writing results to PBO
    render_kernel(gridSize, blockSize, d_output, width, height, density, brightness, transferOffset, transferScale);
    if(frame_index++ >= 10 && IN_TEST) {
        density = 0.005;
        double dAvgTime = sdkGetTimerValue(&timerInternal);
        sdkResetTimer(&timerInternal);
        printf("%f\n", dAvgTime);
        viewRotation.x  += 0.5;
        viewRotation.y  += 0.7;
        if(frame_index == 115) {
            exit_program = true;
        }
    }


    getLastCudaError("kernel failed");

    checkCudaErrors(cudaGraphicsUnmapResources(1, &cuda_pbo_resource, 0));
}

// display results using OpenGL (called by GLUT)
void display()
{
    sdkStartTimer(&timer);
    sdkStartTimer(&timerInternal);

    // use OpenGL to build view matrix
    GLfloat modelView[16];
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRotatef(-viewRotation.x, 1.0, 0.0, 0.0);
    glRotatef(-viewRotation.y, 0.0, 1.0, 0.0);
    glTranslatef(-viewTranslation.x, -viewTranslation.y, -viewTranslation.z);
    glGetFloatv(GL_MODELVIEW_MATRIX, modelView);
    glPopMatrix();

    invViewMatrix[0] = modelView[0];
    invViewMatrix[1] = modelView[4];
    invViewMatrix[2] = modelView[8];
    invViewMatrix[3] = modelView[12];
    invViewMatrix[4] = modelView[1];
    invViewMatrix[5] = modelView[5];
    invViewMatrix[6] = modelView[9];
    invViewMatrix[7] = modelView[13];
    invViewMatrix[8] = modelView[2];
    invViewMatrix[9] = modelView[6];
    invViewMatrix[10] = modelView[10];
    invViewMatrix[11] = modelView[14];

    render();

    // display results
    glClear(GL_COLOR_BUFFER_BIT);

    // draw image from PBO
    glDisable(GL_DEPTH_TEST);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
#if 0
    // draw using glDrawPixels (slower)
    glRasterPos2i(0, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, pbo);
    glDrawPixels(width, height, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
#else
    // draw using texture

    // copy from pbo to texture
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, pbo);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

    // draw textured quad
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);
    glVertex2f(0, 0);
    glTexCoord2f(1, 0);
    glVertex2f(1, 0);
    glTexCoord2f(1, 1);
    glVertex2f(1, 1);
    glTexCoord2f(0, 1);
    glVertex2f(0, 1);
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
#endif

    glutSwapBuffers();
    glutReportErrors();

    sdkStopTimer(&timer);
    sdkStopTimer(&timerInternal);

    computeFPS();
}

void idle()
{
    if(exit_program) {
        glutDestroyWindow(glutGetWindow());    
    }else{
        glutPostRedisplay();
    }
}

void keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27:
            #if defined (__APPLE__) || defined(MACOSX)
                exit(EXIT_SUCCESS);
            #else
                glutDestroyWindow(glutGetWindow());
                return;
            #endif
            break;

        case 'f':
            linearFiltering = !linearFiltering;
            setTextureFilterMode(linearFiltering);
            break;

        case '+':
            density += 0.01f;
            break;

        case '-':
            density -= 0.01f;
            break;

        case ']':
            brightness += 0.1f;
            break;

        case '[':
            brightness -= 0.1f;
            break;

        case ';':
            transferOffset += 0.01f;
            break;

        case '\'':
            transferOffset -= 0.01f;
            break;

        case '.':
            transferScale += 0.01f;
            break;

        case ',':
            transferScale -= 0.01f;
            break;

        default:
            break;
    }

    // printf("density = %.2f, brightness = %.2f, transferOffset = %.2f, transferScale = %.2f\n", density, brightness, transferOffset, transferScale);
    glutPostRedisplay();
}

int ox, oy;
int buttonState = 0;

void mouse(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN)
    {
        buttonState  |= 1<<button;
    }
    else if (state == GLUT_UP)
    {
        buttonState = 0;
    }

    ox = x;
    oy = y;
    glutPostRedisplay();
}

void motion(int x, int y)
{
    float dx, dy;
    dx = (float)(x - ox);
    dy = (float)(y - oy);

    if (buttonState == 4)
    {
        // right = zoom
        viewTranslation.z += dy / 100.0f;
    }
    else if (buttonState == 2)
    {
        // middle = translate
        viewTranslation.x += dx / 100.0f;
        viewTranslation.y -= dy / 100.0f;
    }
    else if (buttonState == 1)
    {
        // left = rotate
        viewRotation.x += dy / 5.0f;
        viewRotation.y += dx / 5.0f;
    }

    ox = x;
    oy = y;
    glutPostRedisplay();
}

int iDivUp(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

void reshape(int w, int h)
{
    width = w;
    height = h;
    initPixelBuffer();

    // calculate new grid size
    gridSize = dim3(iDivUp(width, blockSize.x), iDivUp(height, blockSize.y));

    glViewport(0, 0, w, h);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
}

void cleanup()
{
    sdkDeleteTimer(&timer);

    freeCudaBuffers();

    if (pbo)
    {
        cudaGraphicsUnregisterResource(cuda_pbo_resource);
        glDeleteBuffers(1, &pbo);
        glDeleteTextures(1, &tex);
    }
    // Calling cudaProfilerStop causes all profile data to be
    // flushed before the application exits
    checkCudaErrors(cudaProfilerStop());
}

void initGL(int *argc, char **argv)
{
    putenv( (char *) "__GL_SYNC_TO_VBLANK=0" );
    // initialize GLUT callback functions
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(width, height);
    glutCreateWindow("CUDA volume rendering");

    if (!isGLVersionSupported(2,0) ||
        !areGLExtensionsSupported("GL_ARB_pixel_buffer_object"))
    {
        printf("Required OpenGL extensions are missing.");
        exit(EXIT_SUCCESS);
    }
}

void initPixelBuffer()
{
    if (pbo)
    {
        // unregister this buffer object from CUDA C
        checkCudaErrors(cudaGraphicsUnregisterResource(cuda_pbo_resource));

        // delete old buffer
        glDeleteBuffers(1, &pbo);
        glDeleteTextures(1, &tex);
    }

    // create pixel buffer object for display
    glGenBuffers(1, &pbo);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, pbo);
    glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, width*height*sizeof(GLubyte)*4, 0, GL_STREAM_DRAW_ARB);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

    // register this buffer object with CUDA
    checkCudaErrors(cudaGraphicsGLRegisterBuffer(&cuda_pbo_resource, pbo, cudaGraphicsMapFlagsWriteDiscard));

    // create texture for display
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Load raw data from disk
void *loadRawFile(char *filename, size_t size)
{
    FILE *fp = fopen(filename, "rb");

    if (!fp)
    {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        return 0;
    }

    void *data = malloc(size);
    size_t read = fread(data, 1, size, fp);
    fclose(fp);

#if defined(_MSC_VER_)
    printf("Read '%s', %Iu bytes\n", filename, read);
#else
    printf("Read '%s', %zu bytes\n", filename, read);
#endif

    return data;
}

void runSingleTest(const char *ref_file, const char *exec_path)
{
    bool bTestResult = true;

    uint *d_output;
    checkCudaErrors(cudaMalloc((void **)&d_output, width*height*sizeof(uint)));
    checkCudaErrors(cudaMemset(d_output, 0, width*height*sizeof(uint)));

    float modelView[16] =
    {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 4.0f, 1.0f
    };

    invViewMatrix[0] = modelView[0];
    invViewMatrix[1] = modelView[4];
    invViewMatrix[2] = modelView[8];
    invViewMatrix[3] = modelView[12];
    invViewMatrix[4] = modelView[1];
    invViewMatrix[5] = modelView[5];
    invViewMatrix[6] = modelView[9];
    invViewMatrix[7] = modelView[13];
    invViewMatrix[8] = modelView[2];
    invViewMatrix[9] = modelView[6];
    invViewMatrix[10] = modelView[10];
    invViewMatrix[11] = modelView[14];

    // call CUDA kernel, writing results to PBO
    copyInvViewMatrix(invViewMatrix, sizeof(float4)*3);

    // Start timer 0 and process n loops on the GPU
    int nIter = 10;

    for (int i = -1; i < nIter; i++)
    {
        if (i == 0)
        {
            cudaDeviceSynchronize();
            sdkStartTimer(&timer);
        }

        render_kernel(gridSize, blockSize, d_output, width, height, density, brightness, transferOffset, transferScale);
    }

    cudaDeviceSynchronize();
    sdkStopTimer(&timer);
    // Get elapsed time and throughput, then log to sample and master logs
    double dAvgTime = sdkGetTimerValue(&timer)/(nIter * 1000.0);
    printf("volumeRender, Throughput = %.4f MTexels/s, Time = %.5f s, Size = %u Texels, NumDevsUsed = %u, Workgroup = %u\n",
           (1.0e-6 * width * height)/dAvgTime, dAvgTime, (width * height), 1, blockSize.x * blockSize.y);


    getLastCudaError("Error: render_kernel() execution FAILED");
    checkCudaErrors(cudaDeviceSynchronize());

    unsigned char *h_output = (unsigned char *)malloc(width*height*4);
    checkCudaErrors(cudaMemcpy(h_output, d_output, width*height*4, cudaMemcpyDeviceToHost));

    sdkSavePPM4ub("volume.ppm", h_output, width, height);
    bTestResult = sdkComparePPM("volume.ppm", sdkFindFilePath(ref_file, exec_path), MAX_EPSILON_ERROR, THRESHOLD, true);

    cudaFree(d_output);
    free(h_output);
    cleanup();

    exit(bTestResult ? EXIT_SUCCESS : EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
    pArgc = &argc;
    pArgv = argv;

    char *ref_file = NULL;

#if defined(__linux__)
    setenv ("DISPLAY", ":0", 0);
#endif




    if (checkCmdLineFlag(argc, (const char **)argv, "file"))
    {
        getCmdLineArgumentString(argc, (const char **)argv, "file", &ref_file);
        fpsLimit = frameCheckNumber;
    }

    if (ref_file)
    {
        findCudaDevice(argc, (const char **)argv);
    }
    else
    {
        // First initialize OpenGL context, so we can properly set the GL for CUDA.
        // This is necessary in order to achieve optimal performance with OpenGL/CUDA interop.
        initGL(&argc, argv);

        findCudaDevice(argc, (const char **)argv);
    }

    // parse arguments
    char *filename;

    if (getCmdLineArgumentString(argc, (const char **) argv, "volume", &filename))
    {
        volumeFilename = filename;
    }

    int n;

    if (checkCmdLineFlag(argc, (const char **) argv, "size"))
    {
        n = getCmdLineArgumentInt(argc, (const char **) argv, "size");
        volumeSize.width = volumeSize.height = volumeSize.depth = n;
    }

    if (checkCmdLineFlag(argc, (const char **) argv, "xsize"))
    {
        n = getCmdLineArgumentInt(argc, (const char **) argv, "xsize");
        volumeSize.width = n;
    }

    if (checkCmdLineFlag(argc, (const char **) argv, "ysize"))
    {
        n = getCmdLineArgumentInt(argc, (const char **) argv, "ysize");
        volumeSize.height = n;
    }

    if (checkCmdLineFlag(argc, (const char **) argv, "zsize"))
    {
        n= getCmdLineArgumentInt(argc, (const char **) argv, "zsize");
        volumeSize.depth = n;
    }

    // load volume data
    char *path = sdkFindFilePath(volumeFilename, argv[0]);

    if (path == 0)
    {
        printf("Error finding file '%s'\n", volumeFilename);
        exit(EXIT_FAILURE);
    }


#ifdef CC_LATTICE
    size_t size = volumeSize.width*volumeSize.height*volumeSize.depth*sizeof(VolumeType);
    void *h_volume = malloc(size);

    sampleMarshnerLobb_CC((unsigned char*)h_volume, volumeSize.width,volumeSize.height,volumeSize.depth);
    initCuda(h_volume, h_volume, h_volume,h_volume, volumeSize);
    free(h_volume);
#endif

#ifdef BCC_LATTICE
    // Break this into two lattices
    int x0 = ceil(volumeSize.width),  y0 = ceil(volumeSize.height),  z0=ceil(volumeSize.depth);
    int x1 = floor(volumeSize.width), y1 = floor(volumeSize.height), z1=floor(volumeSize.depth);


    void *h_volume0 = malloc(x0*y0*z0*sizeof(VolumeType));
    void *h_volume1 = malloc(x1*y1*z1*sizeof(VolumeType));


    sampleMarshnerLobb_BCC((unsigned char*)h_volume0, (unsigned char*)h_volume1,volumeSize.width*2,volumeSize.height*2,volumeSize.depth*2);

    initCuda(h_volume0, h_volume1, h_volume1, h_volume1, volumeSize);
    free(h_volume0);
    free(h_volume1);
#endif

#ifdef FCC_LATTICE
    // Break this into two lattices
    int x0 = ceil(volumeSize.width),  y0 = ceil(volumeSize.height),  z0=ceil(volumeSize.depth);
    int x1 = floor(volumeSize.width), y1 = floor(volumeSize.height), z1=ceil(volumeSize.depth);
    int x2 = floor(volumeSize.width), y2 = ceil(volumeSize.height), z2=floor(volumeSize.depth);
    int x3 = ceil(volumeSize.width), y3 = floor(volumeSize.height), z3=floor(volumeSize.depth);


    void *h_volume0 = malloc(x0*y0*z0*sizeof(VolumeType));
    void *h_volume1 = malloc(x1*y1*z1*sizeof(VolumeType));
    void *h_volume2 = malloc(x2*y2*z2*sizeof(VolumeType));
    void *h_volume3 = malloc(x3*y3*z3*sizeof(VolumeType));

    sampleMarshnerLobb_FCC((unsigned char*)h_volume0, (unsigned char*)h_volume1, (unsigned char*)h_volume2, (unsigned char*)h_volume3, volumeSize.width*2,volumeSize.height*2,volumeSize.depth*2);

    initCuda(h_volume0, h_volume1, h_volume2, h_volume3, volumeSize);
    free(h_volume0);
    free(h_volume1);
    free(h_volume2);
    free(h_volume3);
#endif
    dim3 _blockSize(32*8);
    dim3 _gridSize(1024*32/8);

    random_test(_gridSize, _blockSize);
    // initCuda(h_volume, h_volume, h_volume, h_volume, volumeSize);
    // free(h_volume);


    exit(-1);

    sdkCreateTimer(&timer);
    sdkCreateTimer(&timerInternal);

    // printf("Press '+' and '-' to change density (0.01 increments)\n"
    //        "      ']' and '[' to change brightness\n"
    //        "      ';' and ''' to modify transfer function offset\n"
    //        "      '.' and ',' to modify transfer function scale\n\n");

    // calculate new grid size
    gridSize = dim3(iDivUp(width, blockSize.x), iDivUp(height, blockSize.y));

    if (ref_file)
    {
        runSingleTest(ref_file, argv[0]);
    }
    else
    {
        // This is the normal rendering path for VolumeRender
        glutDisplayFunc(display);
        glutKeyboardFunc(keyboard);
        glutMouseFunc(mouse);
        glutMotionFunc(motion);
        glutReshapeFunc(reshape);
        glutIdleFunc(idle);

        initPixelBuffer();

#if defined (__APPLE__) || defined(MACOSX)
        atexit(cleanup);
#else
        glutCloseFunc(cleanup);
#endif

        glutMainLoop();
    }
}

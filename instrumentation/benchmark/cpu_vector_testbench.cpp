
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


#define BATCH_SIZE 1024
#define SAMPLES 1000*BATCH_SIZE // ~1 million
#define BURNIN 10*BATCH_SIZE

#define DIMENSION 3

#if defined(LATTICE_cc) || defined(LATTICE_CC)
    #define CARTESIAN_CUBIC
#elif defined(LATTICE_bcc) || defined(LATTICE_BCC)
    #define BODY_CENTERED_CUBIC
#elif defined(LATTICE_fcc) || defined(LATTICE_FCC)
    #define FACE_CENTERED_CUBIC
#endif

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
extern float4 __reconstruct__(float4 x, float4 y, _lattice_info *);
#elif DIMENSION == 3 && !defined(TEST_TRILINEAR)
extern float4 __reconstruct__(float4 x, float4 y, float4 z, _lattice_info *);
#elif DIMENSION == 4
extern float4 __reconstruct__(float4 x, float4 y, float4 z, float4 w, _lattice_info *);
#endif
}
struct _coset_info *allocate_3d_coset(unsigned int x, unsigned int y, unsigned int z) {
	struct _coset_info *coset = (_coset_info *)calloc(sizeof(_coset_info), 1);
	coset->buffer = (float *)calloc(x*z, y*sizeof(float));
	coset->bounds = (uint32_t *)calloc(3, sizeof(uint32_t));
	coset->bounds[0] = x;
	coset->bounds[1] = y;
	coset->bounds[2] = z;
	int idx = 0;
	for(int i = 0; i < x; i++)
	    for(int j = 0; j < y; j++)
	        for(int k = 0; k < z; k++) {
	            coset->buffer[idx++] = 0.0f;
	        }
	return coset;
}

struct _lattice_info *allocate_cc(unsigned int x, unsigned int y, unsigned int z) {
	struct _coset_info *coset_0 = allocate_3d_coset(x+1, y+1, z+1);
	struct _lattice_info *lattice_instance = (_lattice_info *)calloc(1, sizeof(_lattice_info));


    const int wh = x + 1;
    const int ww = y + 1;
    const int wd = z + 1;

    // Set 10,10,10 to 1.0
    coset_0->buffer[10 + ww * (10 + wd* 10)] = 1.0;

	lattice_instance->coset_count = 1;
	lattice_instance->dimension = 3;
	lattice_instance->cosets = (_coset_info *) calloc(1, sizeof(_coset_info));
	memcpy(&lattice_instance->cosets[0], coset_0, sizeof(_coset_info));


	return lattice_instance;
}

void free_lattice(struct _lattice_info *l) {
    const int ncosets = l->coset_count;
    for(int i = 0; i < ncosets; i++) {
        free(l->cosets[i].bounds);
        free(l->cosets[i].buffer);
    }
    free(l->cosets);
    free(l);
}

struct _lattice_info *allocate_bcc(unsigned int x, unsigned int y, unsigned int z) {
	const unsigned int extend_x = x % 2;
	const unsigned int extend_y = y % 2;
	const unsigned int extend_z = z % 2;

	struct _coset_info *coset_0 = allocate_3d_coset((x/2) + 1, (y/2) + 1, (z/2) + 1);
	struct _coset_info *coset_1 = allocate_3d_coset((x/2) + extend_x, (y/2) + extend_y, (z/2) + extend_z);
	struct _lattice_info *lattice_instance = (_lattice_info *)calloc(1, sizeof(_lattice_info));


    const int wh = x/2 + 1;
    const int ww = y/2 + 1;
    const int wd = z/2 + 1;

    // Set 10,10,10 to 1.0
    coset_0->buffer[5 + ww * (5 + wd* 5)] = 1.0;

	lattice_instance->coset_count = 2;
	lattice_instance->dimension = 3;
	lattice_instance->cosets = (_coset_info *) calloc(lattice_instance->coset_count, sizeof(_coset_info));
	memcpy(&lattice_instance->cosets[0], coset_0, sizeof(_coset_info));
	memcpy(&lattice_instance->cosets[1], coset_1, sizeof(_coset_info));

	return lattice_instance;
}

struct _lattice_info *allocate_fcc(unsigned int x, unsigned int y, unsigned int z) {
	const unsigned int extend_x = x % 2;
	const unsigned int extend_y = y % 2;
	const unsigned int extend_z = z % 2;

	struct _coset_info *coset_0 = allocate_3d_coset((x/2) + 1, (y/2) + 1, (z/2) + 1);
	struct _coset_info *coset_1 = allocate_3d_coset((x/2) + extend_x, (y/2) + extend_y, (z/2));
	struct _coset_info *coset_2 = allocate_3d_coset((x/2) + extend_x, (y/2), (z/2) + extend_z);
	struct _coset_info *coset_3 = allocate_3d_coset((x/2) + extend_x, (y/2) + extend_y, (z/2));

	struct _lattice_info *lattice_instance = (_lattice_info *)calloc(1, sizeof(_lattice_info));
	lattice_instance->coset_count = 4;
	lattice_instance->dimension = 3;
	lattice_instance->cosets = (_coset_info *) calloc(lattice_instance->coset_count, sizeof(_coset_info));

    const int wh = x/2 + 1;
    const int ww = y/2 + 1;
    const int wd = z/2 + 1;

    // Set 10,10,10 to 1.0
    coset_0->buffer[5 + ww * (5 + wd* 5)] = 1.0;

	memcpy(&lattice_instance->cosets[0], coset_0, sizeof(_coset_info));
	memcpy(&lattice_instance->cosets[1], coset_1, sizeof(_coset_info));
	memcpy(&lattice_instance->cosets[2], coset_2, sizeof(_coset_info));
	memcpy(&lattice_instance->cosets[3], coset_3, sizeof(_coset_info));

	return lattice_instance;
}


struct _lattice_info *allocate_lattice(float scale, int *lattice_bound) {
	struct _lattice_info *lattice_instance = NULL;
	*lattice_bound = (int)(1/scale);
//
//#if defined(LATTICE_CP)
//	if(scale < 0) lattice_bound = 128;
//	allocate_cp(lattice_bound, lattice_bound, &coset0);
//#elif defined(LATTICE_QC)
//	lookup_func = _access_qc;
//	if(scale < 0) lattice_bound = 180;
//	allocate_qc(lattice_bound, lattice_bound, &coset0, &coset1);
#if defined(CARTESIAN_CUBIC)
	if(scale < 0) *lattice_bound = 128;
	lattice_instance = allocate_cc(*lattice_bound, *lattice_bound, *lattice_bound);
#elif defined(BODY_CENTERED_CUBIC)
	if(scale < 0) *lattice_bound = 202;
	lattice_instance = allocate_bcc(*lattice_bound, *lattice_bound, *lattice_bound);
#elif defined(FACE_CENTERED_CUBIC)
	if(scale <0) *lattice_bound = 161;
	lattice_instance = allocate_fcc(*lattice_bound, *lattice_bound, *lattice_bound);
//#elif defined(LATTICE_C4)
//	lookup_func = _access_c4;
//	if(scale <0) lattice_bound = 128;
//	allocate_c4(lattice_bound, lattice_bound, lattice_bound, lattice_bound, &coset0);
//#elif defined(LATTICE_D4)
//	lookup_func = _access_d4;
//	if(scale <0) lattice_bound = 214;
//	allocate_d4(lattice_bound, lattice_bound, lattice_bound, lattice_bound, &coset0, &coset1);
#endif //202 BCC, 128 CC, 214 D4S, 128 C4, 161 FCC
	return lattice_instance;
}

float lattice_accessor_3d(int x, int y, int z) {
    if(x == y && y == z && x == 10)
        return 1.0;
    return 0.0;
}


int main(int argc, char const *argv[]) {
	int lattice_bound = 0;
    struct _lattice_info *lattice_instance = allocate_lattice(-1, &lattice_bound);
	cpu_set_t set;

    if(argc < 2) {
        printf("Input to this test bench should include the reconstruction filter width. I.e\n"
               "./valtest [width]\n");
        return -1;
    }

    int idx = 0;
    float fwidth = atof(argv[1]);
    printf("[");
    for(float z = 10.0 - fwidth; z < 10.0 + fwidth; z += (fwidth/50.0) ) {
        printf("[");
        for(float y = 10.0 - fwidth; y < 10.0 + fwidth; y += (fwidth/50.0) ) {
            printf("[");
            for(float x = 10.0 - fwidth; x < 10.0 + fwidth; x += (fwidth/50.0) ) {
                float4 xv = {5.3, 5.1, 5, 5.3};
                float4 yv = {5.4, 5.2, 6, 5.2};
                float4 zv = {5.5, 5.3, 7, 5.1};
//                float4 eval = {5.0, 5.3, 7, 5.1};

                xv[idx] = x;
                yv[idx] = y;
                zv[idx] = z;
//                printf("%f, %f, %f->", x, y, z);
                float4 eval =  __reconstruct__(xv,yv,zv, lattice_instance);

                printf("%f,", eval[idx]);
                idx = (idx+1) % 4;
            }
            printf("],\n");
        }
        printf("],\n");

    }
    printf("]\n");
    return 0;
}
//
//int main(int argc, char const *argv[]) {
//    //
//    float4 x = {10.0, 9.5, 9.25, 9.0};
//    float4 y = {10.0, 9.25, 9.25, 9.0};
//    float4 z = {10.0, 9.5, 9.25, 9.0};
//    float4 res =  __reconstruct__(x, y, z, lattice_accessor_3d);
//
//    printf("%f\n", res[0]);
//    printf("%f\n", res[1]);
//    printf("%f\n", res[2]);
//    printf("%f\n", res[3]);
//}

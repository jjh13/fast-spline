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



float __attribute__ ((noinline)) __reconstruct___(float x, float y, float z, _lattice_info *l) {
	float *base_addr = l->cosets[0].buffer;

	// Note this should be floor -- this works fine for positive values of x,y,z
	// but will fail for negative
	const int i = int(x);
	const int j = int(y);
	const int k = int(z);
	const float tx = x - i;
	const float ty = y - j;
	const float tz = z - k;

	const uint32_t bx = l->cosets[0].bounds[0];
	const uint32_t by = l->cosets[0].bounds[1];
	const uint32_t bz = l->cosets[0].bounds[2];

    /*k
    k*bz + j
    (k*bz + j)*by + i
*/
	const uint32_t index = i + by * (j + bz* k); 
	const float d000 = base_addr[index];
	const float d001 = base_addr[index + 1];
	const float d010 = base_addr[index + by];
	const float d011 = base_addr[index + by + 1];
	base_addr += by*bz; 

	const float d100 = base_addr[index];
	const float d101 = base_addr[index + 1];
	const float d110 = base_addr[index + by];
	const float d111 = base_addr[index + by  + 1];
	
	const float xm1 = 1. - tx, ym1 = 1. - ty, zm1 = 1. - tz;
	return  xm1 * ym1 * zm1 * d000 +
                tx  * ym1 * zm1 * d100 + 
		xm1 * ty  * zm1 * d010 + 
		tx  * ty  * zm1 * d110 +
		xm1 * ym1 * tz  * d001 + 
		tx  * ym1 * tz  * d101 +
	        xm1 * ty  * tz  * d011 + 
		tx  * ty  * tz  * d111;
	

}


extern "C" {
#if DIMENSION == 2
extern float __reconstruct__(float x, float y, _lattice_info *);
#elif DIMENSION == 3 && !defined(TEST_TRILINEAR)
extern float __reconstruct__(float x, float y, float z, _lattice_info *);
#elif DIMENSION == 4
extern float __reconstruct__(float x, float y, float z, float w, _lattice_info *);
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
	            coset->buffer[idx++] = 1.0f;
	        }
	return coset;
}

struct _lattice_info *allocate_cc(unsigned int x, unsigned int y, unsigned int z) {
	struct _coset_info *coset_0 = allocate_3d_coset(x+1, y+1, z+1);
	struct _lattice_info *lattice_instance = (_lattice_info *)calloc(1, sizeof(_lattice_info));

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
	if(scale < 0) *lattice_bound = 32;
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
	

	// Lock this thread to a single core
	CPU_ZERO(&set);
	CPU_SET(0, &set);
	sched_setaffinity(getpid(), sizeof(set), &set);
	// Do some burn in
	float pos[DIMENSION];
	float noopt = 0, last = 0;
//
//	for(int i = 0; i < BURNIN; i++) {
//        float ret = 0.0;
//		#pragma unroll
//		for(int k = 0; k < DIMENSION; k++) { pos[k] = ((float)(lattice_bound))*(float)rand()/(float)(RAND_MAX); }
//
//		#if DIMENSION == 2
//		ret = __reconstruct__(pos[0], pos[1], lattice_instance);
//		#elif DIMENSION == 3
//		ret = __reconstruct__(pos[0], pos[1], pos[2], lattice_instance);
//		#elif DIMENSION == 4
//		ret = __reconstruct__(pos[0], pos[1], pos[2], pos[3], lattice_instance);
//		#endif
//
//		noopt += ret;
//	}

    const int padding = 4;
    std::vector<std::tuple<float, float, float>> eval_points;
    for(int i = 0; i < SAMPLES + BURNIN; i++) {
        eval_points.push_back({
           padding + (((float) (lattice_bound - 2*padding)) * (float) rand() / float(RAND_MAX)),
           padding + (((float) (lattice_bound - 2*padding)) * (float) rand() / float(RAND_MAX)),
           padding + (((float) (lattice_bound - 2*padding)) * (float) rand() / float(RAND_MAX))
        });
    }

    for(int i = 0; i < SAMPLES/BATCH_SIZE + BURNIN/BATCH_SIZE ; i++) {
        const int offset = i * BATCH_SIZE;
        auto const t0 = std::chrono::steady_clock::now();

        for(int j = 0; j < BATCH_SIZE; j++) {
            pos[0] = std::get<0>(eval_points[offset+j]);
            #if DIMENSION == 2
            pos[1] = std::get<1>(eval_points[offset+j]);
            last = __reconstruct__(pos[0], pos[1], lattice_instance);
            #elif DIMENSION == 3
            pos[1] = std::get<1>(eval_points[offset+j]);
            pos[2] = std::get<2>(eval_points[offset+j]);
            last =  __reconstruct__(pos[0], pos[1], pos[2], lattice_instance);
//
//            const uint32_t bx = lattice_instance->cosets[0].bounds[0];
//            const uint32_t by = lattice_instance->cosets[0].bounds[1];
//            const uint32_t bz = lattice_instance->cosets[0].bounds[2];
//
//            int ii = (int)pos[0];
//            int jj = (int)pos[1];
//            int kk = (int)pos[2];
//            const uint32_t index = ii + by * (jj + bz* kk);
//            printf("%f, %f, %f:  %f (%d)\n", pos[0], pos[1], pos[2], last, index);

            #elif DIMENSION == 4
            pos[1] = std::get<1>(eval_points[offset+j]);
            pos[2] = std::get<2>(eval_points[offset+j]);
            pos[3] = std::get<3>(eval_points[offset+j]);
            last = __reconstruct__(pos[0], pos[1], pos[2], pos[3], lattice_instance);
            #endif
            noopt += last;
        }
        double const time0 = std::chrono::duration_cast<std::chrono::duration<double>>(
                std::chrono::steady_clock::now() - t0)
                .count();

        double mean_time = (double)(BATCH_SIZE) / time0;
        if( i > BURNIN/BATCH_SIZE)
	        printf("%f,", mean_time*1e-6);
    }

    free_lattice(lattice_instance);
    return (int) noopt;
}
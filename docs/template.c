

typedef internal_type float;
typedef in_type float;
typedef out_type float;


out_type __reconstruct__(in_type x, in_type y, in_type z, lattice l*) {

    ls_lut[bsp_index?, poly_index?, idx]

    // GEN
    index_type i,j,k = rho(x,y,z);

    // If subregions > 1
    float x00, x01, x02;
    float x10, x11, x12;
    float x20, x21, x22;
    float s0, s1, s2;

    // * number of cosets
	const index_type bx = l->cosets[0].bounds[0];
	const index_type by = l->cosets[0].bounds[1];
	const index_type bz = l->cosets[0].bounds[2];

    // * pipeline depth
    float d[pipeline_depth];

    // * coset
    float *buffer0;

    // if COSET_RAW
    index_type lindex = 0;


    x = x - i;
    y = y - j;
    z = z - k;

    // If subregions > 1
    {
        index_type bsp_index =  0;
        plane_index |= (x * 1.0 + y * 0.5  - 0.5) < 0 ? 1 : 0;
        plane_index |= (x * 1.0 + y * 0.5  - 0.5) < 0 ? 2 : 0;
        plane_index |= (x * 1.0 + y * 0.5  - 0.5) < 0 ? 4 : 0;
        plane_index |= (x * 1.0 + y * 0.5  - 0.5) < 0 ? 8 : 0;

        bsp_index = bsp_index % 100;

        x00 ...  x22 = xform_index[bsp_index];

        tx = x00*x + x01*y + x02*z + s0;
        ty = x10*x + x11*y + x12*z + s1;
        z = x20*x + x21*y + x22*z + s2;
        x = tx;
        y = ty;

        // If ref subregions > 1
        poly_index = polyindex(bsp_index);

    }



    // for idx in range pipeline_depth
    // if ABSTRACT
    ii,jj,kk = ls_lut[bsp_index?, poly_index?, idx]
    d[idx % pipeline_depth] = memfetch(ii + i, jj + j, kk + k);



}
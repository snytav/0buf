/*
 * kernels.h
 *
 *  Created on: Jun 2, 2018
 *      Author: snytav
 */

#ifndef KERNELS_H_
#define KERNELS_H_

#include "cell.h"

__device__ double cuda_atomicAdd(double *address, double val);

template <template <class Particle> class Cell >
__global__
void GPU_getCellEnergy(
		Cell<Particle>  **cells,double *d_ee,
		double *d_Ex,double *d_Ey,double *d_Ez);




#endif /* KERNELS_H_ */

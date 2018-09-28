/*
 * gpucell.cu
 *
 *  Created on: Jul 23, 2018
 *      Author: snytav
 */

#include "cell.h"


void dbgPrintGPUParticleAttribute(Cell *d_c,int n_particle,int attribute,char *name )
{
    double t;
    Cell *h_c;
    int shift = (attribute + n_particle*sizeof(Particle)/sizeof(double));
    int err;

    h_c = new Cell;

    err = (int) MemoryCopy(h_c,d_c,sizeof(Cell),DEVICE_TO_HOST);
    if(err != cudaSuccess)
        {
        	printf("pyFieldsToGPU err %d %s \n",err,getErrorString(err));
        	exit(0);
        }
    double *vec = h_c->doubParticleArray + shift;

    MemoryCopy(&t,vec,sizeof(double),DEVICE_TO_HOST);

    printf("%s %10.3e \n",name,t);
}


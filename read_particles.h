/*
 * read_particles.h
 *
 *  Created on: Apr 9, 2016
 *      Author: snytav
 */

#ifndef READ_PARTICLES_H_
#define READ_PARTICLES_H_

#include <stdio.h>


#define SORTS 3

typedef struct ParticleArrays{
	double *dbg_x,*dbg_y,*dbg_z,
		   *dbg_px,*dbg_py,*dbg_pz,q_m,*m;

    int total;
} ParticleArrays;

typedef struct ParticleFloatArrays{
	float *dbg_x,*dbg_y,*dbg_z,
		   *dbg_px,*dbg_py,*dbg_pz,q_m,m;
    int total;
} ParticleFloatArrays;

typedef ParticleArrays ParticleArraysGroup[SORTS];

typedef ParticleFloatArrays ParticleFloatArraysGroup[SORTS];

      int readFortranBinaryArray(FILE *f, double* d);

      void InitBinaryParticlesArrays(char *fn,int nt,
    		  ParticleArrays *ions,
    		  ParticleArrays *electrons,
    		  ParticleArrays *beam_electrons,
    		  int Nx,int Ny,int Nz,
    		  int beam_plasma);

//    	     double **dbg_el_x, double **dbg_el_y, double **dbg_el_z,
//    	     double **dbg_el_px,double **dbg_el_py,double **dbg_el_pz,double *el_q_m,double *el_m,
//    	     int *total_electrons,
//    	     double **dbg_beam_x, double **dbg_beam_y, double **dbg_beam_z,
//    	     double **dbg_beam_px,double **dbg_beam_py,double **dbg_beam_pz,
//    	     double *beam_q_m,double *beam_m,
//    	     int *total_beam_electrons,int Nx,int Ny,int Nz);

      void AllocateBinaryParticlesArrays(
     			  ParticleArrays *ions,ParticleArrays *electrons,ParticleArrays *beam_electrons);

      void AllocateBinaryParticlesArraysFloat(
      			  ParticleFloatArrays *ions,ParticleFloatArrays *electrons,ParticleFloatArrays *beam_electrons);

      int AllocateDeviceParticleDiagnosticPointers(ParticleFloatArraysGroup **d_pfag,
     			                                       ParticleFloatArraysGroup *host_copy_d_pfag,
     			                                       ParticleFloatArraysGroup *pfag);

#endif /* READ_PARTICLES_H_ */

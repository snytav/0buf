/*
 * diagnose.h
 *
 *  Created on: Jun 6, 2016
 *      Author: snytav
 */

#ifndef DIAGNOSE_H_
#define DIAGNOSE_H_

#include "read_particles.h"

#define DIAGNOSE_NX 50
#define DIAGNOSE_NY 50

typedef float diagnose_matrix[DIAGNOSE_NX][DIAGNOSE_NY];

typedef enum DIAG_OP{CHARGE=0,
	         VELX,
	         VELY,
	         VELZ
            } particle_diagnostic_operation;

            int TransformParticlesToLists(int sorts,
                    ParticleArraysGroup pag,
                    float *particles,
                    ParticleFloatArraysGroup pfag
            		);

            int Particle2D_Distributions(int sorts,char *shot_name,int nt,
            		 ParticleFloatArraysGroup pfag,
            		 double Lx,double Ly
            		 );

int write_particles_component(
		float *x,float *y,float *z,float *pu,float *pv, float *pw,
		double *x0,double *y0, double *z0, double *pu0, double *pv0, double *pw0,
		int num, char *comp_name, char *shot_name, int nt
		);

int write_particles(ParticleFloatArraysGroup *pfag,ParticleArraysGroup *init,char *variant_name,int nt);

//int write_field_component(int nt,double *d_f,char *where,char *name,int size);

#endif /* DIAGNOSE_H_ */

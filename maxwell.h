/*
 * init.h
 *
 *  Created on: Apr 14, 2016
 *      Author: snytav
 */

#ifndef MAXWELL_H_
#define MAXWELL_H_

#include "particle.h"

#include "read_particles.h"

#include <vector>

double rnd_uniform();

double rnd_gaussian(double a,double b);

int getMassCharge(ParticleArrays *ions,ParticleArrays *electrons,ParticleArrays *beam_electrons,
		double ni,double rbd,int lp);

int InitUniformMaxwellianParticles(int beamf,int jmb,
				   double tex0,double tey0,double tez0,
				   double beam_lx, double beam_ly, double beam_lz,int *jmb_real,
                   double lx,double ly,double lz, int meh,double Tb,double rimp,double rbd,
				   double *xi,double *yi, double *zi,double *ui,double *vi, double *wi,
				   double *xb,double *yb, double *zb,double *ub,double *vb, double *wb,
				   double *xf,double *yf, double *zf,double *uf,double *vf, double *wf
				  );

int AddBeamParticles(int jmb,
				   double tex0,double tey0,double tez0,
				   double beam_lx, double beam_ly,int *jmb_real,
                   double lx,double ly,double lz, int meh,double Tb,double rimp,double rbd,
				   double *xb,double *yb, double *zb,double *ub,double *vb, double *wb
				  );

int writeParamsFile(double tex0,double tey0,double tez0,
                         double Tb,double rimp,
                         double rbd,double ni,
	                      double lx,double ly,double lz,
	                      int lp,int nx,int ny,int nz,
	                      double tau,double B0,
	                      double bx,double by,double bz,
	                      double py,double pz,
	                      int beam_plasma,int start_from_file,
	                      int ts,int ms,int phase
	                      );

int getUniformMaxwellianParticles(
		std::vector<Particle>  & ion_vp,
		                           std::vector<Particle>  & el_vp,
		                           std::vector<Particle>  & beam_vp
		                           );

int convertParticleArraysToSTLvector(
    		  double *dbg_x,
    		  double *dbg_y,
			  double *dbg_z,
			  double *dbg_px,
			  double *dbg_py,
			  double *dbg_pz,
			  double q_m,
			  double m,
			  int total_particles,
			  particle_sorts sort,
    		  std::vector<Particle> & vp
    		  );

#endif /* INIT_H_ */

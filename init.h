/*
 * init.h
 *
 *  Created on: Apr 14, 2016
 *      Author: snytav
 */

#ifndef INIT_H_
#define INIT_H_

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


#endif /* INIT_H_ */

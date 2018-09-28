/*
 * params.h
 *
 *  Created on: Sep 14, 2016
 *      Author: snytav
 */

#ifndef PARAMS_H_
#define PARAMS_H_

int readParameterFile(int *beam_plasma,int *start_from_file,
		              double *tex0,double *tey0,double *tez0,
		              double *Tb,double *rimp,
		              double *rbd,double *ni,
			          double *lx,double *ly,double *lz,
			          int *lp,int *nx,int *ny,int *nz,
			          double *tau,double *B0,
			          int *np,double *bx,double *by,double *bz,
			          double *pl_y,double *pl_z,
			          int *ts,int *ms,int *ph);

#endif /* PARAMS_H_ */

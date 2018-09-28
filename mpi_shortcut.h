/*
 * mpi_shortcut.h
 *
 *  Created on: Nov 5, 2014
 *      Author: snytav
 */

#ifndef MPI_SHORTCUT_H_
#define MPI_SHORTCUT_H_

int InitMPI(int argc,char *argv[]);

int sumMPI(int size,double *jx,double *jy,double *jz);

int sumMPIenergy(double *e);

int  CloseMPI();

int getRank();

int getSize();

#endif /* MPI_SHORTCUT_H_ */

/*
 * service_functions.h
 *
 *  Created on: Apr 30, 2018
 *      Author: snytav
 */



#ifndef SERVICE_FUNCTIONS_H_
#define SERVICE_FUNCTIONS_H_

#include <string>

double CheckArraySilent	(double* a, double* dbg_a,int size);

void read3Darray(char* name, double* d);

void read3Darray(std::string name, double* d);

void get_load_data_file_names(
		std::string & t_jxfile,
		std::string & t_jyfile,
		std::string & t_jzfile,
		std::string & t_d_jxfile,
		std::string & t_d_jyfile,
		std::string & t_d_jzfile,
		std::string & t_np_jxfile,
		std::string & t_np_jyfile,
		std::string & t_np_jzfile,
		std::string & t_qxfile,
		std::string & t_qyfile,
		std::string & t_qzfile,int nt);

double get_meminfo(void);
double get_meminfo1(void);

int setPrintfLimit();

void cudaMalloc3D(double **X,double **Y,double**Z,int nx,int ny,int nz);

void InitGPUFields(
		double **d_Ex,double **d_Ey,double **d_Ez,
		double **d_Hx,double **d_Hy,double **d_Hz,
		double **d_Jx,double **d_Jy,double **d_Jz,
		double **d_npJx,double **d_npJy,double **d_npJz,
        double **d_Qx,double **d_Qy,double **d_Qz,
        double *Ex,double *Ey,double *Ez,
		double *Hx,double *Hy,double *Hz,
		double *Jx,double *Jy,double *Jz,
		double *npJx,double *npJy,double *npJz,
		double *Qx,double *Qy,double *Qz,
		int Nx,int Ny,int Nz
        );

void copyFieldsToGPU(
						double *d_Ex,double *d_Ey,double *d_Ez,
						double *d_Hx,double *d_Hy,double *d_Hz,
						double *d_Jx,double *d_Jy,double *d_Jz,
						double *d_npJx,double *d_npJy,double *d_npJz,
						double *d_Qx,double *d_Qy,double *d_Qz,
						double *Ex,double *Ey,double *Ez,
		        		double *Hx,double *Hy,double *Hz,
		        		double *Jx,double *Jy,double *Jz,
		        		double *npJx,double *npJy,double *npJz,
		                double *Qx,double *Qy,double *Qz,
		                int Nx,int Ny,int Nz
		);



#endif /* SERVICE_FUNCTIONS_H_ */

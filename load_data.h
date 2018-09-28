/*
 * load_data.h
 *
 *  Created on: Jun 9, 2018
 *      Author: snytav
 */

#ifndef LOAD_DATA_H_
#define LOAD_DATA_H_


#include "particle.h"
#include <vector>

#include "maxwell.h"

#include <string>



std::string getMumuFileName(int nt);

int readFortranBinaryArray(FILE *f, double* d);

FILE *readPreliminary3Darrays(std::string fn,int nt,int nx,int ny,int nz);

void debugPrintParticleCharacteristicArray(double *p_ch,int np,int nt,std::string name,int sort);

int readBinaryParticleArraysOneSort(
		  FILE *f,
		  double **dbg_x,
		  double **dbg_y,
		  double **dbg_z,
		  double **dbg_px,
		  double **dbg_py,
		  double **dbg_pz,
		  double *qq_m,
		  double *mm,
		  int nt,
		  int sort
		  );

int getParticlesOneSortFromFile(
		                          FILE *f,
                                particle_sorts sort,
                                int nt,
                                std::vector<Particle> & vp,
                                double *q_m,
                                double *m
                                );

void readBinaryParticlesOneSort(FILE *f,std::vector<Particle> & vp,
			                                  particle_sorts sort,int nt);

std::vector<Particle> readBinaryParticlesOneSortSTL(FILE *f, particle_sorts sort,int nt);

int readBinaryParticlesAllSorts(FILE *f,int nt,
			                          std::vector<Particle> & ion_vp,
                                      std::vector<Particle> & el_vp,
                                      std::vector<Particle> & beam_vp);

int LoadParticleData(int nt,
		               std::vector<Particle> & ion_vp,
		               std::vector<Particle> & el_vp,
		               std::vector<Particle> & beam_vp, int nx,int ny,int nz);



#endif /* LOAD_DATA_H_ */

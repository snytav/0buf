#ifndef  CELL_BASE_H
#define CELL_BASE_H

#include "particle.h"

template<class Particle>
class CellZBase
{
public:
  int     k;
  double hz;
  double z0;
  double zm;
  int Nz;

//#ifdef __CUDACC__
//__host__ __device__
//#endif
//  CellZbase(){}
//
//#ifdef __CUDACC__
//__host__ __device__
//#endif
//  ~CellZbase(){}


};

#endif

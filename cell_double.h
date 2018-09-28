/*
 * cell_double.h
 *
 *  Created on: May 28, 2016
 *      Author: snytav
 */

#ifndef CELL_DOUBLE_H_
#define CELL_DOUBLE_H_

#define CellExtent 5

class CellDouble {


    double M[CellExtent][CellExtent][CellExtent];
public:


#ifdef __CUDACC__
__host__ __device__
#endif
    double get(int i,int l,int k){return M[i][l][k];}
#ifdef __CUDACC__
__host__ __device__
#endif
    void put(int i,int l,int k,double t){M[i][l][k] = t;}
#ifdef __CUDACC__
__host__ __device__
#endif
    double *getp(int i,int l,int k){return &(M[i][l][k]);}
#ifdef __CUDACC__
__host__ __device__
#endif
    double *getMp(){return ((double *)M);}
};


#endif /* CELL_DOUBLE_H_ */

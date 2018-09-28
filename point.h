#ifndef POINT_H_
#define POINT_H_

template<typename T,int dims,int unity>
class Point{
  T Z;
  public:

     T x,y;
#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
     T z(){return ((dims == 3)*Z + (dims == 2)*((T)unity));}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
     void setZ(T zz){Z = zz;}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
     int dimz2(){return ((dims == 3)*(Z+2) + (dims == 2)*((T)unity));}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
     int dim3D(){return (dims == 3);}

#ifdef __CUDACC__
__host__ __device__ __forceinline__
#endif
     T size2(){return (((int)x+2)*((int)y+2)*((int)Z+2));}

};

#endif

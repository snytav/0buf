/*
 * tensor.h
 *
 *  Created on: Jul 10, 2018
 *      Author: snytav
 */

#ifndef TENSOR_H_
#define TENSOR_H_

#ifndef __CUDACC__
#define __host__
#define __device__
#define __forceinline__
#endif

class LiteCurrentTensorComponent {
public:
	char i11, i12, i13,
	 i21, i22, i23,
	 i31, i32, i33,
	 i41, i42, i43;
};


class CurrentTensorComponent {
public:
	char i11, i12, i13,
	 i21, i22, i23,
	 i31, i32, i33,
	 i41, i42, i43;
	double t[4];


	__host__ __device__	CurrentTensorComponent & operator=(CurrentTensorComponent b)
	{

		i11 = b.i11;
		i12 = b.i12;
		i13 = b.i13;

		i21 = b.i21;
		i22 = b.i22;
		i23 = b.i23;

		i31 = b.i31;
		i32 = b.i32;
		i33 = b.i33;

		i41 = b.i41;
		i42 = b.i42;
		i43 = b.i43;

		t[0] = b.t[0];
		t[1] = b.t[1];
		t[2] = b.t[2];
		t[3] = b.t[3];

		return *this;
	}

};

class LiteCurrentTensor {
public:
	LiteCurrentTensorComponent Jx,Jy,Jz;
};

class CurrentTensor {
public:
	CurrentTensorComponent Jx,Jy,Jz;


	__host__ __device__	CurrentTensor & operator=(CurrentTensor b)
	{
		Jx = b.Jx;
		Jy = b.Jy;
		Jz = b.Jz;

		return *this;
	}
};

class DoubleCurrentTensor {
public:
	CurrentTensor t1,t2;

__host__ __device__	DoubleCurrentTensor & operator=(DoubleCurrentTensor b)
	{

		t1 = b.t1;
		t2 = b.t2;

		return *this;
	}
};



#endif /* TENSOR_H_ */

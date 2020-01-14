#ifndef MAT_3_H
#define MAT_3_H

#include <iostream>

#include "vec3d.h"

// Matrix 4x4 Implementation - 
// Template and or Inline Def in Header. 

// Row Major Order. Multipcation with Vector is Post Multiplicative. v * M = v. Thus no M */*= v overloads here. 
// Matrix is Explcitlly Stored in 1D Array form, and mapped to 2D Indices using i + N * j, where N = 4. 
// Iterate in Row Major (j,i) order. 

template <class T>
class matrix_4x4
{
	// Constructors 
	matrix_4x4() { clear(); }
	matrix_4x4() (T aa, T ab, T ac, T ad, T ba, T bb, T bc, T bd, T ca, T cb, T cc, T cd, T da, T db, T dc, T dd) :
		comp[0] (aa), comp[1](ab), comp[2](ac), comp[3](ad), comp[4](ba), comp[5](bb), comp[0](bc), comp[0](bd), \
		comp[0](ca), comp[0](cb), comp[0](cc), comp[0](cd), comp[0](da), comp[0](db), comp[0](db), comp[0](dd) {}

	explicit matrix_4x4(T* val)
	{
		for (int i = 0; i < size; i++)
		{
			comp[i] = val[i];
		}
	};

	matrix_4x4() (const matrix_4x4 &copy)
	{
		memcpy(comp, copy.comp, (sizeof(T) * size));
	}

	~matrix_4x4() {};

	// Copy Assign 
	matrix_4x4& operator= (const matrix_4x4 &copy)
	{
		if (&copy == this) return *this; 
		memcpy(comp, copy.comp, (sizeof(T) * size));
	}

	// Utils
	inline void clear();
	inline matrix_4x4& transpose();

	inline matrix_4x4& operator[] (float aa, float ab, float ac, float ad, float ba, float bb, float bc, float bd, float ca, float cb, float cc, float cd, float da, float db, float dc, float dd);

	// Indexers
	static inline int idx2Dto1D(int i, int j);
	static inline vec2<int> idx1Dto2D(int i);

	// Math Operators - 

	// MAT A +/- MAT B = MAT C --> Matrix C. (Return New Instance, dont Modifiy this)  
	inline matrix_4x4 operator+ (const matrix_4x4 &b) const; 
	inline matrix_4x4 operator- (const matrix_4x4 &b) const;
	inline matrix_4x4 operator* (const matrix_4x4 &b) const; 

	// Matrix Arithmetic (Modifiy and return ref to this)
	inline matrix_4x4& operator+= (const matrix_4x4 &b);
	inline matrix_4x4& operator+= (const T);
	inline matrix_4x4& operator-= (const matrix_4x4 &b);
	inline matrix_4x4& operator-= (const T);
	inline matrix_4x4& operator*= (const matrix_4x4 &b);
	inline matrix_4x4& operator*= (const T);
	inline matrix_4x4& operator*= (const matrix_4x4 &b);
	inline matrix_4x4& operator*= (const T);


	// Matrix Members - Default Inclass Initalized. 
	const std::size_t size = 16; 
	T comp[size] = {0.0f}; // Matrix Component Array. 
	T *start = comp; T *end = start + (size - 1); // Start/End Pointers. 
};

#endif
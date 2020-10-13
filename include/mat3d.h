#ifndef MAT_3_H
#define MAT_3_H

#include <iostream>
#include <cstring>

#include "vec3d.h"
#include "vec2d.h"


/*
	Matrix 4x4 Implementation, Template and or Inline Def in Header. 
	Row Major Order. Multipcation with Vector is Post Multiplicative. v * M = v. (Transpose for OGL)
	Matrix is Explcitlly Stored in 1D Array form, and mapped to 2D Indices using i + N * j, where N = 4. Iterate in Row Major (j,i) order. 
*/

template <class T>
class matrix_4x4
{
public:
	// Ctors
	matrix_4x4() { ident(); }

	matrix_4x4(T aa, T ab, T ac, T ad, T ba, T bb, T bc, T bd, T ca, T cb, T cc, T cd, T da, T db, T dc, T dd)
	{
		comp[0] = aa, comp[1] = ab, comp[2] = ac, comp[3] = ad, comp[4] = ba, comp[5] = bb, comp[6] = bc, comp[7] = bd,
			comp[8] = ca, comp[9] = cb, comp[10] = cc, comp[11] = cd, comp[12] = da, comp[13] = db, comp[14] = dc, comp[15] = dd;
	}

	explicit matrix_4x4(T* val)
	{
		for (int i = 0; i < m_size; i++)
		{
			comp[i] = val[i];
		}
	};

	// Copy Ctor
	matrix_4x4(const matrix_4x4 &copy)
	{
        std::memcpy(comp, copy.comp, (sizeof(T) * m_size));
	}

	// Copy Assign 
	matrix_4x4& operator= (const matrix_4x4 &copy)
	{
		if (&copy == this) return *this;
        std::memcpy(comp, copy.comp, (sizeof(T) * m_size));
	}

	~matrix_4x4() {};

	// Utils
	inline void clear();
	inline void print_mat();

	inline matrix_4x4& ident(); 
	inline matrix_4x4 transpose();

	static inline float degtoRad(float deg);
	static inline float radtoDeg(float rad);

	// Indexers
	static inline int idx2Dto1D(int i, int j);
	static inline vec2<int> idx1Dto2D(int i);

	// Math Operators - 

	// MAT A +/- MAT B = MAT C --> Matrix C. (Return New Instance, doesnt modifiy this)  
	inline matrix_4x4 operator+ (const matrix_4x4 &b) const;
	inline matrix_4x4 operator- (const matrix_4x4 &b) const;
	inline matrix_4x4 operator* (const matrix_4x4 &b) const;

	// Matrix Arithmetic (Modifiy and return ref to this)
	inline matrix_4x4& operator+= (const matrix_4x4 &b);
	inline matrix_4x4& operator-= (const matrix_4x4 &b);

	inline matrix_4x4& operator+= (T s);
	inline matrix_4x4& operator-= (T s);
	inline matrix_4x4& operator*= (T s);
	inline matrix_4x4& operator/= (T s);

	// Transformation Operations (Modifiy and return this ref)
	inline matrix_4x4& translate(const vec3<T> &tv);
	matrix_4x4& rotate(const vec3<T> &axis, T angle);
	inline matrix_4x4& scale(const vec3<T> &sv);

	// Matrix Transform Operation Builder MFs (return new instance) - 
	static matrix_4x4 make_rotate(const vec3<T> &axis, T angle);

	// Matrix Camera Builder MFs - 
	static inline matrix_4x4 make_lookAt(const vec3<T> &cam_P, const vec3<T> &look_P, const vec3<T> &up);
	static inline matrix_4x4 make_perspective(T fov, T ar, T near, T far);

	// Matrix Members - DICI Default Inclass Initalized. 
	const static std::size_t m_size = 16;
	T comp[m_size] = { 0.0f }; // Matrix Component Array. 
	T *start = comp; T *end = start + (m_size - 1); // Start/End Pointers. 
	const char* label = nullptr; 
};

// -------------------------------------------------

// Templates/Inline MFs Implemenation

// matrix_4x4 Utils

template <class T>
inline void matrix_4x4<T>::clear()
{
	for (int i = 0; i < m_size; i++) { comp[i] = T(0); }
}

template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::ident()
{
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			if (i == j)
			{
				comp[idx2Dto1D(i, j)] = (T)1;
			}
			else
			{
				comp[idx2Dto1D(i, j)] = (T)0;
			}
		}
	}
	return *this; 
}

// Return Transpose of Current Matrix. Does not Modifiy this matrix elements. 
template <class T>
inline matrix_4x4<T> matrix_4x4<T>::transpose()
{
	matrix_4x4<T> tp; 

	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			tp.comp[idx2Dto1D(i, j)] = this->comp[idx2Dto1D(j, i)];
		}
	}
	return tp; 
}

template <class T>
inline void matrix_4x4<T>::print_mat()
{
	std::cout << "\n--------------------\nDEBUG::MATRIX OUPUT BEGIN : \n";
	if (label != nullptr) std::cout << "MATRIX = " << label << "\n";
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			//std::cout << idx2Dto1D(i, j) << "\n";
			std::cout << comp[idx2Dto1D(i, j)] << "|";
		}
		std::cout << "\n";
	}
	std::cout << "\nDEBUG::MATRIX OUPUT END. \n--------------------\n";
}

// matrix_4x4 Math Operators

// Matrix C = A Operator B. (Dont Modifiy Input Matrices). 
template <class T>
inline matrix_4x4<T> matrix_4x4<T>::operator+ (const matrix_4x4<T> &b) const
{
	matrix_4x4<T> res;
	for (int i = 0; i < 16; i++)
	{
		res.comp[i] = (T)(this->comp[i] + b.comp[i]);
	}
	return res;
}

// Matrix C = Matrix A - Matrix B
template <class T>
inline matrix_4x4<T> matrix_4x4<T>::operator- (const matrix_4x4<T> &b) const
{
	matrix_4x4<T> res;
	for (int i = 0; i < 16; i++)
	{
		res.comp[i] = (T)(this->comp[i] - b.comp[i]);
	}
	return res;
}

// Matrix C = Matrix A * Matrix B (DP of this*b) in RM Order. Return New Result Matrix Instance.
template <class T>
inline matrix_4x4<T> matrix_4x4<T>::operator* (const matrix_4x4<T> &b) const
{
	matrix_4x4<T> res;

	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			res.comp[idx2Dto1D(i, j)] = 0.0f; // Esnure Res Elements are 0 before there mult res is calculated.
			for (int k = 0; k < 4; k++)
			{
				// RM Mult and Additon of each A,B Element. 
				res.comp[idx2Dto1D(i, j)] += this->comp[idx2Dto1D(i, k)] * b.comp[idx2Dto1D(k, j)];
			}
		}
	}

	return res;
}

// this Matrix Addition - 
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator+= (const matrix_4x4<T> &b)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] += b.comp[i];
	}

	return *this;
}
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator+= (T s)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] += s;
	}

	return *this;
}

// this Matrix Subtraction -
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator-= (const matrix_4x4<T> &b)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] -= b.comp[i];
	}

	return *this;
}
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator-= (T s)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] -= s;
	}

	return *this;
}


// this Matrix Multipcation With Scalar s. Modifiy and Return this ref. 
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator*= (T s)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] *= s;
	}

	return *this;
}

// Divide Matrix Componets by Scalar. 
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator/= (T s)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] /= s;
	}

	return *this;
}

// Matrix Transformations - 
// Assume Matrix (this) is in Identity State. 
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::translate(const vec3<T> &tv)
{
	comp[3] += tv.x, comp[7] += tv.y, comp[11] += tv.z;
	return *this; 
}

// Modifies this Matrix to Rotate some Vector by some AngleAxis Rotation. This Result is multiplied ontop of current this matrix values.
// See Make_Rotate() Static MF to build and return a new rotation matrix instead.
template <class T>
matrix_4x4<T>& matrix_4x4<T>::rotate(const vec3<T> &axis, T angle)
{
	// X Rot
	T x_ang = axis.x * angle; 
	matrix_4x4<T> x_r(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, cosf(x_ang), sinf(x_ang), 0.0f, 0.0f, -sinf(x_ang), cosf(x_ang), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	// Y Rot
	T y_ang = axis.y * angle;
	matrix_4x4<T> y_r(cosf(y_ang), 0.0f, -sinf(y_ang), 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, sinf(y_ang), 0.0f, cosf(y_ang), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	// Z Rot
	T z_ang = axis.z * angle;
	matrix_4x4<T> z_r(cosf(z_ang), sinf(z_ang), 0.0f, 0.0f, -sinf(z_ang), cosf(z_ang), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);

	// x*y*z Rot Order. 
	matrix_4x4<T> xy = x_r * y_r; 
	matrix_4x4<T> frot = xy * z_r; 
	return *this = *this * frot; 
}

template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::scale(const vec3<T> &sv)
{
	comp[0] *= sv.x, comp[5] *= sv.y, comp[10] *= sv.z;
	return *this; 
}

// Matrix Camera Operations - 

// LookAt Static, Return new mat_4x4<T> with Lookat Basis Vectors - 
template <class T>
inline matrix_4x4<T> matrix_4x4<T>::make_lookAt(const vec3<T> &cam_P, const vec3<T> &look_P, const vec3<T> &up)
{
	vec3<T> lp_n = look_P, cp_n = cam_P;
	vec3<T> zz = cp_n.normalize() - lp_n.normalize(); // Forw
	vec3<T> xx = vec3<T>::cross(zz, up); // Left
	vec3<T> yy = vec3<T>::cross(xx.normalize(), zz.normalize()); // Up

	//return matrix_4x4<T>(xx.x, xx.y, xx.z, 0.0f, yy.x, yy.y, yy.z, 0.0f, zz.x, zz.y, zz.z, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f); 
	return matrix_4x4<T>(xx.x, xx.y, xx.z, 0.0f, yy.x, yy.y, yy.z, 0.0f, zz.x, zz.y, zz.z, 0.0f, cam_P.x, cam_P.y, cam_P.z, 1.0f); // With Cam/Eye Translation. 
}

// Build and Return a new Perspective Projection Matrix Instance, based on Parms. 
template <class T>
inline matrix_4x4<T> matrix_4x4<T>::make_perspective(T fov, T ar, T near, T far)
{
	float t = 1.0f / tanf(fov * 0.5f * (PI / 180.0f));
	float a = (far + near) / (far - near);
	float b = 2.0f*far*near / (far - near);
	float c = t / ar; // Div by AR if non square. 

	//matrix_4x4<T> M_p(c, 0.0f, 0.0f, 0.0f, 0.0f, t, 0.0f, 0.0f, 0.0f, 0.0f, -a, -b, 0.0f, 0.0f, -1.0f, 0.0f); // CM
	matrix_4x4<T> M_p(c, 0.0f, 0.0f, 0.0f, 0.0f, t, 0.0f, 0.0f, 0.0f, 0.0f, -a, -1.0f, 0.0f, 0.0f, -b, 0.0f);

	return M_p; 
}

// Rotation Matrix Builder Static MF (Return Created RotMat, oppose to multiplying with this, as matrix_4x4<T>::rotate() does)
template<class T>
matrix_4x4<T> matrix_4x4<T>::make_rotate(const vec3<T> &axis, T angle)
{
	// X Rot
	T x_ang = axis.x * angle;
	matrix_4x4<T> x_r(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, cosf(x_ang), sinf(x_ang), 0.0f, 0.0f, -sinf(x_ang), cosf(x_ang), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	// Y Rot
	T y_ang = axis.y * angle;
	matrix_4x4<T> y_r(cosf(y_ang), 0.0f, -sinf(y_ang), 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, sinf(y_ang), 0.0f, cosf(y_ang), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	// Z Rot
	T z_ang = axis.z * angle;
	matrix_4x4<T> z_r(cosf(z_ang), sinf(z_ang), 0.0f, 0.0f, -sinf(z_ang), cosf(z_ang), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);

	// x*y*z Rot Order. 
	matrix_4x4<T> xy = x_r * y_r;
	matrix_4x4<T> frot = xy * z_r;

	return frot; 
}

// Static MFs - Indexers - 

template <class T>
inline int matrix_4x4<T>::idx2Dto1D(int i, int j)
{
	return (int)i + 4 * j;
}

template <class T>
inline vec2<int> matrix_4x4<T>::idx1Dto2D(int i)
{
	int ii = i % 4;
	int jj = i / 4;
	return vec2<int>(ii, jj);
}


// Matrix_4x4 Static MFs - 

// Angle Conversion - Explcitlly uses <T> == float. 

// Degrees to Radians
template <class T>
inline float matrix_4x4<T>::degtoRad(float deg)
{
	return (float)deg * (PI / 180.0f);
}

// Radians to Degrees 
template <class T>
inline float matrix_4x4<T>::radtoDeg(float rad)
{
	return (float)rad * (180.0f / PI);
}


#endif

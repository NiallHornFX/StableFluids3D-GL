#ifndef MAT_3_H
#define MAT_3_H

#include <iostream>

/*	Matrix Info -
	Matrix 4x4 Implementation, Template and or Inline Def in Header. 
	Row Major Order. Multipcation with Vector is Post Multiplicative. v * M = v. Thus no M * / *= v overloads here. 
	Matrix is Explcitlly Stored in 1D Array form, and mapped to 2D Indices using i + N * j, where N = 4. Iterate in Row Major (j,i) order. 
*/

template <class T>
class matrix_4x4
{
public:
	// Constructors 
	matrix_4x4() { clear(); }

	matrix_4x4(T aa, T ab, T ac, T ad, T ba, T bb, T bc, T bd, T ca, T cb, T cc, T cd, T da, T db, T dc, T dd)
	{
		comp[0] = aa, comp[1] = ab, comp[2] = ac, comp[3] = ad, comp[4] = ba, comp[5] = bb, comp[6] = bc, comp[7] = bd,
			comp[8] = ca, comp[9] = cb, comp[10] = cc, comp[11] = cd, comp[12] = da, comp[13] = db, comp[14] = dc, comp[15] = dd;
	}

	explicit matrix_4x4(T* val)
	{
		for (int i = 0; i < size; i++)
		{
			comp[i] = val[i];
		}
	};

	// Copy Ctor
	matrix_4x4(const matrix_4x4 &copy)
	{
		memcpy(comp, copy.comp, (sizeof(T) * m_size));
	}

	// Copy Assign 
	matrix_4x4& operator= (const matrix_4x4 &copy)
	{
		if (&copy == this) return *this;
		memcpy(comp, copy.comp, (sizeof(T) * size));
	}

	~matrix_4x4() {};

	// Utils
	inline void clear();
	inline matrix_4x4& transpose();
	inline void print_mat();

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
	inline matrix_4x4& operator+= (const T s);
	inline matrix_4x4& operator-= (const matrix_4x4 &b);
	inline matrix_4x4& operator-= (const T s);
	inline matrix_4x4& operator*= (const matrix_4x4 &b);
	inline matrix_4x4& operator*= (const T s);
	inline matrix_4x4& operator/= (const matrix_4x4 &b);
	inline matrix_4x4& operator/= (const T s);


	// Matrix Members - Default Inclass Initalized. 
	const static std::size_t m_size = 16;
	T comp[m_size] = { 0.0f }; // Matrix Component Array. 
	T *start = comp; T *end = start + (m_size - 1); // Start/End Pointers. 
};



// Templates/Inline MFs Implemenation \\

// matrix_4x4 Utils \\

template <class T>
inline void matrix_4x4<T>::clear()
{
	for (int i = 0; i < m_size; i++) { comp[i] = 0.0f; }
}

template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::transpose()
{
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			comp[idx2Dto1D(i, j)] = comp[idx2Dto1D(j, i)];
		}
	}

	return *this;
}

template <class T>
inline void matrix_4x4<T>::print_mat()
{
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			std::cout << comp[res.idx2Dto1D(i, j)] << "|";
		}
		std::cout << "\n";
	}
}

// matrix_4x4 Math Operators \\

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

// Matrix A x Matrix B (DP of RxC)
template <class T>
inline matrix_4x4<T> matrix_4x4<T>::operator* (const matrix_4x4<T> &b) const
{
	matrix_4x4<T> res;

	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int k = 0; k < 4; k++)
			{
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
inline matrix_4x4<T>& matrix_4x4<T>::operator+= (const T s)
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
inline matrix_4x4<T>& matrix_4x4<T>::operator-= (const T s)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] -= s;
	}

	return *this;
}

// this Matrix Multipcation - 

template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator*= (const matrix_4x4<T> &b)
{
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int k = 0; k < 4; k++)
			{
				this->comp[idx2Dto1D(i, j)] += this->comp[idx2Dto1D(i, k)] * b.comp[idx2Dto1D(k, j)];
			}
		}
	}

	return *this;

}
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator*= (const T s)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] *= s;
	}

	return *this;
}

// this Matrix Division - 

template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator/= (const matrix_4x4<T> &b)
{

}
template <class T>
inline matrix_4x4<T>& matrix_4x4<T>::operator/= (const T s)
{
	for (int i = 0; i < 16; i++)
	{
		comp[i] /= s;
	}

	return *this;
}

// Static MFs - Indexers - 

template <class T>
inline int matrix_4x4<T>::idx2Dto1D(int i, int j)
{
	return (int)i * 4 + j;
}

template <class T>
inline vec2<int> matrix_4x4<T>::idx1Dto2D(int i)
{
	int ii = i % 4;
	int jj = i / 4;
	return vec2<int>(ii, jj);
}


#endif
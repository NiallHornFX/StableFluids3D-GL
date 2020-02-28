#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

extern const double PI; 

template <class T>
class vec3
{
public:
	vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
	vec3() : x(0), y(0), z(0) {}
	~vec3() {}

	// Copy 
	vec3(const vec3 &copy) : x(copy.x), y(copy.y), z(copy.z) {}
	vec3& operator= (const vec3 &copy)
	{
		if (&copy == this) return *this; 
		this->x = copy.x, this->y = copy.y, this->z = copy.z;
		return *this; 
	}

	// vec3 Vector Math Operations, Return Scalar. 
	inline float length() const; // 
	inline float dot(const vec3 &b) const; 
	inline float angle(const vec3 &b) const; 

	// vec3 Vector Math Operations, Return modified *this. 
	inline vec3& normalize(); 
	inline vec3& clear(); 

	// LHS vec3 Arithmetic Overloads - Return new vec3 with operation of this and B vec3 completed. Do Not Modifiy this.
	inline vec3 operator+ (const vec3 &b) const; 
	inline vec3 operator- (const vec3 &b) const;

	inline float operator*(const vec3 &b) const; 

	// LHS Operator Overloads vec3 - Modifiy this return *this. 
	inline vec3& operator+= (T scalar); 
	inline vec3& operator+= (const vec3 &b); 
	inline vec3& operator-= (T scalar); 
	inline vec3& operator-= (const vec3 &b);
	inline vec3& operator/= (T scalar); 
	inline vec3& operator/= (const vec3 &b); 
	inline vec3& operator*= (T scalar); 

	// Cross Product Return new vec3<T> c = a x b
	static inline vec3 cross(const vec3 &a, const vec3 &b);

	// Angle Conversion
	static inline float degtoRad(float deg);
	static inline float radtoDeg(float rad);

	T x, y, z; 
};


// vec3 Implementation \\

// _vec3 Math Operations \\

// length = Sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z))
template <class T>
inline float vec3<T>::length() const
{
	return std::sqrtf(std::powf(x, 2.0) + std::powf(y, 2.0) + std::powf(z, 2.0));
}

// vec3 Dot Product 
// dot = (this->x * b.x) + (this->y * b.y) + (this->z * b.z)
template <class T>
inline float vec3<T>::dot(const vec3 &b) const
{
	return (float)((x * b.x) + (y * b.y) + (z * b.z));
}
// dot v1*v2 overload = dot = (this->x * b.x) + (this->y * b.y) + (this->z * b.y)
template <class T>
inline float vec3<T>::operator*(const vec3 &b) const
{
	dot(b);
}

// vec3 Angle - Angle Between this & b (Radians) 
// theta = acos(dot(normalize(this), normalize(b)))
template <class T>
inline float vec3<T>::angle(const vec3 &b) const
{
	vec3 v1 = *this, v2 = b; 

	float theta = std::acosf((v1.normalize()).dot(v2.normalize())); 
	return theta;
}

// vec3 Normalize 
// this / ||this|| --> *this
template <class T>
inline vec3<T>& vec3<T>::normalize()
{
	float l = length();
	x = x / l;
	y = y / l;
	z = z / l;

	return *this; 
}

// vec3 Clear
// this->x = 0, this->y = 0, this->z = 0 --> *this
template <class T>
inline vec3<T>& vec3<T>::clear()
{
	x = (T) 0, y = (T) 0, z = (T) 0;
	return *this; 
}


// _LHS vec3 Operand overloads \\

// Vector Additon and Subtraction with LHS and RHS Vector. Return new resulting vec3. 
// v3 = this + v2
template <class T>
inline vec3<T> vec3<T>::operator+(const vec3 &b) const
{
	T n_x = x + b.x;
	T n_y = y + b.y;
	T n_z = z + b.z;

	return vec3<T>(n_x, n_y, n_z);
}

// v3 = v2 - this
template <class T>
inline vec3<T> vec3<T>::operator-(const vec3 &b) const
{
	T n_x = x - b.x;
	T n_y = y - b.y;
	T n_z = z - b.z;

	return vec3<T>(n_x, n_y, n_z);
}


// _LHS vec3<T> RHS <T> Scalar,Vec3 --> Return *this modified vec3<T> \\

// vec3 Addition by Scalar or vec3
// this->x += scalar, this->y += scalar, this->z += scalar ---> *this
template <class T>
inline vec3<T>& vec3<T>::operator+=(T scalar)
{
	x += scalar, y += scalar, z += scalar;
	return *this;
}
template <class T>
inline vec3<T>& vec3<T>::operator+=(const vec3 &b)
{
	x += b.x, y += b.y, z += b.z;
	return *this;
}

// vec3 Subtraction by Scalar or vec3 
// this->x -= scalar, this->y -= scalar, this->z -= scalar ---> *this
template <class T>
inline vec3<T>& vec3<T>::operator-=(T scalar)
{
	x -= scalar, y -= scalar, z -= scalar;
	return *this;
}
template <class T>
inline vec3<T>& vec3<T>::operator-=(const vec3 &b)
{
	x -= b.x, y -= b.y, z -= b.z;
	return *this;
}

// vec3 Multiplicaton by Scalar
// this->x *= scalar, this->y *= scalar, this->z *= scalar ---> *this
template <class T>
inline vec3<T>& vec3<T>::operator*=(T scalar)
{
	x *= scalar, y *= scalar, z *= scalar;
	return *this;
}

// vec3 Division by Scalar or vec3 
// this->x /= scalar, this->y /= scalar, this->z /= scalar ---> *this
template <class T>
vec3<T>& vec3<T>::operator/=(T scalar)
{
	x /= scalar, y /= scalar, z/= scalar;
	return *this;
}
template <class T>
vec3<T>& vec3<T>::operator/=(const vec3 &b)
{
	x /= b.x, y /= b.y, z /= b.z;
	return *this;
}

// Static Cross Product --> c = a x b
template <class T>
inline vec3<T> vec3<T>::cross(const vec3 &a, const vec3 &b)
{
	float cx = (a.y * b.z) - (a.z * b.y);
	float cy = (a.z * b.x) - (a.x * b.z);
	float cz = (a.x * b.y) - (a.y * b.x);
	return vec3<T>(cx, cy, cz);
}

// Static vec3 (Utility) Member Functions \\

// Degrees to Radians
template <class T>
inline float vec3<T>::degtoRad(float deg)
{
	return (float)deg * (PI / 180.0f);
}

// Radians to Degrees 
template <class T>
inline float vec3<T>::radtoDeg(float rad)
{
	return (float)rad * (180.0f / PI);
}


// _RHS vec3 Operand Global Operator Overload Free Functions \\ 

// s1 * v1 ---> v1;
template <class T>
inline vec3<T> operator* (const float mult, const vec3<T> &vec)
{
	return (vec.x *= mult, vec.y *= mult, vec.z *= mult);
}

// s1 + v1 ---> v1;
template <class T>
inline vec3<T> operator+ (const float add, const vec3<T> &vec)
{
	return (vec.x + add, vec.y + add, vec.z + add); 
}

#endif
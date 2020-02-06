#ifndef VEC2_H
#define VEC2_H

extern const double PI;

template <class T>
class vec2
{
public:
	vec2(T xx, T yy) : x(xx), y(yy) {};
	vec2() : x(0), y(0) {};
	~vec2() {};

	vec2(const vec2 &copy) : x(copy.x), y(copy.y) {}
	vec2& operator= (const vec2 &copy)
	{
		if (&copy == this) return *this;
		this->x = copy.x, this->y = copy.y;
		return *this;
	}

	// Vec2 Vector Math Operations, Return Scalar. 
	inline float length() const; 
	inline float dot(const vec2 &b) const; 
	inline float angle(const vec2 &b) const; 

	// Vec2 Vector Math Operations, Return modified *this. 
	inline vec2& normalize();
	inline vec2& clear(); 

	// LHS Vec2 Arithmetic Overloads - Return new vec2 with operation of this and B vec2 completed. Do Not Modifiy this.
	inline vec2 operator+ (const vec2 &b) const;
	inline vec2 operator- (const vec2 &b) const; 

	inline float operator*(const vec2 &b) const; 

	// LHS Operator Overloads Vec2 - Modifiy this return *this. 
	inline vec2& operator+= (T scalar);
	inline vec2& operator+= (const vec2 &b);
	inline vec2& operator-= (T scalar);
	inline vec2& operator-= (const vec2 &b);
	inline vec2& operator/= (T scalar); 
	inline vec2& operator/= (const vec2 &b);
	inline vec2& operator*= (T scalar);

	static inline float degtoRad(float deg);
	static inline float radtoDeg(float rad);

	T x, y;
};


// vec2 Implementation \\

// _vec2 Math Operations \\

// length = Sqrt((this->x * this->x) + (this->y * this->y))
template <class T>
inline float vec2<T>::length() const
{
	return std::sqrtf(std::powf(x, 2.0) + std::powf(y, 2.0));
}

// vec2 Dot Product 
// dot = (this->x * b.x) + (this->y * b.y)
template <class T>
inline float vec2<T>::dot(const vec2 &b) const
{
	return (float)(x * b.x) + (y * b.y);
}
// dot v1*v2 overload = dot = (this->x * b.x) + (this->y * b.y)
template <class T>
inline float vec2<T>::operator*(const vec2 &b) const
{
	dot(b);
}

// vec2 Angle - Angle Between this & b (Radians) 
// theta = acos(dot(normalize(this), normalize(b)))
template <class T>
inline float vec2<T>::angle(const vec2 &b) const
{
	vec2 v1 = *this, v2 = b;

	float theta = std::acosf((v1.normalize()).dot(v2.normalize()));
	return theta;
}

// vec2 Normalize 
// this / ||this|| --> *this
template <class T>
inline vec2<T>& vec2<T>::normalize()
{
	float l = length();
	x = x / l;
	y = y / l;

	return *this;
}

// vec2 Clear
// this->x = 0, this->y = 0 --> *this
template <class T>
inline vec2<T>& vec2<T>::clear()
{
	x = (T)0, y = (T)0;
	return *this;
}


// _LHS vec2 Operand overloads \\

// Vector Additon and Subtraction with LHS and RHS Vector. Return new resulting vec2. 
// v3 = this + v2
template <class T>
inline vec2<T> vec2<T>::operator+(const vec2 &b) const
{
	T n_x = x + b.x;
	T n_y = y + b.y;

	return vec2(n_x, n_y);
}

// v3 = v2 - this
template <class T>
inline vec2<T> vec2<T>::operator-(const vec2 &b) const
{
	T n_x = x - b.x;
	T n_y = y - b.y;

	return vec2(n_x, n_y);
}


// _LHS vec2<T> RHS <T> Scalar,Vec2 --> Return *this modified vec2<T> \\

// vec2 Addition by Scalar or vec2
// this->x += scalar, this->y += scalar ---> *this
template <class T>
inline vec2<T>& vec2<T>::operator+=(T scalar)
{
	x += scalar, y += scalar;
	return *this;
}
template <class T>
inline vec2<T>& vec2<T>::operator+=(const vec2 &b)
{
	x += b.x, y += b.y;
	return *this;
}

// vec2 Subtraction by Scalar or vec2 
// this->x -= scalar, this->y -= scalar ---> *this
template <class T>
inline vec2<T>& vec2<T>::operator-=(T scalar)
{
	x -= scalar, y -= scalar;
	return *this;
}
template <class T>
inline vec2<T>& vec2<T>::operator-=(const vec2 &b)
{
	x -= b.x, y -= b.y;
	return *this;
}

// vec2 Multiplicaton by Scalar
// this->x *= scalar, this->y *= scalar ---> *this
template <class T>
inline vec2<T>& vec2<T>::operator*=(T scalar)
{
	x *= scalar, y *= scalar;
	return *this;
}

// vec2 Division by Scalar or vec2 
// this->x /= scalar, this->y /= scalar ---> *this
template <class T>
vec2<T>& vec2<T>::operator/=(T scalar)
{
	x /= scalar, y /= scalar;
	return *this;
}
template <class T>
vec2<T>& vec2<T>::operator/=(const vec2 &b)
{
	x /= b.x, y /= b.y;
	return *this;
}

// Static vec2 (Utility) Member Functions \\

// Degrees to Radians
template <class T>
static inline float vec2<T>::degtoRad(float deg)
{
	return (float)deg * (PI / 180.0f);
}

// Radians to Degrees 
template <class T>
static inline float vec2<T>::radtoDeg(float rad)
{
	return (float)rad * (180.0f / PI);
}


// _RHS vec2 Operand Global Operator Overload Free Functions \\ 

// s1 * v1 ---> v1;
template <class T>
inline vec2<T>& operator* (float mult, vec2<T> &vec)
{
	vec.x *= mult, vec.y *= mult;
	return vec;
}

// s1 + v1 ---> v1;
template <class T>
inline vec2<T>& operator+ (float add, vec2<T> &vec)
{
	vec.x += add, vec.y += add;
	return vec;
}


#endif


#ifndef GRIDS2_H
#define GRIDS2_H

#include <iostream>
#include <vector>
#include <tuple>

// Cannot Forward Decl due to Class TMP. 
#include "vec2d.h"
#include "vec3d.h"
#include "mat3d.h"

// Enums - 

enum interpType
{
	interp_TriLinear = 0,
	interp_TriCosine,
	interp_TriCubic
};

// ABC grid3<T> Class. To be part of derived grid3_scalar + grid3_vector implementations.
// It is Responsible for Grid_Data<T> (vector) construction. 

// Grid3 Abstract Base Class. 
template <class T>
class grid3
{
	friend class renderobject_3D_OGL;
public:
	grid3() = delete;
	grid3(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s);

	virtual ~grid3();

	// !TO-DO Inline Grid Getter/Setters. 
	void setdata(T data, int i);
	void setdata(T data, int i, int j, int k);
	void adddata(T data, int i);
	void adddata(T data, int i, int j, int k);

	T getdata(int i) const;
	T getdata(int i, int j, int k) const;

	// PV Grid Sampler and MinMax -
	virtual T sampler(const vec3<float> &gs_loc, interpType interp) const = 0;
	virtual std::tuple<T, T> minmax(int i0, int j0, int k0) const = 0; 

	// V Utils
	virtual void swap(const grid3 *B);
	virtual void clear() = 0;
	virtual void printinfo() const = 0;

	// Util Member Getters 
	vec3<std::size_t> get_dimsize() const; 
	std::size_t get_edgesize() const; 
 
	// Unsafe Getters 
	virtual std::vector<T>& griddatavector_getter(); 
	virtual T* getdataarray() const;

	// Indexers 
	__forceinline int idx_3Dto1D(int i, int j, int k) const;
	__forceinline vec3<int> idx_1Dto3D(int k) const;

protected:
	std::vector<T> grid_data;
	std::size_t x_size, y_size, z_size, edge_size, total_size;

	matrix_4x4<float> GStoWS; 

};

// Grid3 Scalar Class - 
template <class T>
class grid3_scalar : public grid3<T>
{
	friend class renderobject_3D_OGL;
public:
	grid3_scalar() = delete;
	grid3_scalar(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s);

	virtual ~grid3_scalar() override;

	virtual T sampler(const vec3<float> &gs_loc, interpType interp) const override; 
	virtual std::tuple<T, T> minmax(int i0, int j0, int k0) const override; 

	virtual void clear() override; 
	virtual void printinfo() const override; 
};

// Grid3 Vector Class - 
// vec_t (vec_t t = Template, vec_t<T>T = float). Vector Type is Templated, but Internal Type is explictlly float. 

template <class T>
class grid3_vector : public grid3<T>
{
	friend class renderobject_3D_OGL;
public:
	grid3_vector() = delete;
	grid3_vector(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s);

	virtual ~grid3_vector() override;

	// Vec Component Wise Data Acess (1D, 3D) -
	void setdata_x(float xx, int i);
	void setdata_x(float xx, int i, int j, int k); 
	void setdata_y(float yy, int i);
	void setdata_y(float yy, int i, int j, int k);
	void setdata_z(float zz, int i);
	void setdata_z(float zz, int i, int j, int k);

	void adddata_x(float xx, int i);
	void adddata_x(float xx, int i, int j, int k);
	void adddata_y(float yy, int i);
	void adddata_y(float yy, int i, int j, int k);
	void adddata_z(float zz, int i);
	void adddata_z(float zz, int i, int j, int k);

	float getdata_x(int i) const;
	float getdata_x(int i, int j, int k) const; 
	float getdata_y(int i) const;
	float getdata_y(int i, int j, int k) const; 
	float getdata_z(int i) const;
	float getdata_z(int i, int j, int k) const;

	virtual T sampler(const vec3<float> &gs_loc, interpType interp) const override; 
	virtual std::tuple<T, T> minmax(int i0, int j0, int k0) const override;

	virtual void clear() override;
	virtual void printinfo() const override;

};

// !TO-DO Move Inlined Defs to grid3d.inl and incl. 

// INLINE GRID MEMBER FUNTIONS \\ 

// Inline (Forced) Indexing - 
template <class T>
vec3<int> grid3<T>::idx_1Dto3D(int i) const
{
	int ii = i / ((y_size + edge_size) * (z_size + edge_size));
	int jj = (i / (z_size + edge_size)) % (y_size + edge_size);
	int kk = i % (z_size + edge_size);

	return vec3<int>(ii, jj, kk);
}

template <class T>
int grid3<T>::idx_3Dto1D(int i, int j, int k) const
{
	return (int) i + (x_size + edge_size) * (j + (z_size + edge_size) * k);
}


#endif
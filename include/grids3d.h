#ifndef GRIDS2_H
#define GRIDS2_H

#include <iostream>
#include <vector>

// Cannot Forward Decl due to Class TMP. 
#include "vec3d.h"
#include "mat3d.h"

// ABC grid3<T> Class. To be part of derived grid3_scalar + grid3_vector implementations.
// It is Responsible for Grid_Data<T> Construction. 

// Grid3 Abstract Base Class. 
template <class T>
class grid3
{
	friend class renderobject_3D_OGL;
public:
	grid3() = delete;
	grid3(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s);

	virtual ~grid3();

	// PVMFs to Override in grid3_[..] derivations. 
	// Data Acess (1D, 3D) -
	void setdata(T data, int i);

	void setdata(T data, int i, int j, int k);
	void adddata(T data, int i);

	void adddata(T data, int i, int j, int k);
	T getdata(int i) const;
	T getdata(int i, int j, int k) const;

	virtual void swap(const grid3 *B);
	virtual void clear() = 0;

	virtual void printinfo() const = 0;
	vec3<std::size_t> get_dimsize() const; 
	std::size_t get_edgesize() const; 
 
	virtual std::vector<T>* griddataptr_getter() const;
	virtual T* getdataarray() const;

	// Non Virtual Indexers 
	int idx_3Dto1D(int i, int j, int k) const;
	vec3<int> idx_1Dto3D(int k) const;

protected:
	std::vector<T> *grid_data;
	std::size_t x_size, y_size, z_size, edge_size, total_size;

	matrix_4x4<float> GStoWS; // Local Grid Space to World Space Transformation Matrix. 

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

	virtual void clear() override; 
	virtual void printinfo() const override; 
};


// Grid3 Vector Class - 
template <class T>
class grid3_vector : public grid3<T>
{
	friend class renderobject_3D_OGL;
public:
	grid3_vector() = delete;
	grid3_vector(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s);

	virtual ~grid3_vector() override;

	// Data Acess (1D, 3D) -
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

	virtual void clear() override;
	virtual void printinfo() const override;

};


#endif
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
	void setdata_x(T xx, int i);
	void setdata_x(T xx, int i, int j, int k); 
	void setdata_y(T yy, int i);
	void setdata_y(T yy, int i, int j, int k); 
	void setdata_z(T zz, int i);
	void setdata_z(T zz, int i, int j, int k);

	void adddata_x(T xx, int i);
	void adddata_x(T xx, int i, int j, int k);
	void adddata_y(T yy, int i);
	void adddata_y(T yy, int i, int j, int k);
	void adddata_z(T zz, int i);
	void adddata_z(T zz, int i, int j, int k);

	T getdata_x(int i) const;
	T getdata_x(int i, int j, int k) const; 
	T getdata_y(int i) const;
	T getdata_y(int i, int j, int k) const; 
	T getdata_z(int i) const;
	T getdata_z(int i, int j, int k) const;

	virtual void clear() override;
	virtual void printinfo() const override;

};







/* !! OLD NON POLYMORPHIC CODE !!

class grid3_scalar
{
friend class renderobject_3D_OGL; // Temp. 

public:
	grid3_scalar() = delete; 
	grid3_scalar(int x_s, int y_s, int z_s, int e_s, int idd, float spc);

	grid3_scalar(grid3_scalar &copy);

	~grid3_scalar();

	// Grid3 MFs - 

	grid3_scalar& operator=(grid3_scalar &copy);

	//Set/Get Data Acess - 
	void setdata(float data, int i);
	void setdata(float data, int i, int j, int k);

	void adddata(float data, int i);
	void adddata(float data, int i, int j, int k);

	float getdata(int i); // 1D Access
	float getdata(int i, int j, int k); 

	// Swap Field Grid_Data With Another Field Passed By Pointer.
	void swap(grid3_scalar *B);

	// Util MFs - 
	void clear();
	void printinfo();
	void printsize();

	// Grid_Data vector and Internal Array Getters - 
	std::vector<float>* griddataptr_getter();
	float* getdataarray();

	void tranpose_gridata(std::vector<float>*ptr); 

	// Indexing - 
	int idx_2Dto1D(int i, int j); 
	vec2<int> idx_1Dto2D(int k);

	// Debug MFs - 
	void print_celldata(); // Do NOT Call this on large Grids ! 
	void disp();

private:

	// Declare Grid Data Raw Ptr, to later define holding Scalar Grid Data as std vector type float. 
	// Size will be determined by Constructor call Dimensions and spacing. 
	std::vector<float> *grid_data;

	int x_size, y_size, z_size, edge_size, total_size; 

	float grid_spacing;

	int ID; // Grid ID Var 
};

// In future vector grid, should be polymorphically derived from grid2_scalar -
// Really makes sense for shared MFs, and members ofc. For now, no polymorhpishm.

class grid2_vector
{
	friend class renderobject_2D_OGL; // Temp 
public:
	grid2_vector() = delete; 
	grid2_vector(int x_s, int y_s, int e_s, int idd, float spc);
	grid2_vector(grid2_vector &copy);

	~grid2_vector();

	// Grid2 MFs - 

	grid2_vector& operator=(grid2_vector &copy);

	//Set/Get Data Acess (With Per Member/Component Acess Incase operations need to be scalar) - 
	// 1D and 2D IDX overloads. 

	// Set Data - vec2 
	void setdata(vec2<float> data, int i);
	void setdata(vec2<float> data, int i, int j); // Keep GridData Private, Acess through MFs. 

	// Set Data - vec2 Components 
	void setdata_x(float xx, int i);
	void setdata_x(float xx, int i, int j); // Set X Component Vec Only at ij 2D Index. 
	void setdata_y(float yy, int i);
	void setdata_y(float yy, int i, int j); // Set Y Component Vec Only, at ij 2D Index.

	// Add Data - vec2
	void adddata(vec2<float> data, int i);
	void adddata(vec2<float> data, int i, int j);

	// Get Data - vec2
	vec2<float> getdata(int i, int j); // Get vec2 data at ij 2D Index. 
	vec2<float> getdata(int i);

	// Get Data - vec2 Components 
	float getdata_x(int i);
	float getdata_x(int i, int j); // Get vec2 float component data at ij 2d index
	float getdata_y(int i);
	float getdata_y(int i, int j); // Get vec2 float component data at ij 2d index

	// Swap Field Grid_Data With Another Field Passed By Pointer.
	void swap(grid2_vector *B);

	// Util MFs - 
	void clear();
	void printinfo();

	// Grid_Data vector and Internal Array Getters - 
	std::vector<vec2<float>>* griddataptr_getter();
	vec2<float>* getdataarray();

	// Indexing - 
	int idx_2Dto1D(int i, int j);
	vec2<int> idx_1Dto2D(int k);

	// Debug MFs - 
	void print_celldata(); // Do NOT Call this on large Grids ! 

private:

	// Declare Grid Data Raw Ptr, to later define holding Scalar Grid Data as std vector type float. 
	// Size will be determined by Constructor call Dimensions and spacing. 
	std::vector<vec2<float>> *grid_data;

	int x_size, y_size, edge_size, total_size;

	float grid_spacing; // VoxelSize Default 1.0f.  

	int ID; // Grid ID Var

};

*/
#endif
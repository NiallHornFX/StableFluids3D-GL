#ifndef GRIDS2_H
#define GRIDS2_H

#include <iostream>
#include <vector>

#include "vec2.h"
#include "rgb.h"

class grid2_scalar
{
friend class renderobject_2D_OGL; // Temp. 

public:
	grid2_scalar() = delete; 
	grid2_scalar(int x_s, int y_s, int e_s, int idd, float spc);

	grid2_scalar(grid2_scalar &copy);

	~grid2_scalar();

	// Grid2 MFs - 

	grid2_scalar& operator=(grid2_scalar &copy);

	//Set/Get Data Acess - 
	// 1D and 2D IDX overloads. 
	void setdata(float data, int i);
	void setdata(float data, int i, int j);

	void adddata(float data, int i);
	void adddata(float data, int i, int j);

	float getdata(int i); // 1D Access
	float getdata(int i, int j); 

	// Swap Field Grid_Data With Another Field Passed By Pointer.
	void swap(grid2_scalar *B);

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

	int x_size, y_size, edge_size, total_size; 

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

class grid2_img
{
public:

	explicit grid2_img(const char *path, int x_s, int y_s, int idd); 
	~grid2_img();

	void setdata(RGB data, int i, int j); // Keep GridData Private, Acess through MFs. 
	void setdata_R(float xx, int i, int j); // Set X Component Vec Only at ij 2D Index. 
	void setdata_G(float yy, int i, int j); // Set Y Component Vec Only, at ij 2D Index.
	void setdata_B(float yy, int i, int j); // Set Y Component Vec Only, at ij 2D Index.

	RGB getdata(int i, int j); // Get vec2 data at ij 2D Index. 
	float getdata_R(int i, int j); // Get vec2 float component data at ij 2d index
	float getdata_G(int i, int j); // Get vec2 float component data at ij 2d index
	float getdata_B(int i, int j); // Get vec2 float component data at ij 2d index

	// Swap Field Grid_Data With Another Field Passed By Pointer.
	void swap(grid2_img *B);

	// Indexing - 
	int idx_2Dto1D(int i, int j);
	vec2<int> idx_1Dto2D(int k);

	// Debug MFuncs.
	void print_rawdata();

private:

	unsigned char *raw_data_bytes; // Currently Unsused. For Storing Pure Binary 8Bit RGB Data. 
	std::string *raw_data; // Keep Raw Data on Heap Always Via Ptr Also.
	std::vector<RGB> *grid_data;
	int x_size;
	int y_size;
	int ID;
	const char *img_path; 
};

#endif
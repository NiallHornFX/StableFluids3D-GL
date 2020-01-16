// Implementation of grids2 

#include "grids3d.h"
#include "vec3d.h"

#include <memory>
#include <fstream>
#include <string>
#include <cassert>

#include <omp.h> 


extern short verbose; // Get Verbose Global from main. 

/*--
Seeing the copied code reminds me to make this more polymorhic on next round, ie derive
vector grid from scalar etc. vec3 from vec2... 
--*/

/*

// grid2_scalar Implementation - \\

// // grid2_scalar implemented Constructor. 
grid2_scalar::grid2_scalar(int x_s, int y_s, int e_s, int idd, float spc) :
	x_size(x_s), y_size(y_s), edge_size(e_s), ID(idd), grid_spacing(spc)
{
	total_size = (x_size + edge_size) * (y_size + edge_size);

	// Alloacte vector<float> Grid Data Via Grid_Data Pointer Member. 
	grid_data = new std::vector<float>(total_size, 0.0f);
}

// grid2_scalar Copy Constructor - 
// Copy Contents from copy grid to this, on initalization call to Cctor. 
// Deep Copy Grid_Data Memory. (Not just ptr).
grid2_scalar::grid2_scalar(grid2_scalar &copy)
	: x_size(copy.x_size), y_size(copy.y_size), edge_size(copy.edge_size), ID(copy.ID), 
	grid_spacing(copy.grid_spacing)
{
	total_size = (x_size + edge_size) * (y_size * edge_size);

	// Copy Grid_Data Directly. 
	grid_data = new std::vector<float>((*copy.grid_data));
}

// grid2_scalar Copy Assignment Operator Overload - 
// Copy Contents from copy grid to this, post initalization. 
// Deep Copy Grid_Data Memory. (Not just ptr).
grid2_scalar&::grid2_scalar::operator=(grid2_scalar &copy)
{
	// Check if this passed. 
	if (&copy == this) return *this; 

	// Dealloc current grid_data vector. 
	delete this->grid_data;

	// Copy grid_data vector from copy grid. 
	grid_data = new std::vector<float>((*copy.grid_data));

	// Copy Local Member Values 
	this->x_size = copy.x_size; this->y_size = copy.y_size; this->edge_size = copy.edge_size; 
	this->total_size = (x_size + edge_size) * (y_size + edge_size);
	this->ID = copy.ID + 5; this->grid_spacing = copy.grid_spacing;

	return *this; 
}

// grid2_scalar Destructor.
grid2_scalar::~grid2_scalar()
{
	// Clear grid_data vector. Set Ptr to null. 
	delete grid_data;
	grid_data = nullptr;

}

// Util Member Function to Print Grid Info 
void grid2_scalar::printinfo()
{
	std::cout << "DEBUG::Grid 2D Scalar " << ID << " Info BEGIN - \n";
	std::cout << "Grid Cell Count = " << total_size << "\n";
	std::cout << "Grid X Row Size = " << x_size << "\n";
	std::cout << "Grid Y Row Size = " << y_size << "\n";
	std::cout << "Grid Edge Cell Size = " << edge_size << "\n";
	std::cout << "DEBUG::Grid 2D Scalar " << ID << " Info END. \n \n";
}

// Member Function to set Grid Data Of Scalar Grid. 
void grid2_scalar::setdata(float data, int i)
{

	(*grid_data)[i] = data;
}

void grid2_scalar::setdata(float data, int i, int j)
{
	// Get 1D Index from 2D Input Index, Deref and float data to grid_data. 
	int index = idx_2Dto1D(i, j);
	(*grid_data)[index] = data; // Dereference Vector (pointer) and then set [index] at [] array offset.
}

void grid2_scalar::adddata(float data, int i)
{
	(*grid_data)[i] += data;
}

void grid2_scalar::adddata(float data, int i, int j)
{
	// Get 1D Index from 2D Input Index, Deref and float data to grid_data. Currently NOT using data * dt. 
	int index = idx_2Dto1D(i, j);
	(*grid_data)[index] += data;
}

//  Get Data With Provided 2D Index (Converted to 1D).
float grid2_scalar::getdata(int i, int j)
{
	return (*grid_data)[idx_2Dto1D(i, j)];
}

// Get Data With Provided 1D Index.
float grid2_scalar::getdata(int i)
{
	return (*grid_data)[i];
}

// Member Function to 0 all grid/vector scalar elements. 
void grid2_scalar::clear()
{
	for (std::size_t i = 0; i < total_size; i++)
	{
		(*grid_data)[i] = 0.0f; 
	}
}

void grid2_scalar::printsize()
{
	std::cout << "grid2_scalar SIZE == " << grid_data->size() << "\n";
}

//#define IX(i,j) ((i)+(N+2) * (j)) N = Size for XY, For this, N should be same/square XY. 
// So use member grid_dim_x as N for now. Possibly will make my code so grid has to be square ..rm comment

// Future Index Function can deal with diffrent ij/XY Sizes (Non Square Grids).

// 2D to 1D Array Index Function - 
int grid2_scalar::idx_2Dto1D(int i, int j)
{
	// Assuming Grid is square for now. (use x size (Has EdgeSize Included)) .
	int NE = x_size + edge_size; // N+Edge_Size (N Is Same for x and y dim (square grid))
	return int((i)+(NE) * (j));
}

// 1D to 2D Array Index Function - 
// Returns 2D Index as vec2 (Cast to float).
vec2<int> grid2_scalar::idx_1Dto2D(int k)
{
	int NE = x_size + edge_size; // N+Edge_Size (N Is Same for x and y dim (square grid))

	int ii = k % (NE); 
	int jj = k / (NE);

	// Cast to float for vec2 return - 
	return vec2<int>(ii, jj);
}

// DEBUG test Dispate Density. (Should be f2obj or solver really). 
void grid2_scalar::disp()
{
	for (int i = 0; i < grid_data->size(); i++)
	{
		(*grid_data)[i] *= 0.9;
	}
}

// SWAP Implemetnation (grid2_scalar)- 

// Even though the swappee object is a sepreate instance, because its same object, will count as
// internal member acess to grid_data member. 

void grid2_scalar::swap(grid2_scalar *B)
{
	// Swap my grid_data vector contents, with Passed B Grid, grid_data vector contents.
	grid_data->swap(*(B->grid_data)); // DeRef B Grid_data vector pointer member.

	// Check for Miss-Allignment of Vector Sizes. All Grids (to be swapped) should be same size. 
	//assert(this->grid_data->size() == B->grid_data->size()); 
	if (grid_data->size() != B->grid_data->size()) { std::cout << "WARN::Grid2_Scalar-Swap::Grid Size Missalignment \n"; }
}

// Getter Function to all external acess to get pointer to Grid_Data Vector. Use with Caution. 
std::vector<float>* grid2_scalar::griddataptr_getter()
{
	return grid_data; 
}

// grid2_scalar print_celldata Implementation, for debugging small grids, to check for Assigned Cell Values. 
void grid2_scalar::print_celldata()
{
	//std::cout << "DEBUG::grid2_scalar " << ID << "Cell Data BEGIN - \n";
	//std::cout << "Xi size = " << grid_dim_x << " Yj size = " << grid_dim_x << " \n \n";

	for (std::size_t j = 1; j <= (std::size_t) y_size; j++)
	{
		for (std::size_t i = 1; i <= (std::size_t) x_size; i++)
		{
			if (i == 1) std::cout << "| ";

			//std::cout << getdata(i, j) << "|  ";
			//std::cout << "* | ";
			//std::cout << "LOOP INDEX::" << "[" << i << "," << j << "]   ";
			//std::cout << "[" << i << "," << j << "]  | ";

			std::cout << getdata(i, j) << " | ";
			
			if (i == (std::size_t) x_size) std::cout << "\n";
		}
	}
	
	//std::cout << "\n" << "DEBUG::grid2_scalar " << ID << " Cell Data END. \n";
}

// TEMP / TESTING MFUNCS  -

// Return Internal Grid_Data std::vector Data Array Ptr.
float* grid2_scalar::getdataarray()
{
	return grid_data->data();
}

// Return Pointer to Newly Allocated grid_data std::vector<float> where grid elements (as flat 1D Array) have been transposed
// This will be costly to lookup grid_data array for column major'd values and set them to cur row major, non cache friendly ! 

// New Implemnetation, pass preallocated grid_data vector by ptr, so not reallocating each kth call. 

void grid2_scalar::tranpose_gridata(std::vector<float>* ptr)
{
	#pragma omp parallel for
	for (int j = 1; j <= x_size; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= y_size; i++)
		{
			(*ptr)[idx_2Dto1D(i, j)] = getdata(j,i); 
		}
	}
	 
}

//----------------------------------------------------------------------\\

// grid2_vector Implementation - \\

// grid2_vector implemented Constructor
grid2_vector::grid2_vector(int x_s, int y_s, int e_s, int idd, float spc) :
	x_size(x_s), y_size(y_s), edge_size(e_s), ID(idd), grid_spacing(spc)
{
	total_size = (x_size + edge_size) * (y_size + edge_size);

	// Alloacte vector<vec2> Grid Data Via Raw Ptr. Init all vec2 elements to (XY 0.0). 
	grid_data = new std::vector<vec2<float>>(total_size, vec2<float>(0.0f, 0.0f));
	//grid_data->resize(total_size); // Resize Vector (Grid_Data) to 1D GridSize x*y.
}

// grid2_vector Copy Constructor - 
// Copy Contents from copy grid to this, on initalization call to Cctor. 
// Deep Copy Grid_Data Memory. (Not just ptr).
grid2_vector::grid2_vector(grid2_vector &copy)
	: x_size(copy.x_size), y_size(copy.y_size), edge_size(copy.edge_size), ID(copy.ID),
	grid_spacing(copy.grid_spacing)
{
	total_size = (x_size + edge_size) * (y_size * edge_size);

	// Copy Grid_Data Directly. 
	grid_data = new std::vector<vec2<float>>((*copy.grid_data));
}

// grid2_scalar Copy Assignment Operator Overload - 
// Copy Contents from copy grid to this, post initalization. 
// Deep Copy Grid_Data Memory. (Not just ptr).
grid2_vector&::grid2_vector::operator=(grid2_vector &copy)
{
	// Check if this passed. 
	if (&copy == this) return *this;

	// Dealloc current grid_data vector. 
	delete this->grid_data;

	// Copy grid_data vector from copy grid. 
	grid_data = new std::vector<vec2<float>>((*copy.grid_data));

	// Copy Local Member Values 
	this->x_size = copy.x_size; this->y_size = copy.y_size; this->edge_size = copy.edge_size;
	this->total_size = (x_size + edge_size) * (y_size + edge_size);
	this->ID = copy.ID + 5; this->grid_spacing = copy.grid_spacing;

	return *this;
}


// grid2_vector destructor. 
grid2_vector::~grid2_vector()
{
	// Clear grid_data member (Vector Memory) - 
	delete grid_data;
	grid_data = nullptr;
}

// Util Member Function to Print Grid Info 
void grid2_vector::printinfo()
{
	std::cout << "DEBUG::Grid 2D Vector " << ID << " Info BEGIN - \n";
	std::cout << "Grid Cell Count = " << total_size << "\n";
	std::cout << "Grid X Row Size = " << x_size << "\n";
	std::cout << "Grid Y Row Size = " << y_size << "\n";
	std::cout << "Grid Edge Cell Size = " << edge_size << "\n";
	std::cout << "DEBUG::Grid 2D Vector " << ID << " Info END. \n \n";
}

// grid2_vector Acess Member Functions - 
// SET DATA - 

void grid2_vector::setdata(vec2<float> data, int i)
{
	(*grid_data)[i] = data;
}

// Member Function to set Grid Data Of Vector Grid. 
void grid2_vector::setdata(vec2<float> data, int i, int j)
{
	// Get 1D Index from 2D Input Index, Deref and vec2data to grid_data. 
	int index = idx_2Dto1D(i, j);
	(*grid_data)[index] = data;
}

void grid2_vector::setdata_x(float xx, int i)
{
	(*grid_data)[i].x = xx;
}

// Member Function to set Grid Data vec2 X Compoennt Only Of Vector Grid. 
void grid2_vector::setdata_x(float xx, int i, int j)
{
	// Get 1D Index from 2D Input Index.
	int index = idx_2Dto1D(i, j);
	// Set X Component of De-Referenced vec2 in grid_data stdvector. 
	(*grid_data)[index].x = xx; 
}

void grid2_vector::setdata_y(float yy, int i)
{
	(*grid_data)[i].y = yy;
}

// Member Function to set Grid Data vec2 X Compoennt Only Of Vector Grid. 
void grid2_vector::setdata_y(float yy, int i, int j)
{
	// Get 1D Index from 2D Input Index.
	int index = idx_2Dto1D(i, j);
	// Set Y Component of De-Referenced vec2 in grid_data stdvector. 
	(*grid_data)[index].y = yy;
}

void grid2_vector::adddata(vec2<float> data, int i)
{
	vec2<float> cur = (*grid_data)[i];

	// Add Cur to passed data vec2 - 
	(*grid_data)[i] = cur + data;
}

// Member Function to Add Grid Data Of Vector Grid. 
void grid2_vector::adddata(vec2<float> data, int i, int j)
{
	// Get 1D Index from 2D Input Index, Deref and vec2data to grid_data. 
	int index = idx_2Dto1D(i, j);
	vec2<float> cur = (*grid_data)[index]; // Bind Cur Vec2 to local vec2

	// Add Cur to passed data vec2 - 
	(*grid_data)[index] = cur + data;
}

// GET DATA - 

vec2<float> grid2_vector::getdata(int i)
{
	return (*grid_data)[i];
}

// Member Function to get vec2 grid_data at 2D Index ij. 
vec2<float> grid2_vector::getdata(int i, int j)
{
	return (*grid_data)[idx_2Dto1D(i, j)]; 
}

float grid2_vector::getdata_x(int i)
{
	return (*grid_data)[i].x;
}

// Member Function to get vec2 X Component Float grid_data at 2D Index ij. 
float grid2_vector::getdata_x(int i, int j)
{
	return (*grid_data)[idx_2Dto1D(i, j)].x; 
}

float grid2_vector::getdata_y(int i)
{
	return (*grid_data)[i].y;
}

// Member Function to get vec2 Y Component Float grid_data at 2D Index ij. 
float grid2_vector::getdata_y(int i, int j)
{
	return (*grid_data)[idx_2Dto1D(i, j)].y;
}

// Member Function to 0 all grid/vector scalar elements. 
void grid2_vector::clear()
{
	for (std::size_t i = 0; i < total_size; i++)
	{
		(*grid_data)[i] = vec2<float>(0.0f, 0.0f);
	}
}

// 2D to 1D Array Index Function - (grid2_vector) 
int grid2_vector::idx_2Dto1D(int i, int j)
{
	// Assuming Grid is square for now. (use x size (Has EdgeSize Included)) .
	int NE = x_size + edge_size; // N+Edge_Size (N Is Same for x and y dim (square grid))
	return int((i)+(NE) * (j));
}

// 1D to 2D Array Index Function - (grid2_vector)
// Returns 2D Index as vec2 (Cast to float). (Uses same vec2 class as grid2_vector grid data). 
vec2<int> grid2_vector::idx_1Dto2D(int k)
{
	int NE = x_size + edge_size; // N+Edge_Size (N Is Same for x and y dim (square grid))

	int ii = k % (NE);
	int jj = k / (NE);

	// Cast to float for vec2 return - 
	return vec2<int>(ii, jj);
}

// SWAP Implemetnation (grid2_vector)- 

// Even though the swappee object is a sepreate instance, because its same object, will count as
// internal member acess to grid_data member.

void grid2_vector::swap(grid2_vector *B)
{
	// Swap my grid_data vector contents, with Passed B Grid, grid_data vector contents.
	grid_data->swap(*(B->grid_data)); // DeRef B Grid_data vector pointer member.

	// Check for Miss-Allignment of Vector Sizes. All Grids (to be swapped) should be same size. 
	//assert(this->grid_data->size() == B->grid_data->size()); 
	if (grid_data->size() != B->grid_data->size()) { std::cout << "WARN::Grid2_Vector-Swap::Grid Size Missalignment \n"; }
}

// grid2_vector print_celldata Implementation, for debugging small grids, to check for Assigned Cell Values. 
void grid2_vector::print_celldata()
{
	//std::cout << "DEBUG::grid2_scalar " << ID << "Cell Data BEGIN - \n";
	//std::cout << "Xi size = " << grid_dim_x << " Yj size = " << grid_dim_x << " \n \n";

	for (std::size_t j = 1; j <= (std::size_t) y_size; j++)
	{
		for (std::size_t i = 1; i <= (std::size_t) x_size; i++)
		{
			if (i == 1) std::cout << "| ";

			//std::cout << getdata(i, j) << "|  ";
			//std::cout << "* | ";
			//std::cout << "LOOP INDEX::" << "[" << i << "," << j << "]   ";
			//std::cout << "[" << i << "," << j << "]  | ";

			std::cout << getdata(i, j).x << "," << getdata(i, j).y << " | ";

			if (i == (std::size_t) x_size) std::cout << "\n";
		}
	}

	//std::cout << "\n" << "DEBUG::grid2_scalar " << ID << " Cell Data END. \n";
}

// Grid_Data Ptr and Grid_Data->data() (internal std::vector data array) Ptr Getters - 

// Grid_Data std::vector<vec2> Ptr Getter - 
std::vector<vec2<float>>* grid2_vector::griddataptr_getter()
{
	return grid_data; 
}

// Grid Data std::vector<vec2>::data Internal Data Array Ptr Getter 
vec2<float>* grid2_vector::getdataarray()
{
	return grid_data->data(); 
}
*/

//----------------------------------------------------------------------\\

// Grid3 ABC Implmentation \\ 

template <class T>
grid3<T>::grid3<T>(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s)
	: x_size(x_s), y_size(y_s), z_size(z_s), edge_size(e_s)
{
	total_size = (x_size + edge_size) * (y_size + edge_size) * (z_size + edge_size);

	grid_data = new std::vector<T>(total_size, T());
}

// ABC Resonsible for Alloc and Dealloc of Grid_Data. 
template <class T>
grid3<T>::~grid3() 
{
	if (grid_data || grid_data != nullptr)
	{
		delete grid_data; grid_data = nullptr; 
	}
}

// Grid 3 ABC Virtual Default Implmenetation - 

// Grid_Data - SETTERS \\

template <class T>
void grid3<T>::setdata(T data, int i)
{
	(*grid_data)[i] = data;
}

template <class T>
void grid3<T>::setdata(T data, int i, int j, int k)
{
	int index = idx_3Dto1D(i, j, k);
	(*grid_data)[index] = data;
}

template <class T>
void grid3<T>::adddata(T data, int i)
{
	(*grid_data)[i] += data;
}

template <class T>
void grid3<T>::adddata(T data, int i, int j, int k)
{
	int index = idx_3Dto1D(i, j, k);
	(*grid_data)[index] += data;
}

template <class T>
T grid3<T>::getdata(int i, int j, int k) const
{
	return (*grid_data)[idx_3Dto1D(i, j, k)];
}
template <class T>
T grid3<T>::getdata(int i) const
{
	return (*grid_data)[i];
}

// NOTE: Could just take std::vector<T> Ptr or Reference to avoid doing the grid3 abc downcast to derived grid_[..] class. 
template <class T>
void grid3<T>::swap(const grid3<T> *B)
{
	// Check Grid_Data Sizes Match -  
	assert(this->grid_data->size() == B->grid_data->size());
	// Check Grid_Data Element Type Sizes Match -
	assert(sizeof(this->grid_data->at(0)) == sizeof(B->grid_data->at(0)));
	// Swap grid_data vectors. 
	this->grid_data->swap(*(B->grid_data));
}

template <class T>
std::size_t grid3<T>::get_edgesize() const
{
	return edge_size; 
}

template <class T>
vec3<std::size_t> grid3<T>::get_dimsize() const
{
	return vec3<std::size_t>(x_size, y_size, z_size);
}

// Grid Data Pointer / Data Array Pointer Getters - 

// Returns grid_data (std::vector<T>) ptr, for external acess 
// ! Use with caution
template <class T>
std::vector<T>* grid3<T>::griddataptr_getter() const
{
	return grid_data;
}

// Returns grid_data std::vector<T>::data() internal data array pointer, for external acess.
// ! Use with caution
template <class T>
T* grid3<T>::getdataarray() const
{
	return grid_data->data();
}

// Indexers - 

template <class T>
int grid3<T>::idx_3Dto1D(int i, int j, int k) const
{
	return (int) i + (x_size + edge_size) * (j + (z_size + edge_size) * k);
}

template <class T>
vec3<int> grid3<T>::idx_1Dto3D(int i) const
{
	int ii = i / ((y_size + edge_size) * (z_size + edge_size));
	int jj = (i / (z_size + edge_size)) % (y_size + edge_size);
	int kk = i % (z_size + edge_size);

	return vec3<int>(ii, jj, kk);
}


//----------------------------------------------------------------------\\

// Grid3_Scalar Implmentation \\

// grid3_scalar shares most implmentation with base grid3 class, apart from implmeneting pure virtual MFs. 

// Note: explicit use's of this-> due to C3861 from Dependemt Base Class Template <T> Members. (No Two Phase Name Lookup)

// Inilzation of Base grid3 Class Constructor. No Grid3_Scalar Specfic Members to initalize.
template <class T>
grid3_scalar<T>::grid3_scalar<T>(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s)
	: grid3<T>(x_s, y_s, z_s, e_s) {}

template <class T>
grid3_scalar<T>::~grid3_scalar()
{
	// Grid Data Dealloc by ABC grid3 Destructor. 
}

// Pure Virtual MFs Implementation - 

template <class T>
void grid3_scalar<T>::clear() 
{
	for (std::size_t i = 0; i < this->total_size; i++)
	{
		(*(this->grid_data))[i] = (T) 0;
	}
}

template <class T>
void grid3_scalar<T>::printinfo() const
{
	std::cout << "DEBUG::Grid 2D Scalar " << " Info BEGIN - \n";
	std::cout << "Grid Cell Count = " << this->total_size << "\n";
	std::cout << "Grid X Row Size = " << this->x_size << "\n";
	std::cout << "Grid Y Row Size = " << this->y_size << "\n";
	std::cout << "Grid Edge Cell Size = " << this->edge_size << "\n";
	std::cout << "DEBUG::Grid 2D Scalar " << " Info END. \n \n";
}

//----------------------------------------------------------------------\\

// Grid3_Vector Implmentation \\

// grid3_vector shares some implementations with Base class, adds per component getter/setter MFs, aswell as defining Pure Virtaul MFs. 

// Note: explicit use's of this-> due to C3861 from Dependemt Base Class Template <T> Members. (No Two Phase Name Lookup)

// Inilzation of Base grid3 Class Constructor. No Grid3_Vector Specfic Members to initalize.
template <class T>
grid3_vector<T>::grid3_vector<T>(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s)
	: grid3<T>(x_s, y_s, z_s, e_s) {}

template <class T>
grid3_vector<T>::~grid3_vector()
{
	// Grid Data Dealloc by ABC grid3 Destructor. 
}

// Pure Virtual MFs Implementation - 

template <class T>
void grid3_vector<T>::clear()
{
	for (std::size_t i = 0; i < this->total_size; i++)
	{
		(*(this->grid_data))[i] = (T)0;
	}
}

template <class T>
void grid3_vector<T>::printinfo() const
{
	std::cout << "DEBUG::Grid 2D Vector " << " Info BEGIN - \n";
	std::cout << "Grid Cell Count = " << this->total_size << "\n";
	std::cout << "Grid X Row Size = " << this->x_size << "\n";
	std::cout << "Grid Y Row Size = " << this->y_size << "\n";
	std::cout << "Grid Edge Cell Size = " << this->edge_size << "\n";
	std::cout << "DEBUG::Grid 2D Scalar " << " Info END. \n \n";
}

// grid3_vector Specfic, Component Wise Setters/Getters - 

// Set X
template <class T>
void grid3_vector<T>::setdata_x(T xx, int i)
{
	(*(this->grid_data))[i].x = xx;
}
template <class T>
void grid3_vector<T>::setdata_x(T xx, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].x = xx;
}
// Set Y
template <class T>
void grid3_vector<T>::setdata_y(T yy, int i)
{
	(*(this->grid_data))[i].y = yy;
}
template <class T>
void grid3_vector<T>::setdata_y(T yy, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].y = yy;
}
// Set Z
template <class T>
void grid3_vector<T>::setdata_z(T zz, int i)
{
	(*(this->grid_data))[i].z = zz;
}
template <class T>
void grid3_vector<T>::setdata_z(T zz, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].z = zz;
}

// Add X
template <class T>
void grid3_vector<T>::adddata_x(T data, int i)
{
	(*(this->grid_data))[i].x += data;
}
template <class T>
void grid3_vector<T>::adddata_x(T data, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i,j,k)].x += data;
}
// Add Y
template <class T>
void grid3_vector<T>::adddata_y(T data, int i)
{
	(*(this->grid_data))[i].y += data;
}
template <class T>
void grid3_vector<T>::adddata_y(T data, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].y += data;
}
// Add Z
template <class T>
void grid3_vector<T>::adddata_z(T data, int i)
{
	(*(this->grid_data))[i].z += data;
}
template <class T>
void grid3_vector<T>::adddata_z(T data, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].z += data;
}

// Get X
template <class T>
T grid3_vector<T>::getdata_x(int i) const
{
	return (*(this->grid_data))[i].x;
}
template <class T>
T grid3_vector<T>::getdata_x(int i, int j, int k) const
{
	return (*(this->grid_data))[this->idx_3Dto1D(i, j, k)].x;
}
// Get Y
template <class T>
T grid3_vector<T>::getdata_y(int i) const
{
	return (*(this->grid_data))[i].y;
}
template <class T>
T grid3_vector<T>::getdata_y(int i, int j, int k) const
{
	return (*(this->grid_data))[this->idx_3Dto1D(i, j, k)].y;
}
// Get X
template <class T>
T grid3_vector<T>::getdata_z(int i) const
{
	return (*(this->grid_data))[i].z;
}
template <class T>
T grid3_vector<T>::getdata_z(int i, int j, int k) const
{
	return (*(this->grid_data))[this->idx_3Dto1D(i, j, k)].z;
}


//----------------------------------------------------------------------\\

// Explicit Instations - (Defintion of Templated Classes is Sepreated into this source file).

// Grid3 ABC
template class grid3<float>;
template class grid3<double>;
template class grid3<vec3<float>>;
template class grid3<vec3<double>>;

// Grid3 Scalar
template class grid3_scalar<float>;
template class grid3_scalar<double>;

/*
// Grid3 Vector
template class grid3_vector<vec3<float>>;
template class grid3_vector<vec3<double>>;
*/

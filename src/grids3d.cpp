// Implementation of grids3d
#include "grids3d.h"

// Project Headers
#include "vec3d.h"

// Std Headers
#include <memory>
#include <fstream>
#include <string>
#include <cassert>

#include <omp.h> 


extern short verbose; // Get Verbose Global from main. 

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

template <class T>
void grid3<T>::swap(const grid3<T> *B)
{
	assert(this->grid_data->size() == B->grid_data->size());
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

// Returns grid_data (std::vector<T>) ptr, for external acess !Use with caution
template <class T>
std::vector<T>* grid3<T>::griddataptr_getter() const
{
	return grid_data;
}

// Returns grid_data std::vector<T>::data() internal data array pointer, for external acess. !Use with caution
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

/* Grid3_Scalar Implmentation 
	grid3_scalar shares most implmentation with base grid3 class, apart from implmeneting pure virtual MFs. 
	Note: explicit use's of this-> due to C3861 from Dependemt Base Class Template <T> Members. (No Two Phase Name Lookup)
*/

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
		(*(this->grid_data))[i] = T {};
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

/*	Grid3_Vector Implmentation 
	grid3_vector shares some implementations with Base class, adds per component getter/setter MFs, aswell as defining/ovr Pure Virtaul MFs.
	Note: explicit use's of this-> due to C3861 from Dependemt Base Class Template <T> Members. (No Two Phase Name Lookup)
*/

// Inilzation of Base grid3 Class Constructor. No Grid3_Vector Specfic Members to initalize.
template <class T>
grid3_vector<T>::grid3_vector<T>(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s)
	: grid3<T>(x_s, y_s, z_s, e_s) {}

template <class T>
grid3_vector<T>::~grid3_vector()
{
	// Grid Data Dealloc by ABC grid3 dtor.
}

// Pure Virtual MFs Implementation - 
template <class T>
void grid3_vector<T>::clear()
{
	for (std::size_t i = 0; i < this->total_size; i++)
	{
		(*(this->grid_data))[i] = T {};
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

// Vec3 Component Setters -
// Set X
template <class T>
void grid3_vector<T>::setdata_x(float xx, int i)
{
	(*(this->grid_data))[i].x = xx;
}
template <class T>
void grid3_vector<T>::setdata_x(float xx, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].x = xx;
}
// Set Y
template <class T>
void grid3_vector<T>::setdata_y(float yy, int i)
{
	(*(this->grid_data))[i].y = yy;
}
template <class T>
void grid3_vector<T>::setdata_y(float yy, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].y = yy;
}
// Set Z
template <class T>
void grid3_vector<T>::setdata_z(float zz, int i)
{
	(*(this->grid_data))[i].z = zz;
}
template <class T>
void grid3_vector<T>::setdata_z(float zz, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].z = zz;
}

// Vec3 Component Adders-
// Add X
template <class T>
void grid3_vector<T>::adddata_x(float data, int i)
{
	(*(this->grid_data))[i].x += data;
}
template <class T>
void grid3_vector<T>::adddata_x(float data, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i,j,k)].x += data;
}
// Add Y
template <class T>
void grid3_vector<T>::adddata_y(float data, int i)
{
	(*(this->grid_data))[i].y += data;
}
template <class T>
void grid3_vector<T>::adddata_y(float data, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].y += data;
}
// Add Z
template <class T>
void grid3_vector<T>::adddata_z(float data, int i)
{
	(*(this->grid_data))[i].z += data;
}
template <class T>
void grid3_vector<T>::adddata_z(float data, int i, int j, int k)
{
	(*(this->grid_data))[this->idx_3Dto1D(i, j, k)].z += data;
}

// Vec3 Component Getters -
// Get X
template <class T>
float grid3_vector<T>::getdata_x(int i) const
{
	return (*(this->grid_data))[i].x;
}
template <class T>
float grid3_vector<T>::getdata_x(int i, int j, int k) const
{
	return (*(this->grid_data))[this->idx_3Dto1D(i, j, k)].x;
}
// Get Y
template <class T>
float grid3_vector<T>::getdata_y(int i) const
{
	return (*(this->grid_data))[i].y;
}
template <class T>
float grid3_vector<T>::getdata_y(int i, int j, int k) const
{
	return (*(this->grid_data))[this->idx_3Dto1D(i, j, k)].y;
}
// Get X
template <class T>
float grid3_vector<T>::getdata_z(int i) const
{
	return (*(this->grid_data))[i].z;
}
template <class T>
float grid3_vector<T>::getdata_z(int i, int j, int k) const
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

// Grid3 Vector
// 3D Grid - 3D Vectors
template class grid3_vector<vec3<float>>;
template class grid3_vector<vec3<double>>;


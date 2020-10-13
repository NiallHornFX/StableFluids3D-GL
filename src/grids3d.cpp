// Implementation of grids3d
#include "grids3d.h"

// Project Headers
#include "vec3d.h"

// Std Headers
#include <algorithm>
#include <memory>
#include <fstream>
#include <string>
#include <cassert>

// Ext
#include <omp.h> 

extern short verbose; // Get Verbose Global from main. 

//----------------------------------------------------------------------

// Grid3 ABC Implmentation

template <class T>
grid3<T>::grid3(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s)
	: x_size(x_s), y_size(y_s), z_size(z_s), edge_size(e_s)
{
	total_size = (x_size + edge_size) * (y_size + edge_size) * (z_size + edge_size);

	grid_data = std::vector<T>(total_size, T());
}

// ABC Resonsible for Alloc and Dealloc of Grid_Data. (Grid_Data now handled by std::vector)
template <class T>
grid3<T>::~grid3() 
{

}

// Grid 3 ABC Virtual Default Implmenetation - 

// Grid_Data - SETTERS

template <class T>
void grid3<T>::setdata(T data, int i)
{
	grid_data[i] = data;
}

template <class T>
void grid3<T>::setdata(T data, int i, int j, int k)
{
	int index = idx_3Dto1D(i, j, k);
	grid_data[index] = data;
}

template <class T>
void grid3<T>::adddata(T data, int i)
{
	grid_data[i] += data;
}

template <class T>
void grid3<T>::adddata(T data, int i, int j, int k)
{
	int index = idx_3Dto1D(i, j, k);
	grid_data[index] += data;
}

template <class T>
T grid3<T>::getdata(int i, int j, int k) const
{
	return grid_data[idx_3Dto1D(i, j, k)];
}
template <class T>
T grid3<T>::getdata(int i) const
{
	return grid_data[i];
}

template <class T>
void grid3<T>::swap(const grid3<T> *B)
{
	assert(this->grid_data.size() == B->grid_data.size());
//	assert(sizeof(this.grid_data.at(0)) == sizeof(B->grid_data.at(0)));

	// Swap grid_data vectors. 
	this->grid_data.swap( const_cast<std::vector<T>&>(B->grid_data));
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
// NOTE Returns Direct Reference to this grids grid_data vector. !NON Const so it can be re-assigned in setcurtoprev. Use With Caution. 
template <class T> 
std::vector<T>& grid3<T>::griddatavector_getter() 
{
	return grid_data;
}


// Returns grid_data std::vector<T>::data() internal data array pointer, for external acess. !Use with caution
template <class T>
T* grid3<T>::getdataarray() const
{
	T *ret = (const_cast<std::vector<T>&>(grid_data)).data();
	return ret;
}

// Indexers (NonStatic Inline) & Grid-IDX-Grid Space Converters (Static Inline, implemented in Header). 


//----------------------------------------------------------------------

/* Grid3_Scalar Implmentation 
	grid3_scalar shares most implmentation with base grid3 class, apart from implmeneting pure virtual MFs. 
	Note: explicit use's of this-> due to C3861 from Dependent Base Class Template <T> Members. (No Two Phase Name Lookup)
*/

// Inilzation of Base grid3 Class Constructor. No Grid3_Scalar Specfic Members to initalize.
template <class T>
grid3_scalar<T>::grid3_scalar(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s)
	: grid3<T>(x_s, y_s, z_s, e_s) {}

template <class T>
grid3_scalar<T>::~grid3_scalar()
{
	// Grid Data Dealloc by ABC grid3 Destructor. 
}

// Pure Virtual MFs Implementation - 

// grid3_scalar Sampler Implementation - 
// Sample Grid Space Location by calculating Indices and Coefficents for select Interoplation method. 
template <class T>
T grid3_scalar<T>::sampler(const vec3<float> &gs_loc, vec3<int> &r_indices, interpType interp) const
{
	// 1D Scalar Interp GLambdas.
	auto lerp = []( auto val_0, auto val_1, float bias) -> auto // decltype(grid_data.at(0)) 
	{
		return (1.0f - bias) * val_0 + bias * val_1;
	};

	auto cosinterp = [](auto val_0, auto val_1, auto bias) -> auto
	{
		auto mu = (1.0f - std::cos(bias*PI)) / 2.0f;
		return  (1.0f - mu) * val_0 + val_1 * mu; // (decltype(this->grid_data.at(0)))
	};

	// Deduce Grid Indices and Coefficents for Interoplation From Traced Location - 
	vec3<float> idx_loc = this->idx_gridToIndex(gs_loc.x, gs_loc.y, gs_loc.z, this->x_size); // Grid-->IdxSpace
	// Cell (i,j,k|+1) Indices 
	int i0 = int(idx_loc.x); int i1 = i0 + 1;
	int j0 = int(idx_loc.y); int j1 = j0 + 1;
	int k0 = int(idx_loc.z); int k1 = k0 + 1;
	r_indices = vec3<int>(i0, j0, k0); // Write Indices to vec3<int> ref
	// Interoplation Coefficents - 
	float r = idx_loc.x - i0;
	float s = idx_loc.y - j0;
	float t = idx_loc.z - k0;

	// Sample Grid Using Interoplation Method - 
	if (interp == interp_TriLinear)
	{
		// Scalar TriLinear Interoplation <T>
		T L_000_001_t = lerp(this->getdata(i0, j0, k0), this->getdata(i0, j0, k1), t); // Z
		T L_010_011_t = lerp(this->getdata(i0, j1, k0), this->getdata(i0, j1, k1), t);
		T L_100_101_t = lerp(this->getdata(i1, j0, k0), this->getdata(i1, j0, k1), t);
		T L_110_111_t = lerp(this->getdata(i1, j1, k0), this->getdata(i1, j1, k1), t);
		T L_A = lerp(L_000_001_t, L_010_011_t, s); // Y
		T L_B = lerp(L_100_101_t, L_110_111_t, s);
		T L_F = lerp(L_A, L_B, r); // X

		return L_F; 
	}
	else if (interp == interp_TriCosine)
	{
		// Scalar TriCosine Interoplation <T>
		T C_000_001_t = cosinterp(this->getdata(i0, j0, k0), this->getdata(i0, j0, k1), t); // Z
		T C_010_011_t = cosinterp(this->getdata(i0, j1, k0), this->getdata(i0, j1, k1), t);
		T C_100_101_t = cosinterp(this->getdata(i1, j0, k0), this->getdata(i1, j0, k1), t);
		T C_110_111_t = cosinterp(this->getdata(i1, j1, k0), this->getdata(i1, j1, k1), t);
		T C_A = cosinterp(C_000_001_t, C_010_011_t, s); // Y
		T C_B = cosinterp(C_100_101_t, C_110_111_t, s);
		T C_F = cosinterp(C_A, C_B, r); // X

		return C_F;
	}
	else
	{
		std::cerr << "ERR::Invalid Interoplation Type Specifed \n \n";
		std::terminate();
	}
}

// grid3 MinMax Scalar Implementation -
// Deduce Scalar MinAndMax of <T> of cell neighbours at passed Index Space Grid Location. 
template <class T>
std::tuple<T, T> grid3_scalar<T>::minmax(int i0, int j0, int k0) const
{
	int i1 = i0 + 1, j1 = j0 + 1, k1 = k0 + 1; 

	// Cell Indices ijk(0|+1) to get Min/Max of (For MacCormack Limiter of Quanitiy (Scalar))
	T C_000 = this->getdata(i0, j0, k0); T C_001 = this->getdata(i0, j0, k1);
	T C_010 = this->getdata(i0, j1, k1); T C_011 = this->getdata(i0, j1, k0);
	T C_100 = this->getdata(i1, j0, k0); T C_101 = this->getdata(i1, j0, k1);
	T C_110 = this->getdata(i1, j1, k1); T C_111 = this->getdata(i1, j1, k0);

	T f_min = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000, C_001), C_010), C_011), C_100), C_101), C_110), C_111);
	T f_max = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000, C_001), C_010), C_011), C_100), C_101), C_110), C_111);

	return std::make_tuple(f_min, f_max);
}

template <class T>
void grid3_scalar<T>::clear() 
{
	for (std::size_t i = 0; i < this->total_size; i++)
	{
		(this->grid_data)[i] = T {};
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

//----------------------------------------------------------------------

/*	Grid3_Vector Implmentation 
	grid3_vector shares some implementations with Base class, adds per component getter/setter MFs, aswell as defining/ovr Pure Virtaul MFs.
	Note: explicit use's of this-> due to C3861 from Dependemt Base Class Template <T> Members. (No Two Phase Name Lookup)
*/

// Inilzation of Base grid3 Class Constructor. No Grid3_Vector Specfic Members to initalize.
template <class T>
grid3_vector<T>::grid3_vector(std::size_t x_s, std::size_t y_s, std::size_t z_s, std::size_t e_s)
	: grid3<T>(x_s, y_s, z_s, e_s) {}

template <class T>
grid3_vector<T>::~grid3_vector()
{
	// Grid Data Dealloc by ABC grid3 dtor.
}

// Pure Virtual MFs Implementation - 

// grid3_vector Sampler Implementation - 
// Sample Grid Space Location by calculating Indices and Coefficents for select Interoplation method. Hardcoded vec3<float> internally for now. 
template <class T>
T grid3_vector<T>::sampler(const vec3<float> &gs_loc, vec3<int> &r_indices, interpType interp) const
{
	auto vec_lerp = [](const vec3<float> &v_a, const vec3<float> &v_b, float bias) -> vec3<float>
	{
		float xx = (1.0f - bias) * v_a.x + bias * v_b.x;
		float yy = (1.0f - bias) * v_a.y + bias * v_b.y;
		float zz = (1.0f - bias) * v_a.z + bias * v_b.z;
		return vec3<float>(xx, yy, zz);
	};
	auto vec_cerp = [](const vec3<float> &v_a, const vec3<float> &v_b, float bias) -> vec3<float>
	{
		float mu, xx, yy, zz; 
		mu = (1.0f - std::cos(bias*PI)) / 2.0f;
		xx = (1.0f - mu) * v_a.x + v_b.x * mu;
		yy = (1.0f - mu) * v_a.y + v_b.y * mu;
		zz = (1.0f - mu) * v_a.z + v_b.z * mu;
		return vec3<float>(xx, yy, zz);
	};

	// Deduce Grid Indices and Coefficents for Interoplation From Traced Location - 
	vec3<float> idx_loc = this->idx_gridToIndex(gs_loc.x, gs_loc.y, gs_loc.z, this->x_size); // Grid-->IdxSpace
	// Cell (i,j,k|+1) Indices 
	int i0 = int(idx_loc.x); int i1 = i0 + 1;
	int j0 = int(idx_loc.y); int j1 = j0 + 1;
	int k0 = int(idx_loc.z); int k1 = k0 + 1;
	r_indices = vec3<int>(i0, j0, k0); // Write Indices to vec3<int> ref
	// Interoplation Coefficents - 
	float r = idx_loc.x - i0;
	float s = idx_loc.y - j0;
	float t = idx_loc.z - k0;

	// Sample Grid Using Interoplation Method - 
	if (interp == interp_TriLinear)
	{
		// Vector TriLinear Interoplation
		vec3<float> L_000_001_t = vec_lerp(this->getdata(i0, j0, k0), this->getdata(i0, j0, k1), t); // Z
		vec3<float> L_010_011_t = vec_lerp(this->getdata(i0, j1, k0), this->getdata(i0, j1, k1), t);
		vec3<float> L_100_101_t = vec_lerp(this->getdata(i1, j0, k0), this->getdata(i1, j0, k1), t);
		vec3<float> L_110_111_t = vec_lerp(this->getdata(i1, j1, k0), this->getdata(i1, j1, k1), t);
		vec3<float> L_A = vec_lerp(L_000_001_t, L_010_011_t, s); // Y
		vec3<float> L_B = vec_lerp(L_100_101_t, L_110_111_t, s);
		vec3<float> L_F = vec_lerp(L_A, L_B, r); // X

		return L_F;
	}
	else if (interp == interp_TriCosine)
	{
		// Vector TriCosine Interoplation 
		vec3<float> C_000_001_t = vec_cerp(this->getdata(i0, j0, k0), this->getdata(i0, j0, k1), t); // Z
		vec3<float> C_010_011_t = vec_cerp(this->getdata(i0, j1, k0), this->getdata(i0, j1, k1), t);
		vec3<float> C_100_101_t = vec_cerp(this->getdata(i1, j0, k0), this->getdata(i1, j0, k1), t);
		vec3<float> C_110_111_t = vec_cerp(this->getdata(i1, j1, k0), this->getdata(i1, j1, k1), t);
		vec3<float> C_A = vec_cerp(C_000_001_t, C_010_011_t, s); // Y
		vec3<float> C_B = vec_cerp(C_100_101_t, C_110_111_t, s);
		vec3<float> C_F = vec_cerp(C_A, C_B, r); // X

		return C_F;
	}
	else
	{
		std::cerr << "ERR::Invalid Interoplation Type Specifed \n \n";
		std::terminate();
	}

}

// grid3 MinMax Vector Implementation -
/* Deduce Vector MinAndMax of <T> of cell neighbours at passed Index Space Grid Location. 
   Due to T<?> vec3<float> hardcoded for now as internal vec3 type. */

template <class T>
std::tuple<T, T> grid3_vector<T>::minmax(int i0, int j0, int k0) const
{
	int i1 = i0 + 1, j1 = j0 + 1, k1 = k0 + 1;

	// Cell Indices ijk(0|+1) to get Min/Max of (For MacCormack Limiter of Quanitiy (Scalar))
	float C_000_x = this->getdata_x(i0, j0, k0); float C_001_x = this->getdata_x(i0, j0, k1);
	float C_010_x = this->getdata_x(i0, j1, k1); float C_011_x = this->getdata_x(i0, j1, k0);
	float C_100_x = this->getdata_x(i1, j0, k0); float C_101_x = this->getdata_x(i1, j0, k1);
	float C_110_x = this->getdata_x(i1, j1, k1); float C_111_x = this->getdata_x(i1, j1, k0);
	float C_000_y = this->getdata_y(i0, j0, k0); float C_001_y = this->getdata_y(i0, j0, k1);
	float C_010_y = this->getdata_y(i0, j1, k1); float C_011_y = this->getdata_y(i0, j1, k0);
	float C_100_y = this->getdata_y(i1, j0, k0); float C_101_y = this->getdata_y(i1, j0, k1);
	float C_110_y = this->getdata_y(i1, j1, k1); float C_111_y = this->getdata_y(i1, j1, k0);
	float C_000_z = this->getdata_z(i0, j0, k0); float C_001_z = this->getdata_z(i0, j0, k1);
	float C_010_z = this->getdata_z(i0, j1, k1); float C_011_z = this->getdata_z(i0, j1, k0);
	float C_100_z = this->getdata_z(i1, j0, k0); float C_101_z = this->getdata_z(i1, j0, k1);
	float C_110_z = this->getdata_z(i1, j1, k1); float C_111_z = this->getdata_z(i1, j1, k0);

	// Optimize this ! 
	float f_min_x = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_x, C_001_x), C_010_x), C_011_x), C_100_x), C_101_x), C_110_x), C_111_x);
	float f_max_x = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_x, C_001_x), C_010_x), C_011_x), C_100_x), C_101_x), C_110_x), C_111_x);
	float f_min_y = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_y, C_001_y), C_010_y), C_011_y), C_100_y), C_101_y), C_110_y), C_111_y);
	float f_max_y = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_y, C_001_y), C_010_y), C_011_y), C_100_y), C_101_y), C_110_y), C_111_y);
	float f_min_z = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_z, C_001_z), C_010_z), C_011_z), C_100_z), C_101_z), C_110_z), C_111_z);
	float f_max_z = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_z, C_001_z), C_010_z), C_011_z), C_100_z), C_101_z), C_110_z), C_111_z);

	vec3<float> v_min(f_min_x, f_min_y, f_min_z);
	vec3<float> v_max(f_max_x, f_max_y, f_max_z);

	// Implict cast of std::tuple<vec3<float>, vec3<float>> to std::tuple<T,T> this is fine for now as vec3<double> is only other grid3_vec type. (If grid3_vec <T> == vec3<double>)
	return std::make_tuple(v_min, v_max); 
}

template <class T>
void grid3_vector<T>::clear()
{
	for (std::size_t i = 0; i < this->total_size; i++)
	{
		this->grid_data[i] = T {};
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
	(this->grid_data)[i].x = xx;
}
template <class T>
void grid3_vector<T>::setdata_x(float xx, int i, int j, int k)
{
	(this->grid_data)[this->idx_3Dto1D(i, j, k)].x = xx;
}
// Set Y
template <class T>
void grid3_vector<T>::setdata_y(float yy, int i)
{
	(this->grid_data)[i].y = yy;
}
template <class T>
void grid3_vector<T>::setdata_y(float yy, int i, int j, int k)
{
	(this->grid_data)[this->idx_3Dto1D(i, j, k)].y = yy;
}
// Set Z
template <class T>
void grid3_vector<T>::setdata_z(float zz, int i)
{
	(this->grid_data)[i].z = zz;
}
template <class T>
void grid3_vector<T>::setdata_z(float zz, int i, int j, int k)
{
	(this->grid_data)[this->idx_3Dto1D(i, j, k)].z = zz;
}

// Vec3 Component Adders-
// Add X
template <class T>
void grid3_vector<T>::adddata_x(float data, int i)
{
	(this->grid_data)[i].x += data;
}
template <class T>
void grid3_vector<T>::adddata_x(float data, int i, int j, int k)
{
	(this->grid_data)[this->idx_3Dto1D(i,j,k)].x += data;
}
// Add Y
template <class T>
void grid3_vector<T>::adddata_y(float data, int i)
{
	(this->grid_data)[i].y += data;
}
template <class T>
void grid3_vector<T>::adddata_y(float data, int i, int j, int k)
{
	(this->grid_data)[this->idx_3Dto1D(i, j, k)].y += data;
}
// Add Z
template <class T>
void grid3_vector<T>::adddata_z(float data, int i)
{
	(this->grid_data)[i].z += data;
}
template <class T>
void grid3_vector<T>::adddata_z(float data, int i, int j, int k)
{
	(this->grid_data)[this->idx_3Dto1D(i, j, k)].z += data;
}

// Vec3 Component Getters -
// Get X
template <class T>
float grid3_vector<T>::getdata_x(int i) const
{
	return (this->grid_data)[i].x;
}
template <class T>
float grid3_vector<T>::getdata_x(int i, int j, int k) const
{
	return (this->grid_data)[this->idx_3Dto1D(i, j, k)].x;
}
// Get Y
template <class T>
float grid3_vector<T>::getdata_y(int i) const
{
	return (this->grid_data)[i].y;
}
template <class T>
float grid3_vector<T>::getdata_y(int i, int j, int k) const
{
	return (this->grid_data)[this->idx_3Dto1D(i, j, k)].y;
}
// Get X
template <class T>
float grid3_vector<T>::getdata_z(int i) const
{
	return (this->grid_data)[i].z;
}
template <class T>
float grid3_vector<T>::getdata_z(int i, int j, int k) const
{
	return (this->grid_data)[this->idx_3Dto1D(i, j, k)].z;
}


//----------------------------------------------------------------------

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


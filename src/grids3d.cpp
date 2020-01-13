// Implementation of grids2 

#include "grids2.h"

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

//----------------------------------------------------------------------\\

/*
// grid2_img Implementation - \\ WIP

grid2_img::grid2_img(const char *path, int x_s, int y_s, int idd)
	: x_size(x_s), y_size(y_s), ID(idd), img_path(path)
{
	// Size Should Be Grid Size, and thus PPM Image Size. 
	int size_1D = x_size * y_size;

	// Allocate String RawData Member Ptr - 
	raw_data = new std::string();
	raw_data->reserve(10000); // Force Max Reserve. 

	// Temp Array (Before going to Grid_Data Vector) - 
	int *temp_1d = new int[size_1D];

	// Input File Stream - 
	std::ifstream in_img(path, std::ios::binary); // In Binary Mode...

	/* ISSUES WITH TRY/CATCH EXCEPTION ENABLED READING (EVEN WITH NO FAILBIT) NoExcept Thrown, no read occurs. 
	// Enable Execeptions on ifstream obj Set Bitmask, Dont use Failbit for LineByLine Read- 
	in_img.exceptions(std::ios::badbit);

	// Catch Exceptions Throw By ifstream exeception paramters. 
	try
	{
		std::string line;
		in_img.open(path);
		
		if (in_img.is_open())
		{
			while (std::getline(in_img, line)) // Read In LineByLine.
			{
				// Do Parsing (Binary P6 PPM - P3 PPM ASCII) Operations...
				// Write into Raw_Data String Member
				raw_data->operator+= (line);
				std::cout << line << "\n";
			}	
		}
		
		in_img.close();
	}
	catch (std::ios::failure err)
	{
		std::cerr << "ERROR::Grid2_img_" << ID << "_ERROR READING IAMGE::" << err.what() << "\n";
		std::cerr << "ERROR::Grid2_img_" << ID << "_CONSTRUCTION NOT SUCESSFUL \n";
	}


	std::string line;
	if (in_img.is_open())
	{
		int line_n = 0;
		// Dont Want to Parse ASCII PPM Header , so Skip lines 0-2. 
		while (std::getline(in_img, line)) // Read In LineByLine.
		{
			// Append Write into Raw_Data String Member
			raw_data->operator+= (line);
			line_n++;
		}
	}
	in_img.close();

	// Do Parsing (Binary P6 PPM - P3 PPM ASCII Decimal Values) Operations.
	// Raw Data String Size / 3 (RGB) Should Match Img Grid Size, Bar the first 3 Line (9 Character) PPM Header.
	// Parse Final ASCII String, to 1D Temp Array - 

	for (int i = 0; i < size_1D; i++)
	{
		if (i > 9) // Ignore Header ASCII Values.
		{
			// Get rid of Negative ASCII (RGB) Values?.. Use Abs.
			*(temp_1d + i) = std::abs(int(raw_data->at(i)));
		}
		else if (i <= 9) // Set Header Values to 0.
		{
			*(temp_1d + i) = 0;
		}

	}

	// NOT RIGHT DUH 3 VALUES (RGB per PIXEL)... (Minus 9 Chars (First 3 lines of PPM Header))
	// ((img_str - 9) / 3) == Size of Img N^N (And thus 1D Array). 

	// Debug Print 1D Array.
	
	for (int i = 0; i < size_1D; i++) 
	{
		//if (i > 9) std::cout << *(temp_1d + i);
		std::cout << temp_1d[i] << "\n";
	}

	// ARRAY AND STRING SIZE DO NOT MATCH...
	std::cout << "\n";
	std::cout << "1D ARRAY SIZE = " << size_1D << "  " << "STRING SIZE = " << raw_data->size() << std::endl;
	// - 9 (PPM Header Chars) 3 Values per pixel (RGB) so String Size  - 9 ) / 3 should == the 1D Array Size (img dim ^2) 
	if (((raw_data->size() - 9) / 3) == size_1D) std::cout << "IMAGE VALUES RGB to 1D ARRAY - SIZE - TRUE \n";
	
	// Still need to debug ifstream exceptions...

	// DEBUG - Print Raw Data - 
	//print_rawdata();
}

// Grid 2 Image Destructor - 
grid2_img::~grid2_img()
{
	// DeAllocate Raw_Data and Grid_Data From Heap.
	delete raw_data; raw_data = nullptr;
	delete grid_data; grid_data = nullptr;
}

// Debug Function, Print Input Image File, Raw Data String. 
void grid2_img::print_rawdata()
{
	std::cout << *raw_data << std::endl;
}
*/
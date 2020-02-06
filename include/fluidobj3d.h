#ifndef FLUIDOBJ_H
#define FLUIDOBJ_H

// Fluid Object Header - Interface of Fluid Object, which contains Fluid Grids, Sourcing/Grid Setting Handling.

#include "grids3d.h" 

// 2D Fluid Object Class - Interface/Decleration. 

class fluidobj_3d
{

friend class fluidsolver_3;
friend class renderobject_3D;
friend class renderobject_3D_OGL; 

using BYTE = unsigned char; 

public:
	fluidobj_3d(int x_size, int y_size, int z_size, int edge_size);
	fluidobj_3d() = delete; 
	~fluidobj_3d();

	void implicit_sphere_source(float dd, const vec3<float> &vv, const vec3<float> &offset, float rad);

	/* ! TO IMPLEMENT
	// Implicit Sphere Source to Specfic Scalar Grid
	void implicit_sphere_source(grid3_scalar<float> *grid, float quant, const vec3<float> &offset, float rad);
	void implicit_sphere_source(grid3_vector<vec3<float>> *grid, const vec3<float> &quant, const vec3<float> &offset, float rad);
	*/

	// For (Density, Velocity) Grid Specfic Cell MFs - 
	void add_density(float d, int i, int j, int k);
	void add_velocity(const vec3<float> &v, int i, int j, int k);

	/* ! TO IMPLEMENT
	// Velocity Inital Condition/Overrides - 
	void radial_vel(const vec3<float> &orig_offset, float speed); // Vel Override - Radial
	void sink_vel(const vec3<float> &orig_offset, float speed); // Vel Override - Sink
	*/

	// FORCES - 
	void integrate_force(const vec3<float> &force, float dt, int i, int j, int k); // Uniform Force - PerCell. 
	/* ! TO IMPLEMENT
	void radial_force(const vec3<float> &orig_offset, float strength, float dt); // Radial Force - Vel Grid
	*/

	// UTILITYS -

	// For Swap Setting Grids (prev and cur) if not using Diffusion.
	void setcurtoprev(grid3_scalar<float> *grid0, grid3_scalar<float> *grid1);
	void setcurtoprev(grid3_vector<vec3<float>> *grid0, grid3_vector<vec3<float>> *grid1);

	void print_info();
	void writeto_txt(int frame); // DBG to view Vel grid values in ASCII. 

	// Size Members (As Reference to All Grids Created Dimnensions)
	std::size_t x_s, y_s, z_s, e_s, t_s; // X Size, Y Size, Z Size, Edge Size, Total Size. 
private:

	// SCALAR GRIDS \\

	// Scalar Grid - Density
	grid3_scalar<float> *prev_dens;
	grid3_scalar<float> *dens;

	// Scalar Grid - Fuel
	grid3_scalar<float> *prev_fuel;
	grid3_scalar<float> *fuel;

	// Scalar Grid - Heat
	grid3_scalar<float> *prev_heat;
	grid3_scalar<float> *heat;

	// Scalar Grids - Cur Step Stored Only - 

	// Scalar Grid - Collide (Marker)
	grid3_scalar<float> *col;

	// VECTOR GRIDS \\

	// Vector Grid - Velocity
	grid3_vector<vec3<float>> *prev_vel;
	grid3_vector<vec3<float>> *vel;

	// Vector Grid - (Vorticity Curl Gradient)
	grid3_vector<vec3<float>> *vc;
	grid3_vector<vec3<float>> *curl;

	// DEBUG GRIDS \\ 

	// Grids for Rendering Debugging Purposes 
	grid3_vector<vec3<float>> *preproj_vel;  
	grid3_scalar<double> *divergence, *pressure; 
};

#endif

#ifndef FLUIDOBJ_H
#define FLUIDOBJ_H

// Fluid Object Header - Interface of Fluid Object, which contains Fluid Grids, Sourcing/Grid Setting Handling.

#include "grids3d.h" 

// 2D Fluid Object Class - Interface/Decleration. 

class fluidobj_3d
{
// FluidObject is Friend of FluidSolver And RenderObj Class to allow private acess. 

friend class fluidsolver_3d;
friend class renderobject_3d;
friend class renderobject_3d_OGL; 

using BYTE = unsigned char; 

public:
	fluidobj_3d(int x_size, int y_size, int edge_size, int idd, float spc);
	fluidobj_3d() = delete; 
	~fluidobj_3d();

	// Public MFuncs to be called from Fluidsolver to modifiy grids. 
	void implicit_source(float dd, const vec3<float> &vv, const vec3<float> &offset, float rad);

	void add_density(float dens, int i, int j);
	void add_velocity(const vec3<float> &vel, int i, int j);

	// Velocity Inital Condition/Overrides - 
	void radial_vel(const vec3<float> &orig_offset, float speed); // Vel Override - Radial
	void sink_vel(const vec3<float> &orig_offset, float speed); // Vel Override - Sink

	// FORCES - 
	void integrate_force(const vec3<float> &force, float dt, int i, int j); // Uniform Force - PerCell. 
	void radial_force(const vec3<float> &orig_offset, float strength, float dt); // Radial Force - Vel Grid

	// UTILITYS -

	// For Swap Setting Grids (prev and cur) if not using Diffusion.
	void setcurtoprev(grid3_scalar<float> *grid0, grid3_scalar<float> *grid1);
	void setcurtoprev(grid3_vector<vec3<float>> *grid0, grid3_vector<vec3<float>> *grid1);

	void print_info();
	void writeto_txt(int frame); // DBG to view Vel grid values in ASCII. 

private:
	// Size Members (As Reference to All Grids Created Dimnensions)
	std::size_t x_s, y_s, z_s, e_s, t_s; // X Size, Y Size, Z Size, Edge Size, Total Size. 

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

	// Scalar Grid - Curl  
	grid3_scalar<float> *curl;

	// VECTOR GRIDS \\

	// Vector Grid - Velocity
	grid3_vector<vec3<float>> *prev_vel;
	grid3_vector<vec3<float>> *vel;

	// Vector Grid - (Vorticity Curl Gradient)
	grid3_vector<vec3<float>> *vc;

	// DEBUG GRIDS \\ 

	// Grids for Rendering Debugging Purposes 
	grid3_vector<vec3<float>> *preproj_vel;  
	grid3_scalar<double> *divergence, *pressure; 
};

#endif

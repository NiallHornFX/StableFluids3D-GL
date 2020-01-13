#ifndef FLUIDOBJ_H
#define FLUIDOBJ_H
// Fluid Object Header - Interface of Fluid Object, which contains Fluid Grids, and Input/Sourcing Handling.

#include "grids2.h" // Also Includes stdiostream,stdmemory,vec2. Saves relying on HGuards and reincluding here. 

// 2D Fluid Object Class - Interface/Decleration. 

class fluidobj_2d
{
// FluidObject is Friend of FluidSolver And RenderObj Class to allow private acess. 
friend class fluidsolver_2;
friend class renderobject_2D;
friend class renderobject_2D_OGL; 

using BYTE = unsigned char; 

public:
	fluidobj_2d(int x_size, int y_size, int edge_size, int idd, float spc);
	fluidobj_2d() = delete; 
	~fluidobj_2d();

	//static fluidobj_2d f2obj_Factory(const fluid_config &cfg); // Factor MF WIP ... 

	// Configuations of Fluid Objects Create via Factory WIP ...
	enum fluid_config
	{
		CONFIG_BASICDENSITY = 0, // Density,Velocity
		CONFIG_ZALESAK,		     // Density,Velocity
		CONFIG_IMAGE,            // Colour, Velocity
		CONFIG_COMBUSTION        // Fuel,Heat,Velocity
	};

	// Public MFuncs to be called from Fluidsolver to modifiy grids. 
	void print_info();
 
	void implicit_source(float dd, const vec2<float> &vv, const vec2<float> &offset, float rad);

	void add_density(float dens, int i, int j);
	void add_velocity(const vec2<float> &vel, int i, int j);

	// Velocity Inital Condition/Overrides - 
	void radial_vel(const vec2<float> &orig_offset, float speed); // Vel Override - Radial
	void sink_vel(const vec2<float> &orig_offset, float speed); // Vel Override - Sink

	// Forces - 
	void integrate_force(const vec2<float> &force, float dt, int i, int j); // Uniform Force - PerCell. 

	void radial_force(const vec2<float> &orig_offset, float strength, float dt); // Radial Force - Vel Grid

	// SET CUR TO PREV -
	// For Swapping Grids (prev and cur) if not using Diffusion.
	void setcurtoprev_dens();  
	void setcurtoprev_vel(); // rm these soon..

	void setcurtoprev(grid2_scalar *grid0, grid2_scalar *grid1);
	void setcurtoprev(grid2_vector *grid0, grid2_vector *grid1);

	// UTILITYS - 
	void del_dens();
	void del_vel();

	int RGB_imageLoad(const char* path); 

	// Initalize Density and Velocity for Zalesak Disk Test. 
	void zalesak_disc(float rot_speedmult);

	// For Later Possible Implmenetation -
	int writeto_img(int frame); // Currently Just Writes Density to pixel RGB. (Returns 1 if sucessful). 

	// WRITE TO IMAGE - 
	// Write Grids to PPM Images. 
	int writeto_img_pressure(int frame, bool pressureexist, grid2_scalar *pressureptr); // Needs FluidSolver to keep Pressure Field to be stored per frame.
	int writeto_img_vel(int frame);
	int writeto_img_all(int frame, grid2_scalar *pres_grid);
	int writeto_img_col();
	int ppm_imgtest(); 

	// DEBUG - 
	void writeto_txt(int frame); // DBG to view Vel grid values in ASCII. 

	// Size Members (As Reference to All Grids Created Dimnensions) Public. 
	int x_s, y_s, e_s; // X Size, Y Size, Edge Size. 
	float spacing;

private:

	// Grid Object Pointer Members - These May or May Not be allocated, based on Configuration of FluidObj. 

	// SCALAR GRIDS \\

	// Scalar Grid - Density
	grid2_scalar *prev_dens;
	grid2_scalar *dens;

	// Scalar Grid - Colour (RGB Components)
	grid2_scalar *c_R_prev, *c_G_prev, *c_B_prev;
	grid2_scalar *c_R, *c_G, *c_B; 

	// Scalar Grid - Fuel
	grid2_scalar *prev_fuel;
	grid2_scalar *fuel;

	// Scalar Grid - Heat
	grid2_scalar *prev_heat;
	grid2_scalar *heat;

	// Scalar Grids - Cur Step Stored Only - 

	// Scalar Grid - Collide (Marker)
	grid2_scalar *col;

	// Scalar Grid - Curl  
	grid2_scalar *curl;

	// VECTOR GRIDS \\

	// Vector Grid - Velocity
	grid2_vector *prev_vel;
	grid2_vector *vel;

	// Vector Grid - (Vorticity Curl Gradient)
	grid2_vector *vc;

	// DEBUG GRIDS \\ 

	// Grids for Rendering Debugging Purposes 
	grid2_vector *preproj_vel;  // Copy of CurStep Pre-Projected Velocity Field Grid. 
	grid2_scalar *divergence, *pressure; // Copy of CurStep Divergence and Pressure Grids.  

	// UTILITY MEMBERS \\

	// Store Loaded Image Data for Colour/R,G,B Scalar Grids. 
	struct image_data
	{
		int width, height, nch; 
		char *path = nullptr; 
		BYTE *img_raw; 

	}img_data;

	int ID;

};

#endif

/*	-- INFO --
	A Realtime CPU Based 3D Fluid-Solver based on Jos Stams Stable Fluids Paper. 
	With a Basic GPU Realtime RayMarcher For Rendering based in OpenGL. 

	- Support for Cubed/Uniform Grids Only for now. 
	- All Grids are currently collacated (Cell Centered) Mac/Staggerd Grids will be added soon (Vel + Pressure). 
	- Use of OpenMP for heavy Data Parallism on Solver/Grid Operations where threadsafe. All LockFree. 

	- Built and Tested on Windows x64 using MSVC. OpenMP 2.0|GLEW|GLFW. Statically Linked Currently. 
	- Use of amd64 (Intel x64) SIMD intrinsics, with hardcoded instruction sets SSE3 and AVX currently used. 
	- Modern OpenGL, Using 4.0+
	-- END --
*/

// Std Headers 
#include <iostream>
#include <vector>
#include <memory>

// Vendor Headers
#define GLEW_STATIC 
#include <GLEW\glew.h>
#include <GLFW\glfw3.h>

// Project Headers 
#include "vec3d.h"
#include "mat3d.h"
#include "grids3d.h"
#include "fluidobj3d.h"
#include "fluidsolver3d.h"
#include "rendercontext3d.h"
#include "renderobject3d.h"

// Macros - 
#define GLMajor 4
#define GLMinor 0

#define DEBUG_MODE 0 
#define USE_ASSERT 0
#define USE_SIMD 1 
#define SIMD_FORCE_INLINE

// Globals - 
short verbose = 0;
double const PI = 3.14159265359; 

int const cube = 64; // Cube Grid Size N (1-N)
int const edge = 2; // Total Edge Cells E (0 | N+1). 
int win_size_xy = 512; 
int const solve_steps = 1000; 
float const timestep = 1.0f / 60.0f; 

int main()
{
	// Create Render Context For OpenGL Context With Window Setup (in main thread). Window Dimensions Incl Edge Cell/Pixels. 
	render_context_OGL render_c (win_size_xy, win_size_xy, short(GLMajor), short(GLMinor)); 

	//renderobject_3D_OGL render_o ("OpenGL", 4, 2, vec2<int>(win_size_xy, win_size_xy), vec3<int>(cube + edge, cube + edge, cube + edge), render_c.get_window());
	//render_o.render_loop(rend_state::RENDER_DEBUG);


	// Create Fluid Object - Containing Fluid Grids and Data. 
	fluidobj_3d test_fluidobj (cube, cube, cube, edge); 

	// Create FluidSolver Instance,  Pass FluidObj Pointer to It. 
	fluidsolver_3 test_fluidsolver (&test_fluidobj, timestep);

	// Pre Solve Parmaters Inital Values Set -
	test_fluidsolver.Parms.p_Do_Dens_Diff = false; 
	test_fluidsolver.Parms.p_useVorticity = false;
	test_fluidsolver.Parms.p_Do_Dens_Disp = true; 
	test_fluidsolver.Parms.p_Do_Vel_Disp = false; 
	test_fluidsolver.Parms.p_ProjectionType = test_fluidsolver.Parms.Project_GaussSeidel_SOR; 
	test_fluidsolver.Parms.p_SOR_alpha = 1.9f;
	test_fluidsolver.Parms.p_GS_Proj_iter = 5; 
	test_fluidsolver.Parms.p_AdvectionType = test_fluidsolver.Parms.Advect_SL_BackTrace_Euler;

	// Pass Window Pointer from RenderContext - 
	test_fluidsolver.set_window(render_c.get_window());

	// Call Solve Step Loop to start simulation and rendering - 
	test_fluidsolver.solve_step(true, solve_steps);

	return 0;
}



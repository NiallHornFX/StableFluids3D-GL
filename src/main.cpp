/* INFO ---- 
	A Realtime CPU Based 3D Fluid-Solver based on Jos Stams Stable Fluids Paper. Along with Fedkiw and Bridson.
	With a Basic GPU Realtime RayMarcher For Rendering based in OpenGL. 

	- All Grids are currently collacated (Cell Centered) Mac/Staggerd Grids will be added soon (Vel + Pressure). 

	- Built and Tested on Windows x64 using MSVC. OpenMP 2.0|GLEW|GLFW. Statically Linked Currently. 
	- Use of amd64 (Intel x64) SIMD intrinsics, with hardcoded instruction sets SSE3 and AVX currently used. 
	- Modern OpenGL, Using 4.0+
	- TO-DO | CMake BuildSystem | Test on Linux 
---*/

// Std Headers 
#include <iostream>
#include <vector>
#include <memory>

// Vendor Headers
#define GLEW_STATIC // Static Linking GLEW and GLFW For now. 
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
int const edg = 2; // Total Edge Cells E (1 For each dim) (0 | N+1). N+E per dim.
int const solve_steps = 1000; 
float const timestep = 1.0f / 60.0f;  // dt = 0.0166f

/*
-- INFO --
- Fluid Object is created to Allocate Grids and Grid Settings, Which is Passed to Fluid Solver for Solving. 
- Render Context setups up Window and Graphics API Context for rendering.
- Solver contains embedded RenderObject for rendering solved FluidObj Grids to Window Framebuffer (within solvestep). 
- Inputs are handled within FluidSolver SolveStep or Encapsulated RenderObject (within solvestep). 
- Plans to Decouple Application/Simulation Loop from Solver And allow user defined App/Sim Loop to Eval Is WIP. 
-- END --
*/

int main()
{
	
	grid3_scalar<float> test(cube, cube, cube, edg);
	auto ptr = test.griddataptr_getter();

	/*
	// Create Render Context For OpenGL Context With Window Setup (in main thread). Window Dimensions Incl Edge Cell/Pixels. 
	render_context_OGL render_c (sqr + edg, sqr + edg, short(GLMajor), short(GLMinor)); 
	
	// Create Fluid Object - Containing Fluid Grids and Data. 
	fluidobj_2d test_fluidobj (sqr, sqr, edg, 0, 1.0f); // Get Rid of Spacing Control... h is always 1/N.
	//test_fluidobj.RGB_imageLoad("Checker_512_b.jpg"); // Load RGB Image

	// Print Fluid Object Debug Info - 
	//test_fluidobj.print_info();

	// Create FluidSolver Instance,  Pass FluidObj Pointer to It. 
	fluidsolver_2 test_fluidsolver (&test_fluidobj, timestep);

	// Inital Density and Velocity if used - 
	//test_fluidobj.implicit_source(0.75f, vec2<float>(0.0f, 0.0f), vec2<float>(0.5f, 0.5f), 0.01f);
	//test_fluidobj.sink_vel(vec2<float>(0.5f, 0.49999f), 0.1f);
	//test_fluidobj.radial_vel(vec2<float> (0.5f, 0.5f), 5.0f);
	
	// Pre Solve Parmaters Inital Values Set. (Can be Changed In Solve Per Step Later) - 
	//test_fluidsolver.Parms.p_useColour = true;

	test_fluidsolver.Parms.p_Do_Dens_Diff = false; 
	test_fluidsolver.Parms.p_useVorticity = false;
	test_fluidsolver.Parms.p_Do_Dens_Disp = false; 

	test_fluidsolver.Parms.p_ProjectionType = test_fluidsolver.Parms.Project_GaussSeidel_SOR; 
	//test_fluidsolver.Parms.p_ProjectionType = test_fluidsolver.Parms.Project_Jacobi;
	//test_fluidsolver.Parms.p_Jacobi_Proj_Iter = 100; 
	test_fluidsolver.Parms.p_SOR_alpha = 1.9f;
	test_fluidsolver.Parms.p_GS_Proj_iter = 10; 
	test_fluidsolver.Parms.p_AdvectionType = test_fluidsolver.Parms.Advect_SL_BackTrace_Euler;
	test_fluidsolver.Parms.p_InteroplationType = test_fluidsolver.Parms.Interoplation_Cosine;
	//test_fluidsolver.Parms.p_useVorticity = true; 

	// Disable Project/Advect for Pref testing. 
//	test_fluidsolver.Parms.p_AdvectionType = test_fluidsolver.Parms.Advect_NONE_DBG;
//	test_fluidsolver.Parms.p_ProjectionType = test_fluidsolver.Parms.Project_NONE_DBG;

	// For now Pass Window Pointer directly from RenderContext (Oppose to RenderContext Obj/ptr Itself) - 
	test_fluidsolver.set_window(render_c.get_window());

	// Call Solve Step Loop to start simulation and rendering - 
	test_fluidsolver.solve_step(true, false, false, 0.001f, 0.001f, 10, 8, solve_steps);
	*/

	return 0;
}



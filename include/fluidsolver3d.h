#ifndef FLUIDSOLVER_2_H
#define FLUIDSOLVER_2_H

#include "fluidobj3d.h"
#include "renderobject3d.h" // Incl Here not source, to stop Inclusion of Glew within FluidSolversrc.

// Vendor Headers - 
#include <GLFW\glfw3.h>

// Std Headers 
#include <immintrin.h> // SIMD Intrinsics. 
#include <functional> // std::function object. 
#include <iostream>

using ushort = unsigned short; 
using avx256 = __m256;

// SIMD - Free Function Utilties -

// SIMD - 8 Float (m256) Addition Together. 
__forceinline float simd256_hAdd(__m256 a) 
{
	__m256 t1 = _mm256_hadd_ps(a, a);
	__m256 t2 = _mm256_hadd_ps(t1, t1);
	__m128 t3 = _mm256_extractf128_ps(t2, 1);
	__m128 t4 = _mm_add_ss(_mm256_castps256_ps128(t2), t3);
	return _mm_cvtss_f32(t4);
}

// SIMD - 4 Float (m128) Addition Together. 
__forceinline float simd128_hAdd(__m128 a)
{
	__m128 t1 = _mm_hadd_ps(a, a);
	__m128 t2 = _mm_hadd_ps(t1, t1);
	//return t2.m128_f32[0];
	return _mm_cvtss_f32(t2);
}

enum class step
{
	STEP_CUR,
	STEP_PREV 
};

// Solver Utils Container, Useful static - Lambdas, Functors, FreeFuncs, Enums etc. 
struct solver_utils
{
	solver_utils() = delete; 

	// Value Remapping - 
	static std::function<float(float, float, float)> clamp;
	static std::function<float(float, float, float, float, float)> fit;

	// Interoplation - 
	// 1D
	static std::function<float(float, float, float)> lerp;
	static std::function<float(float, float, float)> cosinterp; 

};


// 2D Fluid Solver Class - Interface/Decleration. 

class fluidsolver_3
{
public:

	// Enforce no default constructor/implicit construction. FluidObject Pointer Is NEEDED for Construction.
	fluidsolver_3(fluidobj_3d *f3dptr, float dtt);
	~fluidsolver_3();
	fluidsolver_3() = delete; // Delete Default Constructor. 

	// Setter/Getters - 
	void set_window(GLFWwindow *win);

	// Main Solver Step - (Only Public Fluid Solve MFunc). 
	void solve_step(bool solve, bool do_diffdens, bool do_diffvel, float dens_diff, float vel_diff, int proj_iter, int diff_iter, int max_step);

	// PARAMTERS \\  WIP - 
	struct FluidSolver2_Paramters
	{
		enum ProjType
		{
			Project_GaussSeidel = 0,
			Project_GaussSeidel_SOR,
			Project_Jacobi,
			Project_NONE_DBG
		};

		enum AdvType
		{
			Advect_SL_BackTrace_Euler = 0,
			Advect_SL_BackTrace_RK2,
			Advect_NONE_DBG

			// Interpolant Variants
			// BFECC,MacCormack ... 
		};

		enum InterpType
		{
			Interoplation_Linear = 0, 
			Interoplation_Cosine,
			Interoplation_Cubic
		};

		// Using Default inclass Initalization
		bool p_useColour = false, p_useVorticity = false, p_useFuel = false, p_useHeat = false; // Optional Grids, for Factory prep. 

		// Projection Paramters - 
		ProjType p_ProjectionType = Project_GaussSeidel_SOR;
		int p_GS_Proj_iter = 10, p_Jacobi_Proj_Iter = 50;
		float p_SOR_alpha = 1.8f; // 1.0-2.0
		bool p_Do_PreAdvectProject = false;

		// Advection Paramters - 
		AdvType p_AdvectionType = Advect_SL_BackTrace_Euler;
		InterpType p_InteroplationType = Interoplation_Linear; 

		// Diffusion Switches - 
		bool p_Do_Vel_Diff = false, p_Do_Dens_Diff = true;
		int p_Vel_Diff_iter = 5, p_Dens_Diff_iter = 5;
		float p_Vel_Diffuse_Str = 0.0f, p_Dens_Diffuse_Str = 0.001f;

		// Dissipation Swtiches - 
		bool p_Do_Dens_Disp = true, p_Do_Vel_Disp = false;
		float p_Dens_Disp_Mult = 0.980f, p_Vel_Disp_Mult = 1.0f; 

		bool p_spherebounds_killint = false; 

		//p_dt // Constructor bound.
		float emit_sphere_rad, emit_sphere_surfthresh;
		float col_sphere_rad = 0.005f, col_sphere_surfthresh = 0.003f;

		vec2<float> emit_sphere_offset, col_sphere_offset;

		bool p_MT_Global = true, p_MT_Bounds = true, MT_DivergenceCalc, MT_PressureGradSubtraction; // MultiThread Switches. 

	} Parms;

protected:

	// BOUNDARY CONDTIONS \\ - 

	// ! Generic Boundary MFuncs
	// Edge Bounds - 
	void edge_bounds(grid3_scalar<float> *grid);
	void edge_bounds(grid3_vector<vec3<float>> *grid);

	// Sphere/sphere Bounds - 
	void sphere_bounds_set(float radius, float col_iso, const vec3<float> &offset);  // Set SphereBounds SDF Grid. 
	void sphere_bounds_eval(grid3_scalar<float> *grid, float col_iso);
	void sphere_bounds_eval(grid3_vector<vec3<float>> *grid, float col_iso);

	// DIFFUSION \\ - 
	// Gauss-Seidel Relaxation Diffusion Based on Grid Neighbours 

	// ! Gauss-Seidel Relaxation - Diffusion - 
	void diffuse(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1, float diff, ushort iter);
	void diffuse(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1, float diff, ushort iter);

	// ADVECTION \\ - 
	// Semi-Lagraginin Advection Step Single Backtrace Step. Overload Based on Grid Type. 

	// Semi Lagrangian (Single Backwards Euler) Advection - 
	void advect_sl(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1);
	void advect_sl(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1);

	// Semi Lagrangian RungeKutta 2 - MidPoint Advection.  
	void advect_sl_mp(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1);
	void advect_sl_mp(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1);

	// PROJECTION \\ - 
	// Velocity Field Projection/Pressure Solve - 

	// Projection - Gauss-Seidel Relaxation Method, with SOR (Single-Threaded)
	void project(int iter);
	void project_SIMD(int iter);

	void project_GS_RB(int iter); // Gauss-Seidel RedBlack (MultiThreaded). 

	// Projection - Jacobi Solve (Multi-Threaded)
	void project_jacobi(int iter); 

	// DISSIPATION \\ - 
	void dissipate(grid3_scalar<float> *grid, float disp_mult, float dt);
	void dissipate(grid3_vector<vec3<float>> *grid, float disp_mult, float dt);

	// MISC \\ - 
	//void vel_force(vec3<float> ff, float dtt);
	//// Run Passed Lambda Callback Over Each Vector Grid Cell. 
	//void custom_force(grid3_vector<vec3<float>> *grid, std::function <vec2<float>(vec2<int> idx)> &force);

	// VORTICITY CONFINEMENT WIP \\ - 
	//void vorticty_confine(float strength);
	//void vorticty_confine_otf(float strength); // On the fly, without Vort/Curl Grids. 
	//void vorticity_confine_B(float strength);

	// SOLVER - SUB SOLVERS \\ - 
	void density_step();
	void velocity_step();

	void debug_step(); 

	// UTILITY \\ - 

	// Mouse Input - 
	void updt_mousepos(const step step_id); // 0 - N XY Mouse Coords. 
	void updt_mouseposNorm(const step step_id); // 0 -1 XY Mouse Coords
	void updt_mouseposRange(const step step_id); // -1 to 1. XY Mouse Coords (Incorrect, Depreacted). 
	void updt_mousevel();
	void sphere_rad_test();

	// Sourcing/Collison Setters ... 
	//void set_sourcepos(); // WIP Set in FluidObj would still need passing?
	//void set_spherebounds_pos(vec2 &offset);

	// TEMP / DBG \\ - 
	void fill_test(int x);  // DBG - Iteration test.

	// Manually Delete Temp Grids
	void del_pressure();
	void del_divergence();

private:

	fluidobj_3d *f3obj; // Pointer to FluidObject3D this FluidSolver Instance is Solving/Operating on.
	float dt;

	// Hard Coded on Construction (Constant) Not Paramters - (Get from FluidObj Ideally, Kinda Dont need these do we?)
	// Size Members, Get from Fluid Object. 

	int x_s, y_s, z_s, e_s; // Axis Sizes and Edge Cell Size/Count.
	std::size_t total_size; // Total Cell Size (X+E) * (Y+E). 

	// Assume Square Grid - 
	int N_dim; // Size of SINGLE Dimension w/o Edge Cells.
	int NE_dim; // Size of SINGLE Dimension w/ Edge Cells.

	// Input - Mouse Members/Functions - 
	double xpos_0, ypos_0, xpos_0_N, ypos_0_N, xpos_0_R, ypos_0_R; // Prev_Step Mouse Data
	double xpos_1, ypos_1, xpos_1_N, ypos_1_N, xpos_1_R, ypos_1_R;  // Cur_step Mouse Data 
	vec2<float> mouse_vel;

	// Temp/Scratch Grids Ptrs for Solver Use Only (Never Directly Rendered or used back in FluidObj unless copied).
	grid3_scalar<float> *pressure, *pressure_1;
	grid3_scalar<float> *divergence;
	grid3_scalar<float> *vort; // 2D Vort/Curl Grid. 
	grid3_scalar<float> *spherebounds_sdf; // Grid to hold SphereBounds SDF Value For Cur Step. 

	// Render Members - 
	renderobject_3D_OGL *render_obj; // Uses OGL Not Base Ptr for now, nocast. 
	GLFWwindow *winptr; // From Render Contex Passed In Window. 
};

#endif
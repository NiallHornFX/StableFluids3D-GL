#ifndef FLUIDSOLVER_2_H
#define FLUIDSOLVER_2_H

// Project Headers
#include "fluidobj3d.h"
#include "renderobject3d.h" 

// Vendor Headers - 
#include <GLFW\glfw3.h>

// Std Headers 
#include <functional> 
#include <algorithm>
#include <iostream>
#include <tuple>
#include <cmath>

#include <immintrin.h> 

typedef std::tuple<vec3<float>, vec3<float>> vec3f_tuple; 
typedef std::tuple<float, float> f_tuple; 

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


// 3D Fluid Solver Class - Interface/Decleration. 

class fluidsolver_3
{
public:

	fluidsolver_3(fluidobj_3d *f3dptr, float dtt);
	~fluidsolver_3();
	fluidsolver_3() = delete; // Delete Default Constructor. 

	void set_window(GLFWwindow *win);

	void solve_step(bool solve, int max_step);

	// PARAMTERS \\  WIP - 
	// Paramters Use Inclass Init for Default Values. 
	struct FluidSolver3_Paramters
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
			Advect_MC_Euler,
			Advect_NONE_DBG
		};

		const char* AdvType_Key[4] = 
		{
			"Advect_SL_BackTrace_Euler",
			"Advect_SL_BackTrace_RK2",
			"Advect_MacCormack_Euler",
			"Advect_NONE_DBG"
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
		InterpType p_InteroplationType = Interoplation_Cosine; 
		float p_McC_LimiterStrength = 0.25f; 

		// Diffusion Switches - 
		bool p_Do_Vel_Diff = false, p_Do_Dens_Diff = true;
		int p_Vel_Diff_iter = 5, p_Dens_Diff_iter = 5;
		float p_Vel_Diffuse_Str = 0.0f, p_Dens_Diffuse_Str = 0.001f;

		// Dissipation Swtiches - 
		bool p_Do_Dens_Disp = true, p_Do_Vel_Disp = false;
		float p_Dens_Disp_Mult = 0.980f, p_Vel_Disp_Mult = 1.0f; 

		bool p_spherebounds_killint = false; 

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
	// ! Gauss-Seidel Relaxation - Diffusion - 
	void diffuse(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1, float diff, ushort iter);
	void diffuse(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1, float diff, ushort iter);

	// ADVECTION \\ - 

	// Semi Lagrangian (Single Backwards Euler) Advection - 
	void advect_sl(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1);
	void advect_sl(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1);

	// MacCormack (Single Euler) Backwards,Forwards Advection - 
	void advect_mc(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1);
	void advect_mc(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1);

	// COMBUSTION \\ - 
	// TBD. 

	// PROJECTION \\ - 
	// Projection - Gauss-Seidel Relaxation Method, with SOR (Single-Threaded)
	void project(int iter);
	// Projection - Jacobi Solve (Multi-Threaded)
	void project_jacobi(int iter); 

	// DISSIPATION \\ 
	void dissipate(grid3_scalar<float> *grid, float disp_mult, float dt);
	void dissipate(grid3_vector<vec3<float>> *grid, float disp_mult, float dt);

	// VORTICITY CONFINEMENT WIP \\ - 
	//void vorticty_confine(float strength);

	// SOLVER - SUB SOLVERS \\ - 
	void density_step();
	void combustion_step(); // WIP
	void velocity_step();

	// UTILITY \\ - 
	// Vel Additon - 
	void add_velocity(const vec3<float> &vel);
	// Mouse Input - 
	void updt_mousepos(const step step_id); // 0 - N XY Mouse Coords. 
	void updt_mouseposNorm(const step step_id); // 0-1 XY Mouse Coords
	void updt_mouseposRange(const step step_id); // -1 to 1. 
	void updt_mousevel();
	// Advection Interoplation - 
	float       interp_scalarCosine(grid3_scalar<float> *grid_0, int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) const;
	float       interp_scalarLinear(grid3_scalar<float> *grid_0, int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) const;
	vec3<float> interp_vectorLinear(grid3_vector<vec3<float>> *grid_0, int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) const;
	vec3<float> interp_vectorCosine(grid3_vector<vec3<float>> *grid_0, int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) const;
	f_tuple minmaxcell_scalar(grid3_scalar<float> *grid_0, int i0, int i1, int j0, int j1, int k0, int k1) const;
	vec3f_tuple minmaxcell_vector(grid3_vector<vec3<float>> *grid_0, int i0, int i1, int j0, int j1, int k0, int k1) const;

	// Misc Uitl - 
	void sphere_rad_test();

	// STATIC UTILITY \\ 
	__forceinline static float fitRange(float val, float a_min, float a_max, float b_min, float b_max);
	__forceinline static float clamp(float val, float min, float max);
	__forceinline static float lerp(float val_0, float val_1, float bias);
	__forceinline static float cosinterp(float val_0, float val_1, float bias);
	__forceinline static vec3<float> vec_clamp(const vec3<float> &v_min, const vec3<float> &v_max, const vec3<float> v);
	__forceinline static vec3<float> vec_lerp(const vec3<float> &v_a, const vec3<float> &v_b, float bias);
	// Grid-Index-Grid Space Conversion 
	__forceinline static vec3<float> idx_indexToGrid(int i, int j, int k, int N_dim); 
	__forceinline static vec3<float> idx_gridToIndex(float x, float y, float z, int N_dim);

	// Delete Temp Grids
	void del_pressure();
	void del_divergence();

private:

	fluidobj_3d *f3obj; 
	float dt;

	int x_s, y_s, z_s, e_s; 
	std::size_t total_size; 

	// Assume Cubed Grid - 
	int N_dim; 
	int NE_dim; 

	// Input - Mouse Members/Functions - 
	double xpos_0, ypos_0, xpos_0_N, ypos_0_N, xpos_0_R, ypos_0_R; // Prev_Step Mouse Data
	double xpos_1, ypos_1, xpos_1_N, ypos_1_N, xpos_1_R, ypos_1_R;  // Cur_step Mouse Data 
	vec2<float> mouse_vel;

	// Temp/Scratch Grids for Solver Use Only.
	grid3_scalar<float> *pressure, *pressure_1;
	grid3_scalar<float> *divergence;
	grid3_scalar<float> *vort; // 
	grid3_scalar<float> *spherebounds_sdf; 

	// Render Members - 
	renderobject_3D_OGL *render_obj; // HC OGL.
	GLFWwindow *winptr; 
};

// INLINE (__forceinline) fluidsolver_3 Utility Member Functions - 

// LERP: From One Float Value to Another By 0-1 Bias Value. 
float fluidsolver_3::lerp(float val_0, float val_1, float bias)
{
	// Using FMAs Is Even Slower. 
	//return std::fmaf(bias, val_0, std::fmaf(-bias, val_1, val_1));

	return (1.0f - bias) * val_0 + bias * val_1;
}

// CLAMP: Between Min and Max Range 
float fluidsolver_3::clamp(float val, float min, float max)
{
	return std::max(min, std::min(val, max));
}

// FIT RANGE: Value From Some Range A to Some New Range B (Eg 0 to 1, to -1 to 1 0.5 becomes 0)
float fluidsolver_3::fitRange(float val, float a_min, float a_max, float b_min, float b_max)
{
	return b_min + (val - a_min)*(b_max - b_min) / (a_max - a_min);
}

// COSINTERP 1D Cosine Interpolation. 
float fluidsolver_3::cosinterp(float val_0, float val_1, float bias)
{
	float mu = (1.0f - std::cos(bias*PI)) / 2.0f;
	return (float)(val_0*(1.0f - mu) + val_1 * mu);
}

// VEC_CLAMP: Clamp Vector v components to Min/Max Vectors components.
vec3<float> fluidsolver_3::vec_clamp(const vec3<float> &v_min, const vec3<float> &v_max, const vec3<float> v)
{
	vec3<float> temp = v;
	temp.x = std::max(v_min.x, std::min(v.x, v_max.x));
	temp.y = std::max(v_min.y, std::min(v.y, v_max.y));
	temp.z = std::max(v_min.z, std::min(v.z, v_max.z));
	return temp;
}

// VEC_LERP: Linear Interoplate Between v_a and v_b components by Bias. 
vec3<float> fluidsolver_3::vec_lerp(const vec3<float> &v_a, const vec3<float> &v_b, float bias)
{
	float xx = (1.0f - bias) * v_a.x + bias * v_b.x;
	float yy = (1.0f - bias) * v_a.y + bias * v_b.y;
	float zz = (1.0f - bias) * v_a.z + bias * v_b.z;
	return vec3<float>(xx, yy, zz);
}

// Space-Index Conversion MFuncs \\ 
vec3<float> fluidsolver_3::idx_indexToGrid(int i, int j, int k, int N_dim)
{
	return vec3<float>((float)i / (float)N_dim, (float)j / (float)N_dim, (float)k / (float)N_dim);
}

// Return as vec3<float> to keep fractional component intact for Interp coefficents. 
vec3<float> fluidsolver_3::idx_gridToIndex(float x, float y, float z, int N_dim)
{
	float N_dim_f = (float)N_dim;
	return vec3<float>((x * N_dim_f), (y * N_dim_f), (z * N_dim_f));
};

#endif
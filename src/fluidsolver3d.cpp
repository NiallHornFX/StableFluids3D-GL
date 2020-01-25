// Implementation of fluidsolver_3
#include "fluidsolver3d.h"

// C++ Headers - 
#include <memory>
#include <chrono>
#include <thread>
#include <sstream>
#include <algorithm>
#include <cassert>

// Windows Specfic Headers - 
#include <omp.h> // OpenMP 2.0

#define dospherebound 1 // Enable Sphere Collisions
#define JACOBI_DO_SOR 1 // Do Sucessive Over Relaxiation Inside Jacobi Projection.
#define DO_SPHEREBOUNDS_MT 1

#define RENDER_GL 1
#define RENDER_IMG 0

// Replace with Parms ... 
// Sphere Bound Globals (For Sphere Bound MFunc Paramters) Oppose to Pass Per Nested MFunc Call. 
vec3<float> spherebound_offset(0.0f, 0.0f, 0.0f);  
const float spherebound_coliso = 0.0025f; 
float spherebound_radius = 0.005f; 

extern short verbose; // Get Verbose Global from main. 

// !CLEAN
// !TODO 
// !MT MultiThreading Disabled Currently. 
// !SIMD SIMD Disabeled Currently.


/* ====================================================
	fluidsolver_3 CONSTRUCTOR/DESTRUCTOR'S
   ==================================================== */

// fluidsolver_3 CONSTRUCTOR (Explict Only) - 
// Added dt Param & HardCoded Spherebound Initalization List. 

fluidsolver_3::fluidsolver_3(fluidobj_3d *f3dptr, float dtt) 
	: dt(dtt), f3obj(f3dptr) //,spherebound_coliso(0.002f), spherebound_radius(0.01f), spherebound_offset(-0.5f, 0.75f)
{
	// Get Dimensions,Spacing from FluidObject Ensure Values Match. 
	x_s = f3obj->x_s, y_s = f3obj->y_s, z_s = f3obj->z_s, e_s = f3obj->e_s;
	total_size = f3obj->t_s;

	N_dim =  f3dptr->x_s; // Single Dim Size, No Edge. 
	NE_dim = f3dptr->x_s + f3dptr->e_s; // Single Dim Size, With Edge Size. 

	// Init Solver (non temp) grids.
	spherebounds_sdf = new grid3_scalar<float>(N_dim, N_dim, N_dim, 2);

}

// fluidsolver_3 DESTRUCTOR 
fluidsolver_3::~fluidsolver_3()
{
	delete spherebounds_sdf; spherebounds_sdf = nullptr; 

	del_pressure();  del_divergence();

	if (vort || vort != nullptr)
	{
		delete vort; vort = nullptr; 
	}

	delete render_obj; render_obj = nullptr; 

	// FluidObj Is Responsible for its own deletion. Do Not Delete Here. 
}

/* ====================================================
	MANUAL USER FIELD DELTETION
   ==================================================== */

void fluidsolver_3::del_pressure()
{
	// If Initalized, or not set to nullptr, then delete what there pointing to. 
	if (pressure || pressure != nullptr)
	{
		delete pressure;
		pressure = nullptr;
	}

	// if P_useJacobi ...?
	if (pressure_1 || pressure_1 != nullptr)
	{
		delete pressure_1;
		pressure_1 = nullptr;
	}
}

void fluidsolver_3::del_divergence()
{
	// If Initalized, or not set to nullptr, then delete what there pointing to. 
	if (divergence || divergence != nullptr)
	{
		delete divergence;
		divergence = nullptr;
	}
}

/* ====================================================
	VELOCITY AND DENSITY - BOUNDARY CONDTIONS 
   ==================================================== */

// Explcitlly Setting Solid Edge Cell Boundary Condtions, on Edge/Ghost Cells. 

// fluidsolver_3::edge_bounds INFO - 

// Set Edge Cell Reflection of XYZ Velocity And Average Corner Cells of Vel and Dens Field. 

// ** EDGE_BOUNDS_SCALAR-FIELD IMPLEMENTATION ** \\ - 
// Assumes Grid is Cubed (N = X = Y = Z)
void fluidsolver_3::edge_bounds(grid3_scalar<float> *grid)
{
	// bool do_mt = Parms.p_MT_Global & Parms.p_MT_Global;

	// 6 Edge Faces of Grid (X-|X+|Y-|Y+|Z-|Z+)

	//!MT #pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 1; i <= N_dim; i++)
	{
		// X- Edge Boundary 
		float x0i = grid->getdata(1, i, i); // [0,i,i] Boundary Values from Edge [1,i,i] Values.
		grid->setdata(x0i, 0, i, i);

		// X+ Edge Boundary
		float xNi = grid->getdata(N_dim, i, i); // [N+1,i,i] Boundary Values from Edge [N,i,i] Values.
		grid->setdata(xNi, (N_dim + 1), i, i);

		// Y- Edge Boundary
		float y0i = grid->getdata(i, 1, i); // [i, 0, i] Boundary Values from Edge [i, 1, i] Values.
		grid->setdata(y0i, i, 0, i);

		// Y+ Edge Boundary
		float yNi = grid->getdata(i, N_dim, i); // [i, N+1, i] Boundary Values from Edge [i,N,i] Values.
		grid->setdata(yNi, i, N_dim + 1, i);

		// Z- Edge Boundary
		float z0i = grid->getdata(i, i, 1); // [i, i, 0] Boundary Values from Edge [i,i,1] Values.
		grid->setdata(z0i, i, i, 0);

		// Z+ Edge Boundary
		float zNi = grid->getdata(i, i, N_dim); // [i, i, N+1] Boundary Values from Edge [i, i, N] Values.
		grid->setdata(zNi, i, i, N_dim + 1);
	}

	// 8 Corner Cells, ScalarGrid Edge Bounds Corner Adjacent Cell Neighbour Averages -
	// Self + or - XYZ (0(+1) or N+1(-1(N))) (0 and N+1 Axis Ghost/Edge Cells).
	// If At 0 For Coord Offset +1, if At N+1 For Coord Offset -1. Offset Axis, Keep Others Const. Like a PDE. 

	// Corner Cell = Adjacent X + Y + Z (LH CoordSys Z->Back)
	// 0, 0, 0 = 1, 0, 0 | 0, 1, 0 | 0, 0, 1 // LeftBottomFront (Assume Origin)
	// 0,N+1,0 = 1,N+1,0 | 0,N,1 | 0,N+1,1 // LeftTopFront
	// 0, 0, N+1 = 1, 0, N+1 | 0, 1, N+1 | 0, 0, N //LeftBottomBack
	// 0,N+1,N+1 = 1,N+1,N+1| 0,N,N+1 | 0,N+1,N // LeftTopBack
	// N+1,0,0 = N,0,0 | N+1,1,0 | N+1,0,1 // RightBottomFront
	// N+1,N+1,0 = N,N+1,0 | N+1,N,0 | N+1,N+1,1 // RightTopFront
	// N+1,0,N+1 = N, 0, N+1 | N+1,1,N | N+1,0,N // RightBottomBack
	// N+1,N+1,N+1 = N, N+1, N+1 | N+1,N,N+1 | N+1,N+1,N // RightTopBack

	// 3D 0,0,0 = 1,0,0 + 0,1,0 + 0,0,1 
	float c_0_0_0 = 0.33f * (grid->getdata(1, 0, 0) + grid->getdata(0, 1, 0) + grid->getdata(0, 0, 1));
	grid->setdata(c_0_0_0, 0, 0, 0);

	// 3D 0,N+1,0 = 1,N+1,0 + 0,N,1 + 0,N+1,1
	float c_0_N1_0 = 0.33f * (grid->getdata(1, N_dim+1, 0) + grid->getdata(0, N_dim, 1) + grid->getdata(0, N_dim+1, 1));
	grid->setdata(c_0_N1_0, 0, N_dim + 1, 0); 

	// 3D 0,0,N+1 = 1,0,N+1 + 0,1,N+1 + 0,0,N
	float c_0_0_N1 = 0.33f * (grid->getdata(1, 0, N_dim+1) + grid->getdata(0, 1, N_dim+1) + grid->getdata(0, 0, N_dim));
	grid->setdata(c_0_0_N1, 0, 0, N_dim+1);

	// 3D 0,N+1,N+1 = 1,N+1,N+1 + 0,N,N+1 + 0,N+1,N
	float c_0_N1_N1 = 0.33f * (grid->getdata(1, N_dim + 1, N_dim + 1) + grid->getdata(0, N_dim, N_dim + 1) + grid->getdata(0, N_dim + 1, N_dim));
	grid->setdata(c_0_N1_N1, 0, N_dim+1, N_dim+1);

	// 3D N+1,0,0 = N,0,0 + N+1,1,0 + N+1,0,1
	float c_N1_0_0 = 0.33f * (grid->getdata(N_dim, 0, 0) + grid->getdata(N_dim + 1, 1, 0) + grid->getdata(N_dim + 1, 0, 1));
	grid->setdata(c_N1_0_0, N_dim+1, 0, 0); 

	// 3D N+1,N+1,0 = N,N+1,0 + N+1,N,0 + N+1,N+1,1
	float c_N1_N1_0 = 0.33f * (grid->getdata(N_dim, N_dim + 1, 0) + grid->getdata(N_dim + 1, N_dim, 0) + grid->getdata(N_dim + 1, N_dim + 1, 1));
	grid->setdata(c_N1_N1_0, N_dim+1, N_dim+1, 0);

	// 3D N+1,0,N+1 = N,1,N+1 + N+1,1,N+1 + N+1,0,N
	float c_N1_0_N1 = 0.33f * (grid->getdata(N_dim, 0, N_dim + 1) + grid->getdata(N_dim + 1, 1, N_dim + 1) + grid->getdata(N_dim + 1, 0, N_dim));
	grid->setdata(c_N1_0_N1, N_dim+1, 0, N_dim+1);

	// 3D N+1,N+1,N+1 = N,N+1,N+1 + N+1,N,N+1 + N+1,N+1,N
	float c_N1_N1_N1 = 0.33f * (grid->getdata(N_dim, N_dim + 1, N_dim + 1) + grid->getdata(N_dim + 1, N_dim, N_dim + 1) + grid->getdata(N_dim + 1, N_dim + 1, N_dim)); 
	grid->setdata(c_N1_N1_N1, N_dim + 1, N_dim + 1, N_dim + 1); 



}

// ** EDGE_BOUNDS_VECTOR-FIELD IMPLEMENTATION ** \\ - 

// Is ThreadSafe as Face/Edge Cells and Corner Cells only Read and Write to Selves locally.
void fluidsolver_3::edge_bounds(grid3_vector<vec3<float>> *grid)
{
	// Set Grid Edge Bounds Faces, to reflect the perpendicualr velocity component for cells within each six faces of grid. 

	// X +/- Grid Face Cells = Y Compoent Reflection
	// Y +/- Grid Face Cells = X Component Relfection
	// Z +/- Grid Face Cells = Z Component Relfection (Yes Z reflects Z)
	
	// !MT #pragma omp parallel for num_threads(omp_get_max_threads())

	for (int i = 1; i <= N_dim; i++)
	{
		// X -/+ Face Cells, Reflect Y Velocity Component -

		// X- Edge Boundary 
		float x0i = grid->getdata_y(1, i, i); // [0,i,i] Boundary Values from Edge [1,i,i] Values.
		x0i *= -1.0f; 
		grid->setdata_y(x0i, 0, i, i);

		// X+ Edge Boundary
		float xNi = grid->getdata_y(N_dim, i, i); // [N+1,i,i] Boundary Values from Edge [N,i,i] Values.
		xNi *= -1.0f;
		grid->setdata_y(xNi, (N_dim + 1), i, i);

		// Y -/+ Face Cells, Reflect X Velocity Component -

		// Y- Edge Boundary
		float y0i = grid->getdata_x(i, 1, i); // [i, 0, i] Boundary Values from Edge [i, 1, i] Values.
		y0i *= -1.0f; 
		grid->setdata_x(y0i, i, 0, i);

		// Y+ Edge Boundary
		float yNi = grid->getdata_x(i, N_dim, i); // [i, N+1, i] Boundary Values from Edge [i,N,i] Values.
		yNi *= -1.0f; 
		grid->setdata_y(yNi, i, N_dim + 1, i);

		// Z -/+ Face Cells, Reflect Z Velocity Component -

		// Z- Edge Boundary
		float z0i = grid->getdata_z(i, i, 1); // [i, i, 0] Boundary Values from Edge [i,i,1] Values.
		z0i *= -1.0f; 
		grid->setdata_z(z0i, i, i, 0);

		// Z+ Edge Boundary
		float zNi = grid->getdata_z(i, i, N_dim); // [i, i, N+1] Boundary Values from Edge [i, i, N] Values.
		zNi *= -1.0f; 
		grid->setdata_z(zNi, i, i, N_dim + 1);
	}

	// Corner Cell Interoplation/Averaging of Adjacent Neighbours for Ux,Vy,Wz Velocity Components -  
	// Do Each Corner Cell For each Component. Inline Averaging (messy but perf !). 

	// U (X Component)
	float cU_0_0_0 = 0.33f * (grid->getdata_x(1, 0, 0) + grid->getdata_x(0, 1, 0) + grid->getdata_x(0, 0, 1));
	grid->setdata_x(cU_0_0_0, 0, 0, 0);
	float cU_0_N1_0 = 0.33f * (grid->getdata_x(1, N_dim + 1, 0) + grid->getdata_x(0, N_dim, 1) + grid->getdata_x(0, N_dim + 1, 1));
	grid->setdata_x(cU_0_N1_0, 0, N_dim + 1, 0);
	float cU_0_0_N1 = 0.33f * (grid->getdata_x(1, 0, N_dim + 1) + grid->getdata_x(0, 1, N_dim + 1) + grid->getdata_x(0, 0, N_dim));
	grid->setdata_x(cU_0_0_N1, 0, 0, N_dim + 1);
	float cU_0_N1_N1 = 0.33f * (grid->getdata_x(1, N_dim + 1, N_dim + 1) + grid->getdata_x(0, N_dim, N_dim + 1) + grid->getdata_x(0, N_dim + 1, N_dim));
	grid->setdata_x(cU_0_N1_N1, 0, N_dim + 1, N_dim + 1);
	float cU_N1_0_0 = 0.33f * (grid->getdata_x(N_dim, 0, 0) + grid->getdata_x(N_dim + 1, 1, 0) + grid->getdata_x(N_dim + 1, 0, 1));
	grid->setdata_x(cU_N1_0_0, N_dim + 1, 0, 0);
	float cU_N1_N1_0 = 0.33f * (grid->getdata_x(N_dim, N_dim + 1, 0) + grid->getdata_x(N_dim + 1, N_dim, 0) + grid->getdata_x(N_dim + 1, N_dim + 1, 1));
	grid->setdata_x(cU_N1_N1_0, N_dim + 1, N_dim + 1, 0);
	float cU_N1_0_N1 = 0.33f * (grid->getdata_x(N_dim, 0, N_dim + 1) + grid->getdata_x(N_dim + 1, 1, N_dim + 1) + grid->getdata_x(N_dim + 1, 0, N_dim));
	grid->setdata_x(cU_N1_0_N1, N_dim + 1, 0, N_dim + 1);
	float cU_N1_N1_N1 = 0.33f * (grid->getdata_x(N_dim, N_dim + 1, N_dim + 1) + grid->getdata_x(N_dim + 1, N_dim, N_dim + 1) + grid->getdata_x(N_dim + 1, N_dim + 1, N_dim));
	grid->setdata_x(cU_N1_N1_N1, N_dim + 1, N_dim + 1, N_dim + 1);
	// V (Y Component)
	float cV_0_0_0 = 0.33f * (grid->getdata_y(1, 0, 0) + grid->getdata_y(0, 1, 0) + grid->getdata_y(0, 0, 1));
	grid->setdata_y(cV_0_0_0, 0, 0, 0);
	float cV_0_N1_0 = 0.33f * (grid->getdata_y(1, N_dim + 1, 0) + grid->getdata_y(0, N_dim, 1) + grid->getdata_y(0, N_dim + 1, 1));
	grid->setdata_y(cV_0_N1_0, 0, N_dim + 1, 0);
	float cV_0_0_N1 = 0.33f * (grid->getdata_y(1, 0, N_dim + 1) + grid->getdata_y(0, 1, N_dim + 1) + grid->getdata_y(0, 0, N_dim));
	grid->setdata_y(cV_0_0_N1, 0, 0, N_dim + 1);
	float cV_0_N1_N1 = 0.33f * (grid->getdata_y(1, N_dim + 1, N_dim + 1) + grid->getdata_y(0, N_dim, N_dim + 1) + grid->getdata_y(0, N_dim + 1, N_dim));
	grid->setdata_y(cV_0_N1_N1, 0, N_dim + 1, N_dim + 1);
	float cV_N1_0_0 = 0.33f * (grid->getdata_y(N_dim, 0, 0) + grid->getdata_y(N_dim + 1, 1, 0) + grid->getdata_y(N_dim + 1, 0, 1));
	grid->setdata_y(cV_N1_0_0, N_dim + 1, 0, 0);
	float cV_N1_N1_0 = 0.33f * (grid->getdata_y(N_dim, N_dim + 1, 0) + grid->getdata_y(N_dim + 1, N_dim, 0) + grid->getdata_y(N_dim + 1, N_dim + 1, 1));
	grid->setdata_y(cV_N1_N1_0, N_dim + 1, N_dim + 1, 0);
	float cV_N1_0_N1 = 0.33f * (grid->getdata_y(N_dim, 0, N_dim + 1) + grid->getdata_y(N_dim + 1, 1, N_dim + 1) + grid->getdata_y(N_dim + 1, 0, N_dim));
	grid->setdata_y(cV_N1_0_N1, N_dim + 1, 0, N_dim + 1);
	float cV_N1_N1_N1 = 0.33f * (grid->getdata_y(N_dim, N_dim + 1, N_dim + 1) + grid->getdata_y(N_dim + 1, N_dim, N_dim + 1) + grid->getdata_y(N_dim + 1, N_dim + 1, N_dim));
	grid->setdata_y(cV_N1_N1_N1, N_dim + 1, N_dim + 1, N_dim + 1);
	// W (Z Component)
	float cW_0_0_0 = 0.33f * (grid->getdata_z(1, 0, 0) + grid->getdata_z(0, 1, 0) + grid->getdata_z(0, 0, 1));
	grid->setdata_z(cW_0_0_0, 0, 0, 0);
	float cW_0_N1_0 = 0.33f * (grid->getdata_z(1, N_dim + 1, 0) + grid->getdata_z(0, N_dim, 1) + grid->getdata_z(0, N_dim + 1, 1));
	grid->setdata_z(cW_0_N1_0, 0, N_dim + 1, 0);
	float cW_0_0_N1 = 0.33f * (grid->getdata_z(1, 0, N_dim + 1) + grid->getdata_z(0, 1, N_dim + 1) + grid->getdata_z(0, 0, N_dim));
	grid->setdata_z(cW_0_0_N1, 0, 0, N_dim + 1);
	float cW_0_N1_N1 = 0.33f * (grid->getdata_z(1, N_dim + 1, N_dim + 1) + grid->getdata_z(0, N_dim, N_dim + 1) + grid->getdata_z(0, N_dim + 1, N_dim));
	grid->setdata_z(cW_0_N1_N1, 0, N_dim + 1, N_dim + 1);
	float cW_N1_0_0 = 0.33f * (grid->getdata_z(N_dim, 0, 0) + grid->getdata_z(N_dim + 1, 1, 0) + grid->getdata_z(N_dim + 1, 0, 1));
	grid->setdata_z(cW_N1_0_0, N_dim + 1, 0, 0);
	float cW_N1_N1_0 = 0.33f * (grid->getdata_z(N_dim, N_dim + 1, 0) + grid->getdata_z(N_dim + 1, N_dim, 0) + grid->getdata_z(N_dim + 1, N_dim + 1, 1));
	grid->setdata_z(cW_N1_N1_0, N_dim + 1, N_dim + 1, 0);
	float cW_N1_0_N1 = 0.33f * (grid->getdata_z(N_dim, 0, N_dim + 1) + grid->getdata_z(N_dim + 1, 1, N_dim + 1) + grid->getdata_z(N_dim + 1, 0, N_dim));
	grid->setdata_z(cW_N1_0_N1, N_dim + 1, 0, N_dim + 1);
	float cW_N1_N1_N1 = 0.33f * (grid->getdata_z(N_dim, N_dim + 1, N_dim + 1) + grid->getdata_z(N_dim + 1, N_dim, N_dim + 1) + grid->getdata_z(N_dim + 1, N_dim + 1, N_dim));
	grid->setdata_z(cW_N1_N1_N1, N_dim + 1, N_dim + 1, N_dim + 1);

}


// ** SPHERE_BOUNDS_SCALAR-FIELD IMPLEMENTATION ** \\ - 

/* Set SphereBounds At Beginning Of Each Step:
   Calculate SDF Sphere Function from passed Radius, Offset for Each Cell in GridSpace. 
   Write into FluidSolver Spherebounds_SDF Signed Scalar Grid - */ 

void fluidsolver_3::sphere_bounds_set(float radius, float col_iso, const vec3<float> &offset)
{
	float h = 1.0f / N_dim; // Grid Spacing, Recoprical of One Dim Size (N). 

	//!MT #pragma omp parallel for num_threads(omp_get_max_threads())	

	for (int k = 1; k < N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// Implicit Sphere/sphere Function: x^2 + y^2 + z^2 - r . >= 0 <= Thresh == Surface. > Thresh = Exterior. < 0 = Interior. Cells. 
				vec3<float> cell_gridSpace(float(i * h) - offset.x, float(j * h) - offset.y, float(k*h) - offset.z); // Index to Grid Space 0-1N. 
				float sphere_func = ((cell_gridSpace.x * cell_gridSpace.x) + (cell_gridSpace.y * cell_gridSpace.y) + (cell_gridSpace.z * cell_gridSpace.z)) - radius; 

				// Write to  spherebounds_sdf grid for spherebounds_eval() to lookup. 
				spherebounds_sdf->setdata(0.0f, i, j, k); // Zero Out Prev SDF Grid Cell Values. 
				spherebounds_sdf->setdata(sphere_func, i, j, k); // Set New Cur Step SDF Cell Values. 
			}
		}
	}

}

// ** SPHERE_BOUNDS_EVAL - SCALAR-FIELD OVERRIDE ** \\ - 
// Eval SphereBounds - On Scalar Field. Also Set "Col" Viz Grid For Cells Inside Sphere SDF For Render col viz. 
// Unused (Surface and Exterior Condtions disabled for perf). 

void fluidsolver_3::sphere_bounds_eval(grid3_scalar<float> *grid, float col_iso)
{
	float h = 1.0f / N_dim; // Grid Spacing, Recoprical of One Dim Size (N). 
	
	//!MT #pragma omp parallel for num_threads(omp_get_max_threads())	

	for (int k = 1; k < N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// Lookup Sphere/sphere SDF At CurCell i,j. 
				float sphere_func = spherebounds_sdf->getdata(i, j, k);
				f3obj->col->setdata(0.0f, i, j, k); // Zero Out Prev Collison Cell Values. 

				/*
				// Surface -
				if (sphere_func >= 0.0f && sphere_func <= col_iso) // func-radis > 0.0f but <= col_iso Cell is on "surface".
				{
					// On Surface Operations.
				}
				*/

				// Inside - 
				if (sphere_func <= 0.0f)
				{
					// Inside Operations. 

					// Inside Operations -  Input Scalar Grid
					if (Parms.p_spherebounds_killint == true)
					{
						grid->setdata(0.0f, i, j, k); // Zero ANY Density inside Sphere Interior. 
					}

					// Inside Operations - Collision Scalar Grid
					// Set Collison Grid (HardCoded to pointed f3obj member) (Col Grid, probs should be Boolean Grid no? Evneutally MAC). 
					f3obj->col->setdata(1.0f, i, j, k);
				}

				/*
				if (sphere_func > col_iso) // func-radis > 0.0f+col_iso Cell is Outside.
				{
					// Outside Operations.
				}
				*/
			}
		}
	}

}


// ** SPHERE_BOUNDS_EVAL - VECTOR-FIELD OVERRIDE ** \\ - 
// Eval SphereBounds - On Vector Field. Also add Mouse Velocity to Vel Grid, Within SphereBounds. 
// (This Function does not set col grid)
// Unused (Surface and Exterior Condtions disabled for perf). 

void fluidsolver_3::sphere_bounds_eval(grid3_vector<vec3<float>> *grid, float col_iso)
{
	float h = 1.0f / N_dim; // Grid Spacing, Recoprical of One Dim Size (N). 

	//!MT #pragma omp parallel for num_threads(omp_get_max_threads())
	for (int k = 1; k < N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// Lookup Sphere/sphere SDF At CurCell i,j. 
				float sphere_func = spherebounds_sdf->getdata(i, j, k);

				/*
				// Surface
				if (sphere_func >= 0.0f && sphere_func <= col_iso) // func-radis > 0.0f but <= col_iso Cell is on "surface".
				{
					// Surface Operations -

					// Do Sphere Vel Reflection On Surface Cirlce/Sphere Cells -
					// Reflect from offset vector (Sphere Center) to CellPos Diriection, oppose to origin^ -

					//float input_spd = grid->getdata(i, j).length();
					//vec3 refl_dir = vec3(offset.x, offset.y) - vec3(float(i) * h, float(j) * h);
					//refl_dir.normalize();
					//refl_dir *= input_spd * 1.0f;

					// Override Vel With Reflect Dir Vel
					//grid->setdata(refl_dir, i, j);
				}
				*/

				// Inside - 
				if (sphere_func < 0.0f) // func-radius < 0.0f Cell is inside. 
				{
					// Inside Operations - 

					// Do Reflection On Interior Velocites Oppose to Zeroing Them? 
					/* NOT Used, as Causes Relfection on Interior Cells Causes Advection Issues at Surface.
					float input_spd = grid->getdata(i, j).length();
					vec3 refl_dir = vec3(float(i) * h, float(j) * h) - vec3(offset.x, offset.y) ;
					refl_dir.normalize();
					refl_dir *= input_spd * 1.0f;
					// Override Interior Col Vel With Reflect Dir Vel
					grid->setdata(refl_dir, i, j); // ISSUES WITH ADDED GRID VELOCITY ? */

					// Simply Zero Out Vector Grid Values In Interior Sphere Cells -   
					grid->setdata(vec3<float>(0.0f, 0.0f, 0.0f), i, j, k);
				}

				/*
				// Outside -
				if (sphere_func > col_iso) // func-radis > 0.0f+col_iso Cell is Outside.
				{
					// Outside Operations.
				}
				*/

				// !TODO 3D MOUSE VEL - 
				/* ADD MOUSE VELOCITY - Also On Surface & Interior Cells -
				Eventually add this as a Bool Param Option for Sphere_Bounds Vector, and expose Vel Multipler Param.*/
				if (sphere_func <= col_iso) // func-radius <= col_iso is Ethier Inside or on Surface. 
				{
					float mouse_velMult = 1.0f;

					// Get Stored Current Mouse Vel + to Current Grid Cell Velocity. 
					//vec2<float> cur_vel = grid->getdata(i, j);
					//grid->setdata(cur_vel + (mouse_vel *= mouse_velMult), i, j);

					//grid->setdata(grid->getdata(i, j, k) + (mouse_vel *= mouse_velMult), i, j, k);
				}
			}
		}
	}

}


/* ====================================================
	DIFFUSION
   ==================================================== */

// ** DIFFUSION-SCALAR-FIELD-LINSOLVE IMPLEMENTATION ** \\ - 

void fluidsolver_3::diffuse(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1, float diff, ushort iter)
{
	// Diffusion Value for "mat" A - 
	float a = dt * diff * powf(N_dim, 2.0f);

	// Scalar Field Diffusion - 
	for (int l = 0; l < iter; l++)
	{
		for (int k = 1; k <= N_dim; k++)
		{
			for (int j = 1; j <= N_dim; j++)
			{
				for (int i = 1; i <= N_dim; i++)
				{
					// Calc New Scalar Val via Diffusion Linear Solver. 
					float b_val = (grid_0->getdata(i, j, k) + a * 
						(grid_1->getdata(i - 1, j, k) + grid_1->getdata(i + 1, j, k) 
						+ grid_1->getdata(i, j - 1, k) + grid_1->getdata(i, j + 1, k)
						+ grid_1->getdata(i, j, k - 1) + grid_1->getdata(i, j, k + 1))
						) / (1 + 6 * a);

					// Set New Dens Val.
					grid_1->setdata(b_val, i, j, k);
				}
			}
		}

		// Call (Re-Inforce Boundary Condtions) Boundary Calc MFs on each Relaxation Iter - 
		edge_bounds(grid_1); // Generic Func Call, Pass Grid_1 (n+1). 

		#if dospherebound == 1
		sphere_bounds_eval(grid_1, spherebound_coliso);
		#endif
	}
}

// ** DIFFUSION-VECTOR-FIELD-LINSOLVE IMPLEMENTATION ** \\ - 

void fluidsolver_3::diffuse(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1, float diff, ushort iter)
{
	// Diffusion Value for "mat" A - 
	float a = dt * diff * powf(N_dim, 2);

	for (int l = 0; l < iter; l++)
	{
		for (int k = 1; k <= N_dim; k++)
		{
			for (int j = 1; j <= N_dim; j++)
			{
				for (int i = 1; i <= N_dim; i++)
				{
					// Per Scalar (float) Vector component diffusion - 

					// U (x)
					float U_val = (grid_0->getdata_x(i, j, k) + a * 
					(grid_1->getdata_x(i - 1, j, k) + grid_1->getdata_x(i + 1, j, k)
					+ grid_1->getdata_x(i, j - 1, k) + grid_1->getdata_x(i, j + 1, k)
					+ grid_1->getdata_x(i, j, k - 1) + grid_1->getdata_x(i, j, k + 1))
					 ) / (1 + 6 * a);

					// V (y)
					float V_val = (grid_0->getdata_y(i, j, k) + a *
					(grid_1->getdata_y(i - 1, j, k) + grid_1->getdata_y(i + 1, j, k)
					+ grid_1->getdata_y(i, j - 1, k) + grid_1->getdata_y(i, j + 1, k)
					+ grid_1->getdata_y(i, j, k - 1) + grid_1->getdata_y(i, j, k + 1))
					) / (1 + 6 * a);

					// W (z)
					float W_val = (grid_0->getdata_z(i, j, k) + a *
					(grid_1->getdata_z(i - 1, j, k) + grid_1->getdata_z(i + 1, j, k)
					+ grid_1->getdata_z(i, j - 1, k) + grid_1->getdata_z(i, j + 1, k)
					+ grid_1->getdata_z(i, j, k - 1) + grid_1->getdata_z(i, j, k + 1))
					) / (1 + 6 * a);

					vec3<float> new_vec(U_val, V_val, W_val);

					// Set New Vector Val.
					grid_1->setdata(new_vec, i, j, k);

				}
			}
		}
		// Call Boundary Calc MFs on each Relaxation Iter - 
		edge_bounds(grid_1);

		#if dospherebound == 1
		sphere_bounds_eval(grid_1, spherebound_coliso);
		#endif
	}
}

/* ====================================================
	ADVECTION
   ==================================================== */

// Semi Lagragain Advection - Overloads depend on Grid Type Passed. (Scalar,Vector Currently).
// Linear Backtrace along Velocity, Interoplate Neighbours at Backtraced Postion, to get our new advected value. 
// Single Step. MidPoint Coming soon. 

// ** ADVECT_SL(Semi-Lagragagin)_Scalar Overload ** \\ - 
// Assume ALL Advection is done using Main Vel Field. (Not using Input VelGrid Params). 
// BackTrace Done in Index Space
// Dt Scaled to Index Space Dim Size - 

void fluidsolver_3::advect_sl(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1)
{

	float dt0 = dt*N_dim; // Step Distance (Velocity Scaling) of DeltaTime * Grid Length, Hence GridSpace Dt Coeff (dt0). 

	// Density (Scalar Field Advection) - 

	//!MT#pragma omp parallel for num_threads(omp_get_max_threads()) 

	for (int k = 1; k <= N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				float u = f3obj->vel->getdata_x(i, j, k);
				float v = f3obj->vel->getdata_y(i, j, k);
				float w = f3obj->vel->getdata_z(i, j, k);

				// BackTrace Along U Component for X i 
				float x = i - dt0 * u;
				if (x < 0.5) x = 0.5;
				if (x > N_dim + 0.5) x = N_dim + 0.5;
				// Interp Indices i 
				int i0 = int(x); int i1 = i0 + 1;

				// BackTrace Along V Component for Y j 
				float y = j - dt0 * v;
				if (y < 0.5) y = 0.5;
				if (y > N_dim + 0.5) y = N_dim + 0.5;
				// Interp Indices j 
				int j0 = int(y); int j1 = j0 + 1;

				// BackTrace Along W Component for Z k
				float z = k - dt0 * w;
				if (z < 0.5) z = 0.5;
				if (z > N_dim + 0.5) z = N_dim + 0.5;
				// Interp Indices K
				int k0 = int(z); int k1 = j0 + 1;

				// Interoplation Coefficents - 
				float r1 = x - i0; float r0 = 1 - r1;
				float s1 = y - j0; float s0 = 1 - s1;
				float t1 = z - j0; float t0 = 1 - t1;

				/* //!TODO - Implement Cosine / LERP Switces -
				// Calc New Density from BackTrace Advection Neighbours Interoplation Of PrevFrame Grid - 
				float new_dens;

				if (Parms.p_InteroplationType == Parms.Interoplation_Linear)
				{
					// BiLinear Interoplation - 

					//new_dens = s0 * (t0*grid_0->getdata(i0, j0) + t1 * grid_0->getdata(i0, j1))
					//+ s1 * (t0 * grid_0->getdata(i1, j0) + t1 * grid_0->getdata(i1, j1));

					float lerpi_A = solver_utils::lerp(grid_0->getdata(i0, j0), grid_0->getdata(i0, j1), t1);
					float lerpi_B = solver_utils::lerp(grid_0->getdata(i1, j0), grid_0->getdata(i1, j1), t1);
					new_dens = solver_utils::lerp(lerpi_A, lerpi_B, s1);
				}
				else if (Parms.p_InteroplationType == Parms.Interoplation_Cosine)
				{
					// Cosine Interoplation - 

					float cosi_A = solver_utils::cosinterp(grid_0->getdata(i0, j0), grid_0->getdata(i0, j1), t1);
					float cosi_B = solver_utils::cosinterp(grid_0->getdata(i1, j0), grid_0->getdata(i1, j1), t1);
					float cosi_C = solver_utils::cosinterp(cosi_A, cosi_B, s1); // Interp Between Interp A,B.
					new_dens = cosi_C;
				}
				*/

				// TriLinear Interoplation Of Sampled Scalar Field Neighbours at BackTraced Postion - 

				float L_000_001_t = solver_utils::lerp(grid_0->getdata(i0, j0, k0), grid_0->getdata(i0, j0, k1), t1);
				float L_010_011_t = solver_utils::lerp(grid_0->getdata(i0, j1, k0), grid_0->getdata(i0, j1, k1), t1);
				float L_100_101_s = solver_utils::lerp(grid_0->getdata(i1, j0, k0), grid_0->getdata(i1, j0, k1), t1);
				float L_110_111_t = solver_utils::lerp(grid_0->getdata(i1, j1, k0), grid_0->getdata(i1, j1, k1), t1);
				float L_A = solver_utils::lerp(L_000_001_t, L_010_011_t, s1);
				float L_B = solver_utils::lerp(L_100_101_s, L_110_111_t, s1);
				float L_F = solver_utils::lerp(L_A, L_B, r1);
				float new_dens = L_F; 

				// Set New Cur Density to Current Frame Density Grid cell value - 
				grid_1->setdata(new_dens, i, j, k);

			}
		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso); 
	#endif
}

// ** ADVECT_SL(Semi-Lagragagin)_Vector Overload ** \\ - 
// Assume ALL Advection is done using Main Vel Field. (Not using Input VelGrid Params). 
// Backtrace Along Grid_1 (Cur), Sample Grid_0 (prev), Set New Grid_1 (Cur). 

void fluidsolver_3::advect_sl(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1)
{
	float dt0 = dt*N_dim; // Step Distance (Velocity Scaling) of DeltaTime * Grid Length, Hence GridSpace Dt Coeff (dt0). 

	//!MT#pragma omp parallel for num_threads(omp_get_max_threads())

	for (int k = 1; k <= N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				float u0 = f3obj->prev_vel->getdata_x(i, j, k);
				float u = f3obj->vel->getdata_x(i, j, k);
				float v0 = f3obj->prev_vel->getdata_y(i, j, k);
				float v = f3obj->vel->getdata_y(i, j, k);
				float w0 = f3obj->prev_vel->getdata_z(i, j, k);
				float w = f3obj->vel->getdata_z(i, j, k);

				// BackTrace Along U Component for X i 
				float x = i - dt0 * u;
				if (x < 0.5) x = 0.5;
				if (x > N_dim + 0.5) x = N_dim + 0.5;
				// Interp Indices i 
				int i0 = int(x); int i1 = i0 + 1;

				// BackTrace Along V Component for Y j 
				float y = j - dt0 * v;
				if (y < 0.5) y = 0.5;
				if (y > N_dim + 0.5) y = N_dim + 0.5;
				// Interp Indices j 
				int j0 = int(y); int j1 = j0 + 1;

				// BackTrace Along W Component for Z k
				float z = k - dt0 * w;
				if (z < 0.5) z = 0.5;
				if (z > N_dim + 0.5) z = N_dim + 0.5;
				// Interp Indices K
				int k0 = int(z); int k1 = j0 + 1;

				// Interoplation Coefficents - 
				float r1 = x - i0; float r0 = 1 - r1;
				float s1 = y - j0; float s0 = 1 - s1;
				float t1 = z - j0; float t0 = 1 - t1;

				/*
				// Interoplate Velocity U,V Components.  
				float new_u, new_v;
				if (Parms.p_InteroplationType == Parms.Interoplation_Linear)
				{

					//new_u = s0 * (t0*grid_0->getdata_x(i0, j0) + t1 * grid_0->getdata_x(i0, j1))
					//	+ s1 * (t0 * grid_0->getdata_x(i1, j0) + t1 * grid_0->getdata_x(i1, j1));
					//new_v = s0*(t0*grid_0->getdata_y(i0, j0) + t1 * grid_0->getdata_y(i0, j1))
					//	+ s1 * (t0 * grid_0->getdata_y(i1, j0) + t1 * grid_0->getdata_y(i1, j1));

					// Bilinear Interp (U/x Component)
					float U_lerpi_A = solver_utils::lerp(grid_0->getdata_x(i0, j0), grid_0->getdata_x(i0, j1), t1);
					float U_lerpi_B = solver_utils::lerp(grid_0->getdata_x(i1, j0), grid_0->getdata_x(i1, j1), t1);
					new_u = solver_utils::lerp(U_lerpi_A, U_lerpi_B, s1);

					// Bilinear Interp (V/y Component)
					float V_lerpi_A = solver_utils::lerp(grid_0->getdata_y(i0, j0), grid_0->getdata_y(i0, j1), t1);
					float V_lerpi_B = solver_utils::lerp(grid_0->getdata_y(i1, j0), grid_0->getdata_y(i1, j1), t1);
					new_v = solver_utils::lerp(V_lerpi_A, V_lerpi_B, s1);
				}
				else if (Parms.p_InteroplationType == Parms.Interoplation_Cosine)
				{
					// BiCosine Interp (U/x Component)
					float U_cosi_A = solver_utils::cosinterp(grid_0->getdata_x(i0, j0), grid_0->getdata_x(i0, j1), t1);
					float U_cosi_B = solver_utils::cosinterp(grid_0->getdata_x(i1, j0), grid_0->getdata_x(i1, j1), t1);
					new_u = solver_utils::cosinterp(U_cosi_A, U_cosi_B, s1);

					// BiCosine Interp (V/y Component)
					float V_cosi_A = solver_utils::cosinterp(grid_0->getdata_y(i0, j0), grid_0->getdata_y(i0, j1), t1);
					float V_cosi_B = solver_utils::cosinterp(grid_0->getdata_y(i1, j0), grid_0->getdata_y(i1, j1), t1);
					new_v = solver_utils::cosinterp(V_cosi_A, V_cosi_B, s1);

				}

				// Set New Cur vel vec3 to Velocity Grid cell value - 
				vec3<float> new_vec(new_u, new_v);
				grid_1->setdata(new_vec, i, j);
				*/

				// Interoplate Neighbours - for Velocity comp (U/x). 
				float U_000_001_t = solver_utils::lerp(grid_0->getdata_x(i0, j0, k0), grid_0->getdata_x(i0, j0, k1), t1);
				float U_010_011_t = solver_utils::lerp(grid_0->getdata_x(i0, j1, k0), grid_0->getdata_x(i0, j1, k1), t1);
				float U_100_101_s = solver_utils::lerp(grid_0->getdata_x(i1, j0, k0), grid_0->getdata_x(i1, j0, k1), t1);
				float U_110_111_t = solver_utils::lerp(grid_0->getdata_x(i1, j1, k0), grid_0->getdata_x(i1, j1, k1), t1);
				float U_A = solver_utils::lerp(U_000_001_t, U_010_011_t, s1);
				float U_B = solver_utils::lerp(U_100_101_s, U_110_111_t, s1);
				float U_F = solver_utils::lerp(U_A, U_B, r1);

				// Interoplate Neighbours - for Velocity comp (V/y). 
				float V_000_001_t = solver_utils::lerp(grid_0->getdata_y(i0, j0, k0), grid_0->getdata_y(i0, j0, k1), t1);
				float V_010_011_t = solver_utils::lerp(grid_0->getdata_y(i0, j1, k0), grid_0->getdata_y(i0, j1, k1), t1);
				float V_100_101_s = solver_utils::lerp(grid_0->getdata_y(i1, j0, k0), grid_0->getdata_y(i1, j0, k1), t1);
				float V_110_111_t = solver_utils::lerp(grid_0->getdata_y(i1, j1, k0), grid_0->getdata_y(i1, j1, k1), t1);
				float V_A = solver_utils::lerp(V_000_001_t, V_010_011_t, s1);
				float V_B = solver_utils::lerp(V_100_101_s, V_110_111_t, s1);
				float V_F = solver_utils::lerp(V_A, V_B, r1);

				// Interoplate Neighbours - for Velocity comp (W/z). 
				float W_000_001_t = solver_utils::lerp(grid_0->getdata_z(i0, j0, k0), grid_0->getdata_z(i0, j0, k1), t1);
				float W_010_011_t = solver_utils::lerp(grid_0->getdata_z(i0, j1, k0), grid_0->getdata_z(i0, j1, k1), t1);
				float W_100_101_s = solver_utils::lerp(grid_0->getdata_z(i1, j0, k0), grid_0->getdata_z(i1, j0, k1), t1);
				float W_110_111_t = solver_utils::lerp(grid_0->getdata_z(i1, j1, k0), grid_0->getdata_z(i1, j1, k1), t1);
				float W_A = solver_utils::lerp(W_000_001_t, W_010_011_t, s1);
				float W_B = solver_utils::lerp(W_100_101_s, W_110_111_t, s1);
				float W_F = solver_utils::lerp(W_A, W_B, r1);

				// Set New Cur Velocity (grid_1) - 
				grid_1->setdata(vec3<float>(U_F, V_F, W_F), i, j, k);
			}
		}
	}

	// Call Boundary Condtions Post Advection (Density)- 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}

/* MID POINT ADVECTION (RK2) WIP - 
BackTrace to MidPoint, Sample (interoplate) Velocity at MidPoint, Then BackTrace from Cell using MidPoint Vel, 
to sample (interoplate) Final BackTraced Advection Quanitiy
--
Advection in Grid Index Space. Dt Scaled to N_Dim size.
Velocity Components split for Advection Sampling and Interoplation.
*/

// ** ADVECT_SL_RK2(Semi-Lagragagin_MidPoint)_Scalar Overload ** \\ - 
void fluidsolver_3::advect_sl_mp(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1)
{
	float dt0 = dt * N_dim; // Scale DeltaTime to Grid Dimension Size. 

	// Density (Scalar Field Advection) - 

	//!MT#pragma omp parallel for num_threads(omp_get_max_threads()) 

	for (int k = 1; k <= N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// Vel at Cur Cell Postion. 
				float u_P = f3obj->vel->getdata_x(i, j, k);
				float v_P = f3obj->vel->getdata_y(i, j, k);
				float w_P = f3obj->vel->getdata_z(i, j, k);

				// BackTrace Along Negative CurCell Vel - XG - dt0 * u(CurCell)
				// XG -> Midpoint = XG - dt0 * u(XG)
				// BackTrace U,V,W Velocity Components to Midpoint.
				float xxm = i - (0.5 * dt0) * u_P; 
				float yym = j - (0.5 * dt0) * v_P; 
				float zzm = k - (0.5 * dt0) * w_P; 
				if (xxm < 0.5) xxm = 0.5; if (xxm > N_dim + 0.5) xxm = N_dim + 0.5;
				if (yym < 0.5) yym = 0.5; if (yym > N_dim + 0.5) yym = N_dim + 0.5;
				if (zzm < 0.5) zzm = 0.5; if (zzm > N_dim + 0.5) zzm = N_dim + 0.5;

				// MidPoint - Mid Indices 
				int i_mid = int(xxm); int i_mid_1 = i_mid + 1; 
				int j_mid = int(yym); int j_mid_1 = j_mid + 1;
				int k_mid = int(zzm); int k_mid_1 = k_mid + 1;

				// MidPoint - Mid Interp Coefficents - 
				float rm1 = xxm - i_mid; float rm = 1 - rm1;
				float sm1 = yym - j_mid; float sm = 1 - sm1;
				float tm1 = yym - j_mid; float tm = 1 - tm1;

				// Get Mid Point Velocity (Trilinear Interoplation of Velocity Components u,v,w At Midpoint Postion) - 
				// Interoplate Neighbours - for Velocity comp (U/x). 
				float Um_000_001_t = solver_utils::lerp(f3obj->vel->getdata_x(i_mid, j_mid, k_mid), f3obj->vel->getdata_x(i_mid, j_mid, k_mid_1), tm1);
				float Um_010_011_t = solver_utils::lerp(f3obj->vel->getdata_x(i_mid, j_mid_1, k_mid), f3obj->vel->getdata_x(i_mid, j_mid_1, k_mid), tm1);
				float Um_100_101_s = solver_utils::lerp(f3obj->vel->getdata_x(i_mid_1, j_mid, k_mid), f3obj->vel->getdata_x(i_mid_1, j_mid, k_mid_1), tm1);
				float Um_110_111_t = solver_utils::lerp(f3obj->vel->getdata_x(i_mid_1, j_mid_1, k_mid), f3obj->vel->getdata_x(i_mid_1, j_mid_1, k_mid_1), tm1);
				float Um_A = solver_utils::lerp(Um_000_001_t, Um_010_011_t, sm1);
				float Um_B = solver_utils::lerp(Um_100_101_s, Um_110_111_t, sm1);
				// Interoplate Neighbours - for Velocity comp (V/y). 
				float Vm_000_001_t = solver_utils::lerp(f3obj->vel->getdata_y(i_mid, j_mid, k_mid), f3obj->vel->getdata_y(i_mid, j_mid, k_mid_1), tm1);
				float Vm_010_011_t = solver_utils::lerp(f3obj->vel->getdata_y(i_mid, j_mid_1, k_mid), f3obj->vel->getdata_y(i_mid, j_mid_1, k_mid_1), tm1);
				float Vm_100_101_s = solver_utils::lerp(f3obj->vel->getdata_y(i_mid_1, j_mid, k_mid), f3obj->vel->getdata_y(i_mid_1, j_mid, k_mid_1), tm1);
				float Vm_110_111_t = solver_utils::lerp(f3obj->vel->getdata_y(i_mid_1, j_mid_1, k_mid), f3obj->vel->getdata_y(i_mid_1, j_mid_1, k_mid_1), tm1);
				float Vm_A = solver_utils::lerp(Vm_000_001_t, Vm_010_011_t, sm1);
				float Vm_B = solver_utils::lerp(Vm_100_101_s, Vm_110_111_t, sm1);
				// Interoplate Neighbours - for Velocity comp (W/z). 
				float Wm_000_001_t = solver_utils::lerp(f3obj->vel->getdata_z(i_mid, j_mid, k_mid), f3obj->vel->getdata_z(i_mid, j_mid, k_mid_1), tm1);
				float Wm_010_011_t = solver_utils::lerp(f3obj->vel->getdata_z(i_mid, j_mid_1, k_mid), f3obj->vel->getdata_z(i_mid, j_mid_1, k_mid_1), tm1);
				float Wm_100_101_s = solver_utils::lerp(f3obj->vel->getdata_z(i_mid_1, j_mid, k_mid), f3obj->vel->getdata_z(i_mid_1, j_mid, k_mid_1), tm1);
				float Wm_110_111_t = solver_utils::lerp(f3obj->vel->getdata_z(i_mid_1, j_mid_1, k_mid), f3obj->vel->getdata_z(i_mid_1, j_mid_1, k_mid_1), tm1);
				float Wm_A = solver_utils::lerp(Wm_000_001_t, Wm_010_011_t, sm1);
				float Wm_B = solver_utils::lerp(Wm_100_101_s, Wm_110_111_t, sm1);
				
				// Interoplated Mid Point Velocity -
				float u_mid = solver_utils::lerp(Um_A, Um_B, rm1);
				float v_mid = solver_utils::lerp(Vm_A, Vm_B, rm1);
				float w_mid = solver_utils::lerp(Wm_A, Wm_B, rm1);

				// BackTrace Along Negative Midpoint Vel - XG - dt0 * u(midpoint)
				float xxp = i - dt0 * u_mid; 
				float yyp = j - dt0 * v_mid;
				float zzp = k - dt0 * v_mid;
				if (xxp < 0.5) xxp = 0.5; if (xxp > N_dim + 0.5) xxp = N_dim + 0.5;
				if (yyp < 0.5) yyp = 0.5; if (yyp > N_dim + 0.5) yyp = N_dim + 0.5;
				if (zzp < 0.5) zzp = 0.5; if (zzp > N_dim + 0.5) zzp = N_dim + 0.5;

				// MidPoint - Mid Indices 
				int i_f = int(xxp); int i_f_1 = i_f + 1;
				int j_f = int(yyp); int j_f_1 = j_f + 1;
				int k_f = int(zzp); int k_f_1 = k_f + 1;

				// MidPoint - Mid Interp Coefficents  
				float rf1 = xxm - i_mid; float rf = 1 - rf1;
				float sf1 = yym - j_mid; float sf = 1 - sf1;
				float tf1 = yym - j_mid; float tf = 1 - tf1;

				// Trilinearly Sample Scalar Field at backtraced postion (via MidPoint vel) -
				float L_000_001_t = solver_utils::lerp(grid_0->getdata(i_f, j_f, k_f), grid_0->getdata(i_f, j_f, k_f_1), tf1);
				float L_010_011_t = solver_utils::lerp(grid_0->getdata(i_f, j_f_1, k_f), grid_0->getdata(i_f, j_f_1, k_f_1), tf1);
				float L_100_101_s = solver_utils::lerp(grid_0->getdata(i_f_1, j_f, k_f), grid_0->getdata(i_f_1, j_f, k_f_1), tf1);
				float L_110_111_t = solver_utils::lerp(grid_0->getdata(i_f_1, j_f_1, k_f), grid_0->getdata(i_f_1, j_f_1, k_f_1), tf1);
				float L_A = solver_utils::lerp(L_000_001_t, L_010_011_t, sf1);
				float L_B = solver_utils::lerp(L_100_101_s, L_110_111_t, sf1);
				float L_F = solver_utils::lerp(L_A, L_B, rf1);

				// Set Backtraced Scalar Value to Cur Cell in grid_1 -
				grid_1->setdata(L_F, i, j, k);

			}
		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}

// ** ADVECT_SL_RK2(Semi-Lagragagin_MidPoint)_Vector Overload ** \\ - 
void fluidsolver_3::advect_sl_mp(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1)
{
	float dt0 = dt * N_dim; // Scale DeltaTime to Grid Dimension Size. 

	//!MT#pragma omp parallel for num_threads(omp_get_max_threads()) 

	for (int k = 1; k <= N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// Vel at Cur Cell Postion. 
				float u_P = f3obj->vel->getdata_x(i, j, k);
				float v_P = f3obj->vel->getdata_y(i, j, k);
				float w_P = f3obj->vel->getdata_z(i, j, k);

				// BackTrace Along Negative CurCell Vel - XG - dt0 * u(CurCell)
				// XG -> Midpoint = XG - dt0 * u(XG)
				// BackTrace U,V,W Velocity Components to Midpoint.
				float xxm = i - (0.5 * dt0) * u_P;
				float yym = j - (0.5 * dt0) * v_P;
				float zzm = k - (0.5 * dt0) * w_P;
				if (xxm < 0.5) xxm = 0.5; if (xxm > N_dim + 0.5) xxm = N_dim + 0.5;
				if (yym < 0.5) yym = 0.5; if (yym > N_dim + 0.5) yym = N_dim + 0.5;
				if (zzm < 0.5) zzm = 0.5; if (zzm > N_dim + 0.5) zzm = N_dim + 0.5;

				// MidPoint - Mid Indices 
				int i_mid = int(xxm); int i_mid_1 = i_mid + 1;
				int j_mid = int(yym); int j_mid_1 = j_mid + 1;
				int k_mid = int(zzm); int k_mid_1 = k_mid + 1;

				// MidPoint - Mid Interp Coefficents - 
				float rm1 = xxm - i_mid; float rm = 1 - rm1;
				float sm1 = yym - j_mid; float sm = 1 - sm1;
				float tm1 = yym - j_mid; float tm = 1 - tm1;

				// Get Mid Point Velocity (Trilinear Interoplation of Velocity Components u,v,w At Midpoint Postion (m)) - 
				// Interoplate Neighbours - for Velocity comp (U/x). 
				float Um_000_001_t = solver_utils::lerp(f3obj->vel->getdata_x(i_mid, j_mid, k_mid), f3obj->vel->getdata_x(i_mid, j_mid, k_mid_1), tm1);
				float Um_010_011_t = solver_utils::lerp(f3obj->vel->getdata_x(i_mid, j_mid_1, k_mid), f3obj->vel->getdata_x(i_mid, j_mid_1, k_mid), tm1);
				float Um_100_101_s = solver_utils::lerp(f3obj->vel->getdata_x(i_mid_1, j_mid, k_mid), f3obj->vel->getdata_x(i_mid_1, j_mid, k_mid_1), tm1);
				float Um_110_111_t = solver_utils::lerp(f3obj->vel->getdata_x(i_mid_1, j_mid_1, k_mid), f3obj->vel->getdata_x(i_mid_1, j_mid_1, k_mid_1), tm1);
				float Um_A = solver_utils::lerp(Um_000_001_t, Um_010_011_t, sm1);
				float Um_B = solver_utils::lerp(Um_100_101_s, Um_110_111_t, sm1);
				// Interoplate Neighbours - for Velocity comp (V/y). 
				float Vm_000_001_t = solver_utils::lerp(f3obj->vel->getdata_y(i_mid, j_mid, k_mid), f3obj->vel->getdata_y(i_mid, j_mid, k_mid_1), tm1);
				float Vm_010_011_t = solver_utils::lerp(f3obj->vel->getdata_y(i_mid, j_mid_1, k_mid), f3obj->vel->getdata_y(i_mid, j_mid_1, k_mid_1), tm1);
				float Vm_100_101_s = solver_utils::lerp(f3obj->vel->getdata_y(i_mid_1, j_mid, k_mid), f3obj->vel->getdata_y(i_mid_1, j_mid, k_mid_1), tm1);
				float Vm_110_111_t = solver_utils::lerp(f3obj->vel->getdata_y(i_mid_1, j_mid_1, k_mid), f3obj->vel->getdata_y(i_mid_1, j_mid_1, k_mid_1), tm1);
				float Vm_A = solver_utils::lerp(Vm_000_001_t, Vm_010_011_t, sm1);
				float Vm_B = solver_utils::lerp(Vm_100_101_s, Vm_110_111_t, sm1);
				// Interoplate Neighbours - for Velocity comp (W/z). 
				float Wm_000_001_t = solver_utils::lerp(f3obj->vel->getdata_z(i_mid, j_mid, k_mid), f3obj->vel->getdata_z(i_mid, j_mid, k_mid_1), tm1);
				float Wm_010_011_t = solver_utils::lerp(f3obj->vel->getdata_z(i_mid, j_mid_1, k_mid), f3obj->vel->getdata_z(i_mid, j_mid_1, k_mid_1), tm1);
				float Wm_100_101_s = solver_utils::lerp(f3obj->vel->getdata_z(i_mid_1, j_mid, k_mid), f3obj->vel->getdata_z(i_mid_1, j_mid, k_mid_1), tm1);
				float Wm_110_111_t = solver_utils::lerp(f3obj->vel->getdata_z(i_mid_1, j_mid_1, k_mid), f3obj->vel->getdata_z(i_mid_1, j_mid_1, k_mid_1), tm1);
				float Wm_A = solver_utils::lerp(Wm_000_001_t, Wm_010_011_t, sm1);
				float Wm_B = solver_utils::lerp(Wm_100_101_s, Wm_110_111_t, sm1);

				// Interoplated Mid Point Velocity -
				float u_mid = solver_utils::lerp(Um_A, Um_B, rm1);
				float v_mid = solver_utils::lerp(Vm_A, Vm_B, rm1);
				float w_mid = solver_utils::lerp(Wm_A, Wm_B, rm1);

				// BackTrace Along Negative Midpoint Vel - XG - dt0 * u(midpoint) 
				float xxp = i - dt0 * u_mid;
				float yyp = j - dt0 * v_mid;
				float zzp = k - dt0 * v_mid;
				if (xxp < 0.5) xxp = 0.5; if (xxp > N_dim + 0.5) xxp = N_dim + 0.5;
				if (yyp < 0.5) yyp = 0.5; if (yyp > N_dim + 0.5) yyp = N_dim + 0.5;
				if (zzp < 0.5) zzp = 0.5; if (zzp > N_dim + 0.5) zzp = N_dim + 0.5;

				// MidPoint - Mid Indices 
				int i_f = int(xxp); int i_f_1 = i_f + 1;
				int j_f = int(yyp); int j_f_1 = j_f + 1;
				int k_f = int(zzp); int k_f_1 = k_f + 1;

				// MidPoint - Mid Interp Coefficents  
				float rf1 = xxm - i_mid; float rf = 1 - rf1;
				float sf1 = yym - j_mid; float sf = 1 - sf1;
				float tf1 = yym - j_mid; float tf = 1 - tf1;

				// Trilinearly Sample Velocity Components (u,v,w) at Back Traced Postion (f) 
				// Interoplate Neighbours - for Velocity comp (U/x). 
				float Uf_000_001_t = solver_utils::lerp(grid_0->getdata_x(i_f, j_f, k_f), grid_0->getdata_x(i_f, j_f, k_f_1), tf1);
				float Uf_010_011_t = solver_utils::lerp(grid_0->getdata_x(i_f, j_f_1, k_f), grid_0->getdata_x(i_f, j_f_1, k_f), tf1);
				float Uf_100_101_s = solver_utils::lerp(grid_0->getdata_x(i_f_1, j_f, k_f), grid_0->getdata_x(i_f_1, j_f, k_f_1), tf1);
				float Uf_110_111_t = solver_utils::lerp(grid_0->getdata_x(i_f_1, j_f_1, k_f), grid_0->getdata_x(i_f_1, j_f_1, k_f_1), tf1);
				float Uf_A = solver_utils::lerp(Um_000_001_t, Um_010_011_t, sm1);
				float Uf_B = solver_utils::lerp(Um_100_101_s, Um_110_111_t, sm1);
				// Interoplate Neighbours - for Velocity comp (V/y). 
				float Vf_000_001_t = solver_utils::lerp(grid_0->getdata_y(i_f, j_f, k_f), grid_0->getdata_y(i_f, j_f, k_f_1), tf1);
				float Vf_010_011_t = solver_utils::lerp(grid_0->getdata_y(i_f, j_f_1, k_f), grid_0->getdata_y(i_f, j_f_1, k_f_1), tf1);
				float Vf_100_101_s = solver_utils::lerp(grid_0->getdata_y(i_f_1, j_f, k_f), grid_0->getdata_y(i_f_1, j_f, k_f_1), tf1);
				float Vf_110_111_t = solver_utils::lerp(grid_0->getdata_y(i_f_1, j_f_1, k_f), grid_0->getdata_y(i_f_1, j_f_1, k_f_1), tf1);
				float Vf_A = solver_utils::lerp(Vm_000_001_t, Vm_010_011_t, sm1);
				float Vf_B = solver_utils::lerp(Vm_100_101_s, Vm_110_111_t, sm1);
				// Interoplate Neighbours - for Velocity comp (W/z). 
				float Wf_000_001_t = solver_utils::lerp(grid_0->getdata_z(i_f, j_f, k_f), grid_0->getdata_z(i_f, j_f, k_f_1), tf1);
				float Wf_010_011_t = solver_utils::lerp(grid_0->getdata_z(i_f, j_f_1, k_f), grid_0->getdata_z(i_f, j_f_1, k_f_1), tf1);
				float Wf_100_101_s = solver_utils::lerp(grid_0->getdata_z(i_f_1, j_f, k_f), grid_0->getdata_z(i_f_1, j_f, k_f_1), tf1);
				float Wf_110_111_t = solver_utils::lerp(grid_0->getdata_z(i_f_1, j_f_1, k_f), grid_0->getdata_z(i_f_1, j_f_1, k_f_1), tf1);
				float Wf_A = solver_utils::lerp(Wm_000_001_t, Wm_010_011_t, sf1);
				float Wf_B = solver_utils::lerp(Wm_100_101_s, Wm_110_111_t, sf1);

				// Final Sampled Interoplated Back Traced Point Velocity -
				float U_f = solver_utils::lerp(Uf_A, Uf_B, rf1);
				float V_f = solver_utils::lerp(Vf_A, Vf_B, rf1);
				float W_f = solver_utils::lerp(Wf_A, Wf_B, rf1);

				// Set Final Advected Velocity to cur grid - 
				grid_1->setdata(vec3<float>(U_f, V_f, W_f), i, j, k);
			}
		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval (grid_1, spherebound_coliso);
	#endif
}


/* ====================================================
				FORCES AND MISC - 
	==================================================== */
/* !TODO LATER

// Loop Through Vel Grid, and Add passed force to velocity * dt. 
void fluidsolver_3::vel_force(vec3<float> ff, float dtt)
{
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			//vec3<float> temp (ff.x, ff.y); // UN NEGATE AXIS for coord fix/dbg ?
			//f3obj->integrate_force(temp, dtt, i, j);   // Call add_force MF in fluid_obj.
		}
	}
}

// Pass Custom Force Callback (Lambda) to run over each vector grid cell and integrate as force. 
void fluidsolver_3::custom_force(grid3_vector<vec3<float>> *grid, std::function <vec3<float>(vec3<int> idx)> &force)
{
	#pragma omp parallel for 
	for (int j = 1; j <= N_dim; j++)
	{
		#pragma omp parallel for 
		for (int i = 1; i <= N_dim; i++)
		{
			// Call Callback - Get Force Vector. 
			vec3<float> force_vec = force(vec3<int>(i, j));
			// Integrate Force -
			f3obj->integrate_force(force_vec, dt, i, j);
		}
	}
}
*/

/*	====================================================
	VORTICITY CONFINEMENT
==================================================== */

// Compute and Integrate Vorticity Confinement Force, from/to Velocity Field.

// NOTES/Learnt - 
// CURL in 2D is a Scalar, not a vector. As theres no Cross Product for single XY Plane Vel.

/*
void fluidsolver_3::vorticty_confine(float strength)
{
	// VCForce = strength * cell size * (N cross W) ... Intergrate to Vel Field over Dt. 
	float dx = 1.0f / N_dim; // (h / cell size)

	// Allocate Vorticity Grids (Per Call)
	vort = new grid3_vector<vec3<float>>(X_dim, Y_dim, edge_size, 6, 1.0f); // Vorticity Field Grid
	deltavort = new grid3_vector<vec3<float>>(X_dim, Y_dim, edge_size, 7, 1.0f); // Vorticity Gradient Field Grid. 

	// Cap f3obj ptr by ref [&]. 
	std::function<float(int, int)> curl2d = [&](int i , int j) -> float 
	{
		return float(f3obj->vel->getdata_x(i, j + 1) - f3obj->vel->getdata_x(i, j - 1)) 
		+ (f3obj->vel->getdata_y(i + 1, j) - f3obj->vel->getdata_x(i, j - 1));
	};

	/*
		// Curl From Velocity (3D) -  (Curl From Velocity Field Function I wrote in VEX in Houdini DOPs)

		// a(j - l) 
		vector zy = (volumeindex(0, "vel",set(i,j+1,k)) - volumeindexv(0, "vel", set(i,j-1,k)) / h) - (volumeindexv(0, "vel", set(i,j,k+1)) - volumeindexv(0, "vel", set(i,j,k-1)) / h); 
		// a(k - i)
		vector xz = (volumeindexv(0, "vel", set(i,j,k+1)) - volumeindexv(0, "vel", set(i,j,k-1)) / h) - (volumeindexv(0, "vel", set(i+1,j,k)) - volumeindexv(0, "vel", set(i-1,j,k)) / h); 
		// a(i - j)
		vector yx = (volumeindexv(0, "vel", set(i+1,j,k)) - volumeindexv(0, "vel", set(i-1,j,k)) / h) - (volumeindexv(0, "vel", set(i,j+1,k)) - volumeindexv(0, "vel", set(i,j-1,k)) / h); 
		vector curl_v = zy + xz + yx;
		v@curl = curl_v;

		// Curl From Velocity (2D) - Ofc 2D we only have a single plane XY, so just need ^ for XY Partial Derivative Subtraction Terms (Facing Vector in 2D) - 
		// a(x - y)
		vector xy = (volumeindexv(0, "vel", set(i+1,j,k)) - volumeindexv(0, "vel", set(i-1,j,k)) / h) - (volumeindexv(0, "vel", set(i,j+1,k)) - volumeindexv(0, "vel", set(i,j-1,k)) / h);
		v@curl = xy; 

		I'm Calling Curl (w) Vorticity Here ... Bridson dentotes (w') as Voriticty/Curl, because he uses w as vel axis of Z, as u,v,w reprenseting velocity components for 3Dimensions. 

		Apprently 2D Length of Curl Directly Calculating is - 
		curl_length = (f3obj->vel->getdata_x(i, j + 1) - f3obj->vel->getdata_x(i, j - 1)) - (f3obj->vel->getdata_y(i + 1, j) - f3obj->vel->getdata_x(i, j - 1));

	*/

	/*
	// Caluclate Vorticity (Curl w) From Velocity Grid, Write to Vorticity Grid. // Exclude Ghost/Edge Cells. 
	for (int j = 1; j < N_dim; j++)
	{
		for (int i = 1; i < N_dim; i++)
		{
		    // Central Diffrence To Calc Curl from XY Velocity - / i +/- (j) constant, then viceversa. Divide by h (dx) *2 Because Over 2 Cells +/- Distance h. 
		//	vec3 w_xy = (f3obj->vel->getdata(i + 1, j) - f3obj->vel->getdata(i - 1, j) /= dx*2) - (f3obj->vel->getdata(i, j + 1) - f3obj->vel->getdata(i, j - 1) /= dx*2);
			vec3 w_xy = (f3obj->vel->getdata(i + 1, j) - f3obj->vel->getdata(i - 1, j)) - (f3obj->vel->getdata(i, j + 1) - f3obj->vel->getdata(i, j - 1));
			w_xy /= dx * 2; 
			vort->setdata(w_xy, i, j); // Set Resulting Curl/Vorticity (w) Vector Per Voxel. 
		}
	}
	*/

	/* 
	Calculate Vorticity Gradient from ||Vorticity|| HAS to be done sepereate Loop, because All Values (+/- i,J) 
	need to be calculated and stored to vort/curl grid, to then calc gradient. 
	


	for (int j = 2; j < N_dim-1; j++)
	{
		for (int i = 2; i < N_dim-1; i++)
		{
			/*
			float curl_l = (vort->getdata(i,j)).length(); // Length of Curl/Vorticity Value to Calc Gradient Of. 

			// Calc N ie Normalized Gradient/Delta Vort of Curl Length Central Diffrence of Scalar Curl Length Per Partial D, Dimension X,Y. 
			float xx = ((vort->getdata(i+1, j)).length() - (vort->getdata(i-1, j)).length()) / dx * 2;
			float yy = ((vort->getdata(i + 1, j)).length() - (vort->getdata(i - 1, j)).length()) / dx * 2;
			vec3 grad(xx, yy); 
			//grad.normalize(); // Normalize Gradient N

			// Resulting Vorticity Confinement Vector W X N
			vec3 vc_vec = vort->getdata(i, j).cross(grad); 
			vec3 vc_force = vc_vec *= (strength * dx);

			// VCForce = Strength * dx * (W X N)
			// Vel += VCForce * dt; // Integrate VC Force to vel. 

			// Debuggin currently so messing with what I Intergrate to vel field... 

			// Get Current Cell Velocity From VelField - 
			vec3 vel = f3obj->vel->getdata(i, j);
			// vel = vel + (vc_force *= dt) ; // Integrate vc force. 
			vec3 vcurl = vort->getdata(i, j);
			//vcurl.normalize();
			vel = vel + (vcurl *= 0.002f); // Integrate vc force. 
			f3obj->vel->setdata(vel, i, j); // Set to new Vel Value. 
		//	f3obj->vel->setdata(vel + vort->getdata(i,j), i, j); // Set to new Vel Value. 
			

			//float xx = std::abs(curl2d(i, j+1)) - std::abs(curl2d(i, j-1));
			float xx = curl2d(i, j - 1) - curl2d(i, j + 1);
			float yy = curl2d(i+1, j) - curl2d(i-1, j);
			vec3 dir(xx, yy);

		//	float f = (10.0f / dir.length()) + 1e-05f;
		//	dir *= f; 
			vec3 old_vel = f3obj->vel->getdata(i, j);
		//	float dtc = dt * curl2d(i, j);
		//	vec3 vc = dir *= dtc; 
		//	vec3 new_vel = old_vel + vc; 
			f3obj->vel->setdata(old_vel + ((dir += 1e-05f) *= 5.0f * dt), i, j);


		}
	}



	// Set Force Value (Should be done in ^ Loop and intergrated within, oppose to writing to third final
	// vorticity confinement grid). 
	// Intergrate Force to Velocity ^^

	// Delete Vorticity Grids - 
	delete vort; vort = nullptr;
	delete deltavort; deltavort = nullptr; 

	// In Future Might Make Vorticity Grid a FLuidObj Grid member, so Its lifespan/scope is same as fluidobj, 
	// and can be passed to GPU as texture for debug drawing/viz like Vel Field. 

}
*/

/*! TODO Later 
// Vorticity Confinement On The Fly (No Writing to Vorticity/Curl Grids, Curl, Curl Grad and Resulting VC Force
// calculated all within one loop, via Curl Lambda to calc  Curl (w) and Curl Gradient (w'). 
void fluidsolver_3::vorticty_confine_otf(float strength)
{
	// Calculate 2D Velocity Curl Scalar - Without Central Diff Offsets, for Calling within Curl gradient CentralDiff calc lines. 
	std::function<float(int, int)> curl2d = [&](int i, int j) -> float // Cap f3obj ptr by ref [&]. 
	{
		// i j offsets will be set on call. 
		return float(f3obj->vel->getdata_x(i, j) - f3obj->vel->getdata_x(i, j))
			+ (f3obj->vel->getdata_y(i, j) - f3obj->vel->getdata_x(i, j));
		// Should be / dx*2 ie central diff over h. 
	};

	// Curl With ij central diff offsets, for calling when eval curl only. 
	std::function<float(int, int)> curl2d_wo = [&](int i, int j) -> float // Cap f3obj ptr by ref [&]. 
	{
		// i j set with offset.
		return float(f3obj->vel->getdata_x(i, j+1) - f3obj->vel->getdata_x(i, j-1))
			+ (f3obj->vel->getdata_y(i+1, j) - f3obj->vel->getdata_x(i-1, j));
		// Should be / dx*2 ie central diff over h. 
	};

	for (int j = 1; j < N_dim - 1; j++)
	{
		for (int i = 1; i < N_dim - 1; i++)
		{
			// w' from w Calc both curl and gradient curl per axis using central diff. 
			float xx = curl2d(i, j-1) - curl2d(i, j+1);
			float yy = curl2d(i+1, j) - curl2d(i-1, j);
			vec3<float> dir(xx, yy); // resulting w' grad vector. 
			f3obj->vc->setdata(dir, i, j); // Pass To VC Grid for viz. Pre Normalization (will be clampped vals on GPU).
			dir += 1e-05; // Prevent /0 errors on normalization.
			dir.normalize();

			// Integrate VC to new vel value.
			vec3<float> old_vel = f3obj->vel->getdata(i, j);
			dir *= strength; // Mult By Strength Scalar Coeff. 
			//f3obj->vel->setdata(old_vel + (dir *= (dt * curl2d(i, j))), i, j);
			// v + strength*(Curl (w) and GradCurl (w')*dt) = v1 (v with VC). 
			f3obj->vel->setdata(old_vel + ((dir *= dt) *= curl2d_wo(i,j)), i, j);
		}
	}
}
// End of Vorticity_Confine_otf Member Func Implementation.

void fluidsolver_3::vorticity_confine_B(float strength)
{
	// Alloc Vorticity Grid - 
	vort = new grid3_scalar<float>(x_s, y_s, e_s, 5, 1.0f);

	float dx = 1.0f / N_dim; // cellsize dx (h). 

	// Calculate 2D Velocity Curl Scalar - At Given Cell Indices- 
	std::function<float(int, int)> curl2d = [&](int i, int j) -> float // Cap f3obj ptr by ref [&]. 
	{
		return float(
			(f3obj->vel->getdata_x(i, j+1) - f3obj->vel->getdata_x(i, j-1) / (2 * dx))
			- 
			(f3obj->vel->getdata_y(i+1, j) - f3obj->vel->getdata_y(i-1, j) / (2 * dx))
			);
	};

	// Calculate 2D Vorticity Gradient Vector - At Given Cell Indices- - 
	std::function<vec3<float>(int, int)> grad2d = [&](int i, int j) -> vec3<float>
	{
		float xx = (vort->getdata(i + 1, j) - vort->getdata(i - 1, j)) / (2 * dx); // 2.0f * dx;
		float yy = (vort->getdata(i, j + 1) - vort->getdata(i, j - 1)) / (2 * dx); // 2.0f * dx;
		vec3<float> grad(xx, yy);
		grad += 1e-05; // Add Small Constant, to prevent divide by 0 error on Grad Vector Normalization. 
		grad.normalize(); 
		return grad; 
	};

	// Calculate Curl and write to vort grid- 
	for (int j = 1; j < N_dim; j++)
	{
		for (int i = 1; i < N_dim; i++)
		{
			vort->setdata(curl2d(i, j), i, j);
		}
	}

	// Calc and Apply Vortex Confienement - 
	for (int j = 1; j < N_dim; j++)
	{
		for (int i = 1; i < N_dim; i++)
		{
			// Calc Gradient and Final VortCon Vector. 
			float W = vort->getdata(i, j);
			vec3<float> N = grad2d(i, j);

			// Write W (as vec3(W,W)) OR N (VC Gradient) to FluidObj vc grid, for Viz/Rendering via RenderObject. 
		   // f3obj->vc->setdata(vec3(W, W), i, j); // W VIZ
			f3obj->vc->setdata(vec3<float>(N) *= 1.0f , i, j); // N VIZ

			vec3<float> VC = N; VC *= W; 
			//f3obj->vc->setdata(VC, i, j); // N VIZ
			//VC *= (strength * dx);
			VC *= strength; // Final Vort Con Vector

			// Intergrate - 
			vec3<float> vel_0 = f3obj->vel->getdata(i, j);
			//vec3 vel_1 = vel_0 + (VC *= dt);
			N *= strength; //10.0f; // strength; // Mult By Stength Coeff. 
			vec3<float> vel_1 = vel_0 + (N *= (curl2d(i, j) * dt));

			// Set Vel - 
			f3obj->vel->setdata(vel_1, i, j);

		}
	}

	// Delete/Dealloc Vort Grid. 
	delete vort; vort = nullptr; 
}
*/

/* ====================================================
	VELOCITY PROJECTION / PRESSURE SOLVE 	
	==================================================== */

// Velocity Projection / Pressure Solve 
// Calculate Pressure Field, whose Gradient when subtraced from velocity field, will result in divergence ~0.
// Mass Conserving and Maintain Incompresilbitly, by subtracing divergence gradient from velocity field.
// Use an Iterative Linear System, to converge to the projected velocites. 
// Use HH Decomposition to Decompose Velocity Field.
// Uses Matrixless Laplcian A to solve uknown Pressure Vector x by caluclated diverence vector b. 

// Residual Vector = b - Ax

/* PROJECTION - GAUSS-SEIDEL + SUCESSIVE OVER RELAXATION (GS+SOR) - 
	Solved Pressure Values are re-injected into current iterations grids, hence this alogrithim has to be single threaded to avoid 
	data races across threads and allow correct synchroized reads of updated cells within the same iteration for reamaing cells. 
	Divergence and Pressure Gradient Subtraction loops are MT, these are thread safe (Only read or write from same grids). */

void fluidsolver_3::project(int iter)
{
	float h = 1.0f/N_dim; // CellSize 1.0f/N (N_dim); 

	// Ensure previous Divergence and Pressure Temp Grids have been deleted before new grid ptr assignement/Alloc. 
	del_divergence(); del_pressure(); 

	// Init Solvers Own Divergence and Pressure Fields, Only needed for Scope of this Function.  
	divergence = new grid3_scalar<float>(x_s, y_s, z_s, e_s);
	pressure = new grid3_scalar<float>(x_s, y_s, z_s, e_s);

	// DIVERGENCE FIELD CALC \\ -
	// Compute Divergence Field, from Velocity Field - 

	//!MT#pragma omp parallel for num_threads(omp_get_max_threads())

	for (int k = 1; k <= N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// Init to 0 
				divergence->setdata(0.0f, i, j, k);

				// Compute Divergence Cell Value. (0.5 * h oppose to / N) 
				float div = -0.5 * h * (f3obj->vel->getdata_x(i + 1, j, k) - f3obj->vel->getdata_x(i - 1, j, k)
					+ f3obj->vel->getdata_y(i, j + 1, k) - f3obj->vel->getdata_y(i, j - 1, k)
					+ f3obj->vel->getdata_z(i, j,k + 1) - f3obj->vel->getdata_z(i, j, k - 1));

				// Set Divergence Cell Value. 
				divergence->setdata(div, i, j, k);

				// Zero Out Pressure Grid, as Inital Value PreLinSolve. (Index Based oppose to calling grid3_scalar<float>->clear()).
				pressure->setdata(0.0f, i, j, k);

				// Write Inital PreProjected VelField to Grid For dbg - 
				//f3obj->preproj_vel->setdata(f3obj->vel->getdata(i, j, k), i, j, k);
			}
		}
	}
	// Call Boundary Condtions on Divergence Field and Inital Pressure Field - 
	edge_bounds(divergence); edge_bounds(pressure);
	#if dospherebound == 1
	sphere_bounds_eval(divergence, spherebound_coliso);
	sphere_bounds_eval(pressure, spherebound_coliso);
	#endif

	// PRESSURE FIELD LINEAR SOLVE \\ 

	// (Iterativly Compute Inverse of Divergence Grid/Matrix) Gauss-Seidel, Solve Discrete Poission Of Pressure Field, using Linear System. 
	for (int l = 0; l < iter; l++)
	{
		double error = 0.0f; 
		double total_error = 0.0f; 

		for (int k = 1; k <= N_dim; k++)
		{
			for (int j = 1; j <= N_dim; j++)
			{
				for (int i = 1; i <= N_dim; i++)
				{
					float n0 = pressure->getdata(i, j, k); // n (Pressure (l(n)) for SOR)

					float pres = (divergence->getdata(i, j, k) 
						+ pressure->getdata(i - 1, j, k) + pressure->getdata(i + 1, j, k) 
						+ pressure->getdata(i, j - 1, k) + pressure->getdata(i, j + 1, k)
						+ pressure->getdata(i, j, k - 1) + pressure->getdata(i, j, k + 1)
						) / 6.0f;

					float n1 = pres; // n+1 (Pressure (l(n+1)) for SOR)

					// Use SOR or not - 
					if (Parms.p_ProjectionType == Parms.Project_GaussSeidel_SOR)
					{
						// SOR 
						float alpha = Parms.p_SOR_alpha;
						float sor_pres = alpha * n1 + (1 - alpha) * n0;
						pressure->setdata(sor_pres, i, j, k);
					}
					else if (Parms.p_ProjectionType == Parms.Project_GaussSeidel)
					{
						pressure->setdata(pres, i, j, k);
					}

					// Caluclate Resdiual - 
					if (l == iter - 1)
					{
						// Check Error Value of Resulting Pressure Solve Convergence (b-Ax) wip -
						error += std::fabsf(divergence->getdata(i, j, k)) - std::fabsf(pressure->getdata(i, j, k));
					}

				}
			}
		}

		// Call Boundary Condtion Functions On Pressure After Each Pressure Field Iteration.
		edge_bounds(pressure);
		#if dospherebound == 1
		//sphere_bounds_scalar(pressure, spherebound_radius, spherebound_coliso, spherebound_offset);
		sphere_bounds_eval(pressure, spherebound_coliso);
		#endif	

		// Get Avg Residual Error on Last Iter. 
		if (l == iter - 1)
		{
			// Average Total Error (Over num of Cells) and Print. 
			total_error = std::fabsf(error / (float(pow(N_dim, 2))));
			std::cout << "DBG::PRESSURE SOLVER ERROR = " << std::scientific << std::fabsf(total_error) << "\n";
		}

	}
	
	// SUBTRACT PRESSURE GRADEINT FROM VELOCITY FIELD \\  

	//!MT#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int k = 1; k <= N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// 3 Methods Equivalent. 
				//1 Partial Derivatves for Each Pressure Gradient Components (2003 StableFluids)-
				float grad_x = 0.5f * (pressure->getdata(i + 1, j, k) - pressure->getdata(i - 1, j, k)) / h;
				float grad_y = 0.5f * (pressure->getdata(i, j + 1, k) - pressure->getdata(i, j - 1, k)) / h;
				float grad_z = 0.5f * (pressure->getdata(i, j, k + 1) - pressure->getdata(i, j, k - 1)) / h;

				//2 Typical Central Diffrence (/2h) - 
				//grad_x = (pressure->getdata(i + 1, j) - pressure->getdata(i - 1, j)) / (2.0f*h);
				//grad_y = (pressure->getdata(i, j + 1) - pressure->getdata(i, j - 1)) / (2.0f*h);

				//3 Stams Book Implmenetation - 
				//grad_x = 0.5f * N_dim * (pressure->getdata(i + 1, j) - pressure->getdata(i - 1, j));
				//grad_y = 0.5f * N_dim * (pressure->getdata(i, j + 1) - pressure->getdata(i, j - 1));

				// Subtract Gradient Components from Velocity Field Components and set to Velocity-
				float new_vel_x = f3obj->vel->getdata_x(i, j, k) - grad_x;
				f3obj->vel->setdata_x(new_vel_x, i, j, k);
				float new_vel_y = f3obj->vel->getdata_y(i, j, k) - grad_y;
				f3obj->vel->setdata_y(new_vel_y, i, j, k);
				float new_vel_z = f3obj->vel->getdata_z(i, j, k) - grad_z;
				f3obj->vel->setdata_z(new_vel_z, i, j, k);

			}
		}
	}
	// Call and Enforce Boundary Condtions After Projection on Vel Field - 
	edge_bounds(f3obj->vel);
	#if dospherebound == 1
	sphere_bounds_eval(f3obj->vel, spherebound_coliso);
	#endif
		
	// Delete Solver Temp Pressure and Divergence Fields - 
	del_divergence(); del_pressure();

}
// End of Velocity Projection Implementation (GS + SOR).


/* PROJECTION - JACOBI to Solve Pressure Poission Equation -
	Allows Multithreading as Cells (inner i,j presure loop) Only lookup Previous Pressure Values, but results in slower Convergence. 
	Divergence and Pressure Gradient Subtraction loops are MT, these are thread safe.*/

void fluidsolver_3::project_jacobi(int iter)
{
	float h = 1.0f/N_dim; // CellSize 1.0f/N (N_dim); 

	// Ensure previous Divergence and Pressure Temp Grids have been deleted before new grid ptr assignement/Alloc. 
	del_divergence(); del_pressure();
	
	// Alloc New DP Grids. 
	divergence = new grid3_scalar<float>(x_s, y_s, z_s, e_s);
	pressure = new grid3_scalar<float>(x_s, y_s, z_s, e_s);
	pressure_1 = new grid3_scalar<float>(x_s, y_s, z_s, e_s);

	// DIVERGENCE FIELD CALC \\ - 

	// Compute Divergence Field, from Velocity Field - 
	//!MTpragma omp parallel for num_threads(omp_get_max_threads())

	#pragma omp parallel for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				// Init to 0 
				divergence->setdata(0.0f, i, j, k);

				// Compute Divergence Cell Value. 
				float div = -0.5 * h * (f3obj->vel->getdata_x(i + 1, j, k) - f3obj->vel->getdata_x(i - 1, j, k)
					+ f3obj->vel->getdata_y(i, j + 1, k) - f3obj->vel->getdata_y(i, j - 1, k)
					+ f3obj->vel->getdata_z(i, j, k + 1) - f3obj->vel->getdata_z(i, j, k - 1));

				// Set Divergence Cell Value. 
				divergence->setdata(div, i, j, k);

				// Zero Out Pressure Grid, as Inital Value PreLinSolve. (Index Based oppose to calling grid3_scalar<float>->clear()).
				pressure->setdata(0.0f, i, j, k);
			}
		}
	}
	// Call Boundary Condtions on Divergence Field and Inital Pressure Field - 
	edge_bounds(divergence); edge_bounds(pressure); // Edge Bounds
	#if dospherebound == 1 // Sphere Bounds
	sphere_bounds_eval(divergence, spherebound_coliso);
	sphere_bounds_eval(pressure, spherebound_coliso);
	#endif

	// PRESSURE FIELD LINEAR SOLVE \\ 

	// (Iterativly Compute Inverse of Divergence Grid/Matrix) 
	// Jacobi For Possion Pressure Solve. (Write to Seperate "Sratch Grid", Only Swap at end of lth Iteration.
	// pressure == main pressure grid (to read from), pressure_1 == scratch pressure grid (to write to and then swap). 

	for (int l = 0; l < iter; l++)
	{
		#pragma omp parallel for
		for (int k = 1; k <= N_dim; k++)
		{
			#pragma omp parallel for
			for (int j = 1; j <= N_dim; j++)
			{
				#pragma omp parallel for
				for (int i = 1; i <= N_dim; i++)
				{
					float pres_n0 = pressure->getdata(i, j, k);

					float pres = (divergence->getdata(i, j, k) 
						+ pressure->getdata(i - 1, j, k) + pressure->getdata(i + 1, j, k)
						+ pressure->getdata(i, j - 1, k) + pressure->getdata(i, j + 1, k)
						+ pressure->getdata(i, j, k - 1) + pressure->getdata(i, j, k + 1)
						) / 6.0f;

					pressure_1->setdata(pres, i, j, k);
				}
			}
		}

		// Call Boundary Condtion Functions On Pressure After Each Pressure Field Iteration.
		edge_bounds(pressure_1);
		#if dospherebound == 1
		//sphere_bounds_scalar(pressure_1, spherebound_radius, spherebound_coliso, spherebound_offset);
		sphere_bounds_eval(pressure_1,spherebound_coliso); // Optimized SphereBounds for pressure calc. 
		#endif	

		// Swap Pressure Grid with Scratch Grid at end of k Iter After internal i,j,k MT'd Jacobi Projection Iteration is complete.
		// (This is Single Threaded to ensure whole grid is swapped correctly together, and not within a Multithreaded Outer loop).
		pressure_1->swap(pressure);
		// Hence removal of Outer MT Kth Loop and OMP Crticial which was not correct, as Pressure grid was swapping atomically per x threads, not as a singlethread. 
	}
	
	// SUBTRACT PRESSURE GRADEINT FROM VELOCITY FIELD \\ -
	#pragma omp parallel for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				// Partial Derivatves for Each Pressure Gradient Components -
				float grad_x = 0.5 * (pressure->getdata(i + 1, j, k) - pressure->getdata(i - 1, j, k)) / h;
				float grad_y = 0.5 * (pressure->getdata(i, j + 1, k) - pressure->getdata(i, j - 1, k)) / h;
				float grad_z = 0.5 * (pressure->getdata(i, j, k + 1) - pressure->getdata(i, j, k - 1)) / h;

				// Subtract Gradient Components from Velocity Field Components and set to Velocity-
				float new_vel_x = f3obj->vel->getdata_x(i, j, k) - grad_x;
				f3obj->vel->setdata_x(new_vel_x, i, j, k);
				float new_vel_y = f3obj->vel->getdata_y(i, j, k) - grad_y;
				f3obj->vel->setdata_y(new_vel_y, i, j, k);
				float new_vel_z = f3obj->vel->getdata_z(i, j, k) - grad_z;
				f3obj->vel->setdata_z(new_vel_z, i, j, k);
			}
		}
	}
	// Call and Enforce Boundary Condtions After Projection on Vel Field - 
	edge_bounds(f3obj->vel);

	#if dospherebound == 1
	sphere_bounds_eval(f3obj->vel, spherebound_coliso);
	#endif

	// TEMP FIELD DELETION/DEALLOCATION - 
	del_divergence(); del_pressure();
	
}
// End of Velocity Projection (Jacobi) Implementation.

/////////////////////////////////////////////////////////////////////////////

/* !TODO 
// PROJECT - Gauss-Seidel + SOR - SIMD TESTING - 

void fluidsolver_3::project_SIMD(int iter)
{
	float h = 1.0f / N_dim; // CellSize 1.0f/N (N_dim); 

	// Ensure previous Divergence and Pressure Temp Grids have been deleted before new grid ptr assignement/Alloc. 
	del_divergence(); del_pressure();

	// Init Solvers Own Divergence and Pressure Fields, Only needed for Scope of this Function.  
	divergence = new grid3_scalar<float>(x_s, y_s, e_s, 4, 1.0f);
	pressure = new grid3_scalar<float>(x_s, y_s, e_s, 5, 1.0f);

	// DIVERGENCE FIELD CALC \\ - 

	// Compute Divergence Field, from Velocity Field - 
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Init to 0 
			divergence->setdata(0.0f, i, j);

			// Compute Divergence Cell Value. 
			float div = -0.5 * h * (f3obj->vel->getdata_x(i + 1, j) - f3obj->vel->getdata_x(i - 1, j)
				+ f3obj->vel->getdata_y(i, j + 1) - f3obj->vel->getdata_y(i, j - 1));

			// Set Divergence Cell Value. 
			divergence->setdata(div, i, j);

			// Zero Out Pressure Grid, as Inital Value PreLinSolve. (Index Based oppose to calling grid3_scalar<float>->clear()).
			pressure->setdata(0.0f, i, j);

			// Write Inital PreProjected VelField to Grid For dbg - 
			f3obj->preproj_vel->setdata(f3obj->vel->getdata(i, j), i, j);

		}
	}
	// Call Boundary Condtions on Divergence Field and Inital Pressure Field - 
	edge_bounds(divergence); edge_bounds(pressure);
	#if dospherebound == 1
	sphere_bounds_eval(divergence, spherebound_coliso);
	sphere_bounds_eval(pressure, spherebound_coliso);
	#endif

	/*
	// PRESSURE FIELD LINEAR SOLVE \\ 

	std::vector<float>* pres_transpose = new std::vector<float>(total_size, 0.0f); // To Store Transposed Pressure grid_data std::vector<float>

	// (Iterativly Compute Inverse of Divergence Grid/Matrix) Gauss-Seidel, Solve Discrete Poission Of Pressure Field, using Linear System. 
	for (int k = 0; k < iter; k++)
	{
		double error = 0.0f;
		double total_error = 0.0f;

		// Get Transposed Grid with Jth Neighbour Elements in Contigous RowMajor (per k)-
		pressure->tranpose_gridata(pres_transpose);

		// Iterate over +2 Cells Per dim, i,i-1,j,j-1 And then Do i-1, i-1-1 etc...  
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				
				//float n0 = pressure->getdata(i, j); // n (Pressure (k(n)) for SOR)
				// NON SIMD Pressure Solve - 

				//std::size_t idx1d = (std::size_t) pressure->idx_2Dto1D(i, j);
				// Fast Loading to Reigster - IF Grid Mem is Contigous for cells needed (its not !).
				//avx256 t_a = _mm256_load_ps((pressure->getdataarray()) + idx1d); 

				// i-1 | i | i+1 | i+2
				__m128 ii = _mm_load_ps(pressure->getdataarray() + pressure->idx_2Dto1D(i-1, j) );
				// i-1 | 0.0f | i+1 | 0.0f
				ii.m128_f32[1] = 0.0f, ii.m128_f32[3] = 0.0f; 
				
				
				// Use Transposed j (So use same 2d IDX bias i) (So Get j-1 -> j+2 Contigously rowmajor transposed jth (j,i)). 
				// j-1 | j | j+1 | j+2
				__m128 jj = _mm_load_ps((pres_transpose->data()) + pressure->idx_2Dto1D(i-1, j) );
				// j-1 | 0.0f | j+1 | Divergence(i,j)
				jj.m128_f32[1] = 0.0f; ii.m128_f32[3] = divergence->getdata(i, j); // Oppose to zeroiing [3] add Divergence Cell Here. 
				

				// i-1 | 0.0f | i+1 | 0.0f
				//                                      +  
				// j-1 | 0.0f | j+1 | Divergence(i,j)

				// Cast 2 m128 vectors to single _m256 vector. (3 0 unsused elements)  
				// 0        1       2      3      4       5      6           7
				// i-1,j | 0.0f | i+1,j | 0.0f | i,j-1 | 0.0f | i,j+1 | Divergence(i,j)|  With Both Dim Indices for reference. 

				// Naivie Way - 
				//__m256 iijj = _mm256_set_ps(ii.m128_f32[0], ii.m128_f32[1], ii.m128_f32[2], ii.m128_f32[3],
				//	jj.m128_f32[0], jj.m128_f32[1], jj.m128_f32[2], jj.m128_f32[3]);

				// Intrinsics Way-
				__m256 iijj = _mm256_castps128_ps256(ii); // Cast ii to m256
				 iijj = _mm256_insertf128_ps(iijj, jj, 1); // inster jj into latter element (idx 1) for resulting iijj comb. 

				// 256 bit reg Hadditon - 
				// Call SIMD _256 HADD and div over 4. 
				float pres = (simd256_hAdd(iijj)) / 4.0f;

				// Set Pressure to resulting float. 
				pressure->setdata(pres, i, j); 

				// Slow Loading Instructions - 
				//avx256 t_a = _mm256_set_ps(divergence->getdata(i, j), pressure->getdata(i - 1, j), pressure->getdata(i + 1, j), pressure->getdata(i, j - 1), pressure->getdata(i, j + 1),
				//0.0f, 0.0f, 0.0f); 

				// Add Press
				//float pres_ij = (divergence->getdata(i, j) + pressure->getdata(i - 1, j) + pressure->getdata(i + 1, j) +
				//	pressure->getdata(i, j - 1) + pressure->getdata(i, j + 1)) / 4.0f;

				//float n1 = pres; // n+1 (Pressure (k(n+1)) for SOR)

				/*
				// Use SOR or not - 
				if (Parms.p_ProjectionType == Parms.Project_GaussSeidel_SOR)
				{
					// SOR 
					float alpha = Parms.p_SOR_alpha;
					float sor_pres = alpha * n1 + (1 - alpha) * n0;
					pressure->setdata(sor_pres, i, j);
				}
				*/
				//else if (Parms.p_ProjectionType == Parms.Project_GaussSeidel) {}
		
				// pressure->setdata(t_b, i, j);
				

				/*
				// Caluclate Resdiual - 
				if (k == iter - 1)
				{
					// Check Error Value of Resulting Pressure Solve Convergence (b-Ax) wip -
					error += std::fabsf(divergence->getdata(i, j)) - std::fabsf(pressure->getdata(i, j));
				}
				
				

			}
		}


		// Call Boundary Condtion Functions On Pressure After Each Pressure Field Iteration.
		edge_bounds(pressure);
		#if dospherebound == 1
		//sphere_bounds_scalar(pressure, spherebound_radius, spherebound_coliso, spherebound_offset);
		sphere_bounds_eval(pressure, spherebound_coliso);
		#endif	

		// Get Avg Residual Error on Last Iter. 
		if (k == iter - 1)
		{
			// Average Total Error (Over num of Cells) and Print. 
			total_error = std::fabsf(error / (float(pow(N_dim, 2))));
			std::cout << "DBG::PRESSURE SOLVER ERROR = " << std::scientific << std::fabsf(total_error) << "\n";
		}
	}
	*/
	// SUBTRACT PRESSURE GRADEINT FROM VELOCITY FIELD \\  

	// Delete Transposed Grid_Data
	//delete pres_transpose; pres_transpose = nullptr;
/* !TODO
	#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// 3 Methods Equivalent. 
			//1 Partial Derivatves for Each Pressure Gradient Components (2003 StableFluids)-
			float grad_x = 0.5 * (pressure->getdata(i + 1, j) - pressure->getdata(i - 1, j)) / h;
			float grad_y = 0.5 * (pressure->getdata(i, j + 1) - pressure->getdata(i, j - 1)) / h;

			//2 Typical Central Diffrence (/2h) - 
			//grad_x = (pressure->getdata(i + 1, j) - pressure->getdata(i - 1, j)) / (2.0f*h);
			//grad_y = (pressure->getdata(i, j + 1) - pressure->getdata(i, j - 1)) / (2.0f*h);

			//3 Stams Book Implmenetation - 
			//grad_x = 0.5f * N_dim * (pressure->getdata(i + 1, j) - pressure->getdata(i - 1, j));
			//grad_y = 0.5f * N_dim * (pressure->getdata(i, j + 1) - pressure->getdata(i, j - 1));

			// Subtract Gradient Components from Velocity Field Components and set to Velocity-
			float new_vel_x = f3obj->vel->getdata_x(i, j) - grad_x;
			f3obj->vel->setdata_x(new_vel_x, i, j);
			float new_vel_y = f3obj->vel->getdata_y(i, j) - grad_y;
			f3obj->vel->setdata_y(new_vel_y, i, j);

		}
	}

	// Call and Enforce Boundary Condtions After Projection on Vel Field - 
	edge_bounds(f3obj->vel);
	#if dospherebound == 1
	sphere_bounds_eval(f3obj->vel, spherebound_coliso);
	#endif

	// Delete Solver Temp Pressure and Divergence Fields - 
	del_divergence(); del_pressure();

}
// End of Velocity Projection Implementation (GS + SOR) SIMD.
*/

/* ====================================================
	DESNITY SOLVE STEP - 
	==================================================== */
// Implementation of Simulation Solve step of Density solve operations. 

void fluidsolver_3::density_step(int diff_iter, float diffA, bool dodiff)
{
	// DIFFUSE Density - 
	// If Using Diffusion, If Not Need to Manually set Cur to Prev, rather than swapping.
	if (Parms.p_Do_Dens_Diff == true)
	{
		f3obj->dens->swap(f3obj->prev_dens); // Swap Density With Prev_Density. 
		// Gauss-Seidel Iterative Density (Scalar) Diffusion - 
		diffuse(f3obj->prev_dens, f3obj->dens, diffA, diff_iter);

		// Finite Diffrence Unstable Density (Scalar) Diffusion - 
		//diffuse_scalar_FDM(f3obj->prev_dens, f3obj->dens, diffA);
		f3obj->dens->swap(f3obj->prev_dens); // Re-Swap Density With Prev_Density. 
	}
	else if (Parms.p_Do_Dens_Diff == false)
	{
		// Use "SetCurToPrev" Funcs to Copy Grids, Oppose to via Diffusion - 
		f3obj->setcurtoprev(f3obj->prev_dens, f3obj->dens);
	}

	// ADVECT Density - 

	if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_Euler)
	{
		advect_sl(f3obj->prev_dens, f3obj->dens); 
	}
	else if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_RK2)
	{
		advect_sl_mp(f3obj->prev_dens, f3obj->dens);
	}

	if (Parms.p_Do_Dens_Disp)
	{
		// DISSIPATE Density by multipler - 
		dissipate(f3obj->dens, Parms.p_Dens_Disp_Mult, this->dt); // Density Dissipation. 
	}

}
// End of Density Step Implemetnation.

/*	====================================================
	Velocity SOLVE STEP -
	==================================================== */

// Implementation of Simulation Solve step of Velocity solve operations. 
// Excluding Sourcing which will be done via user in FluidObject.
// Removed Pre Advection Projection Calls, So I can use One Call Post Advection, With Higher Iter Counts.

void fluidsolver_3::velocity_step(int diff_iter, int proj_iter, float diffA, bool dodiff)
{
	// If Using Diffusion, If Not Need to Manually set Cur to Prev, rather than swapping. 
	if (dodiff == true)
	{
		f3obj->vel->swap(f3obj->prev_vel); // Swap Vel Field with Prev_VelField.
		diffuse(f3obj->prev_vel, f3obj->vel, diffA, diff_iter);
		f3obj->vel->swap(f3obj->prev_vel); // Re-Swap Vel With Prev_Vel. 
	}
	else if (dodiff == false)
	{
		f3obj->setcurtoprev(f3obj->prev_vel, f3obj->vel);
	}

	// ADVECT VELOCITY FIELD (Self Advect) \\

	if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_Euler)
	{
		advect_sl(f3obj->prev_vel, f3obj->vel);
	}
	else if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_RK2)
	{
		advect_sl_mp(f3obj->prev_vel, f3obj->vel);
	}

	if (Parms.p_Do_Vel_Disp)
	{
		// DISSIPATE Density by multipler - 
		dissipate(f3obj->vel, Parms.p_Vel_Disp_Mult, this->dt); // Density Dissipation. 
	}

	// VORTICITY CONFINEMENT WIP \\

	if (Parms.p_useVorticity)
	{
		/* WIP
		vorticty_confine(1.0f);
		vorticty_confine_otf(5.0f);
		*/
		vorticity_confine_B(5.0f);

	}

	// PROJECT VELOCITY FIELD \\  (Post Advect Only) 

	if (Parms.p_ProjectionType == Parms.Project_GaussSeidel || Parms.p_ProjectionType == Parms.Project_GaussSeidel_SOR)
	{
		// GS + SOR (ST) - (SOR May or May not be enabled by FluidSolver Member Paramter). 
		project(Parms.p_GS_Proj_iter);

		// SIMD Projection WIP 
		//project_SIMD(Parms.p_GS_Proj_iter);
	}
	else if (Parms.p_ProjectionType == Parms.Project_Jacobi)
	{
		// Jacobi (MT) - 
		project_jacobi(Parms.p_Jacobi_Proj_Iter); 
	}

}
// End of Velocity Step Implemetnation.


/*	====================================================
	SOLVE STEP -
	==================================================== */

// Implementation of A Single solve Step of fluidsolver_3. Once User Calls this
// Simulation will begin. And for now occur Untill program is closed. 

void fluidsolver_3::solve_step(bool solve, bool do_diffdens, bool do_diffvel, float dens_diff, float vel_diff, int proj_iter, int diff_iter, int max_step)
{
	// Init Step/time vals. 
	int step_count = 0;
	float dt_acc = 0.0f; // Accumulated DeltaTime (SimTime)
	float total_elapsed = 0.0f; // Chrono Time

	// Render Object Creation/Setup - 
	// Pass Render Context Window Ptr Now FluidSovler2Mem, to RenderObject. (Passed in To Solver in Main via set_window() MF.
	int render_mode = 1; // 0 = Density, 1 = Vel. // May be overriden by Input in SolveStep (RenderObj::ShaderPipe()).
	render_obj = new renderobject_3D_OGL("OpenGL", 4, 0, x_s + e_s, y_s + e_s, this->winptr, render_mode); 
	// Ideally Move render_obj setup to main/outside, and then pass into fluidsolver...?

	// Move this into RenderObject Initalization - 
	// Set Render_Object Texture Units to Sampler Uniforms Pre Solve Loop - 
	glUseProgram(render_obj->shader_prog); // Call Use (Shader) Program First. 

	// Set 3D Tex Sampler Uniforms, to matching Texture Units.
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "d_tex"), 0); // Density = 0. 


	//Solve Step And Render - EXEC LOOP 
	while (solve == true && step_count <= max_step) // Infinite Solve Loop For Now. Will need to link this to drawing/waiting etc. 
	{
		// Time Start - 
		std::chrono::system_clock::time_point timestep_start = std::chrono::system_clock::now();
		std::stringstream log_out;

		// STEP INPUT OPERATIONS \\ ----------------
		// Interactive Sourcing ... 
		// Interactive Sphere Bounds Radius Eval. 
		sphere_rad_test(); // Can Cause Pressure Solver Crashes. 

		// Interp Switching.
		//if (glfwGetKey(winptr, GLFW_KEY_I) == GLFW_PRESS) { Parms.p_InteroplationType == Parms.Interoplation_Linear; };
		//if (glfwGetKey(winptr, GLFW_KEY_C) == GLFW_PRESS) { Parms.p_InteroplationType == Parms.Interoplation_Cosine; };

		// Get CurFrame Mouse Pos And Update Mouse Vel. Repos this call to avoid sleep delay. 
		std::this_thread::sleep_for(std::chrono::milliseconds(5)); // Hold main thread to allow for delta mouse position (VelCalc). 5ms works fine  10ms
		updt_mousepos(step::STEP_CUR); updt_mouseposNorm(step::STEP_CUR); updt_mouseposRange(step::STEP_CUR);
		updt_mousevel(); 

		// Pass to Override SphereBound_offset Values For Sphere Cols. 
		spherebound_offset.x = xpos_1_N; spherebound_offset.y = ypos_1_N;
		//std::cout << "DEBUG:: MOUSE_X = " << xpos_1 << " MOUSE_Y = " << ypos_1 << "\n";

		// PRESTEP OPERATIONS \\ ----------------
		// Eval SphereBounds_SDF - 
		sphere_bounds_set(spherebound_radius, spherebound_coliso, spherebound_offset);

		// STEP SOURCING OPERATIONS \\ ----------------
		//vec3 anim_offset(0.4 + (float(sin(float(step_count) / float(max_step) * 50.0f))) * 0.1f, 0.4 + (float(cos(float(step_count) / float(max_step) * (float(step_count) / float(max_step) * 10.0f)))) * 0.2f);
		f3obj->implicit_sphere_source(0.5f, vec3<float>(0.0f, 0.5f, 0.0f), vec3<float>(0.5f, 0.5f, 0.5f), 0.01f);

		// Forces- 
		//if (step_count <= 20) f3obj->radial_force(vec3<float>(0.499f, 0.499f), 0.8f, this->dt);

		// STEP - SUB - SOLVER STEP OPERATIONS \\ -------------- 
		velocity_step(diff_iter, proj_iter, vel_diff, do_diffvel);
		density_step(diff_iter, dens_diff, do_diffdens);

		// STEP - RENDER CALLS \\ ------------------
		// Pass Cur Step - 
		render_obj->et = step_count; 
		// Uniform/Texture Set Calls (ShaderPipe) -
		render_obj->shader_pipe(f3obj); // Pass Grids Via Friend Acess to shader_pipe() MF
		// Render Operations Call - 
		render_obj->call_ren(rend_state::RENDER_ACTIVE); // Call Render Step Ops here within solver loop, ie NON Debug Mode (Pass RENDER_ACTIVE).
		#endif	


		// STEP DEBUG CONSLE OPERATIONS \\ -------------------- 
		std::chrono::system_clock::time_point timestep_end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed = timestep_end - timestep_start;
		double elapsed_sec = elapsed.count();
		total_elapsed += elapsed_sec;

		// Store Mouse Pos as Prev Frame Pos.. 
		updt_mousepos(step::STEP_PREV); updt_mouseposNorm(step::STEP_PREV); updt_mouseposRange(step::STEP_PREV);

		// Differintate between SimTime Passed of accumlated DeltaTime, Vs Duration of SolveStep Caluclation (WallClock) !
		log_out << "INFO::Solve_Step::SimTime Passed = " << dt_acc << " Seconds of " << dt << " timestep \n";
		log_out << "DEBUG::Solve_Step::Step_Count =  " << step_count << "\n";
		log_out << "DEBUG::Solve_Step::Calc_Duration = " << elapsed_sec << " Seconds" << "\n";

		if (step_count == max_step)
		{
			log_out << "\n *** INFO::Solve_Step::TOTAL_DELTA-TIME = " << dt_acc << " Seconds ***\n";

			//std::cout << std::fixed; 
			log_out << " \n *** DEBUG::Solve_Step::TOTAL_CALC_DURATION = " << total_elapsed << " Seconds *** \n";
			log_out << " \n *** DEBUG::Solve_Step::TOTAL_FPS = ~" << double(max_step) / total_elapsed << " FPS *** \n \n";

			// Write Last Frame To Img 
			//f3obj->writeto_img(step_count);
			//std::cout << "\n"; f3obj->print_info();
			//f3obj->writeto_img_vel(step_count);
		}

		// Print Log to stdout - 
		std::cout << log_out.str() << std::endl; 

		// STEP INCREMENT \\ ------------
		dt_acc += dt; 
		step_count++; // Inc Step_Count
	}
}
// End of Solve Step Implementation.


/*	====================================================
	Utility Member Functions - 
	==================================================== */

// Window Setter. Passed from Render Context Object.  
void fluidsolver_3::set_window(GLFWwindow *win)
{
	assert(win != nullptr); // Check for Passing a nullptr FluidSolver GLFW WinPtr. Assert this.
	this->winptr = win; 
}

// Call To Update Mouse Pos - (Pixel (Cell Index) Space)
void fluidsolver_3::updt_mousepos(const step step_id)
{
	glfwPollEvents();

	if (step_id == step::STEP_CUR)
	{
		glfwGetCursorPos(winptr, &xpos_1, &ypos_1);
		// Window Bounds - 
		if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_1 = 0.0f, ypos_1 = 0.0f;
	}
	else if (step_id == step::STEP_PREV)
	{
		glfwGetCursorPos(winptr, &xpos_0, &ypos_0);
		// Window Bounds - 
		if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_0 = 0.0f, ypos_0 = 0.0f;
	}
}

// Call To Update Mouse Pos Norm (0-1 XY) - (Normalized Window / GridSpace).
void fluidsolver_3::updt_mouseposNorm(const step step_id)
{
	glfwPollEvents();
	if (step_id == step::STEP_CUR)
	{
		glfwGetCursorPos(winptr, &xpos_1_N, &ypos_1_N);

		ypos_1_N = (y_s + e_s) - ypos_1_N; // FLIP Y Axis (So Matches GridSpace Postive Y = Down) (Up When drawn). 

		xpos_1_N /= N_dim; ypos_1_N /= N_dim;

		// Window Bounds - 
		if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_1_N = 0.0f, ypos_1_N = 0.0f;
	}
	else if (step_id == step::STEP_PREV)
	{
		glfwGetCursorPos(winptr, &xpos_0_N, &ypos_0_N);

		ypos_0_N = (y_s + e_s) - ypos_0_N; // FLIP Y Axis (So Matches GridSpace Postive Y = Down) (Up When drawn). 

		xpos_0_N /= N_dim; ypos_0_N /= N_dim;

		// Window Bounds - 
		if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_0_N = 0.0f, ypos_0_N = 0.0f;
	}
}

// Call To Update Mouse Pos Norm Within Clamped Range
// NOT USED NOW Y AXIS issues are fixed within updt_mouseposNorm. Will keep for now. 

void fluidsolver_3::updt_mouseposRange(const step step_id)
{
	glfwPollEvents();
	if (step_id == step::STEP_CUR)
	{
		glfwGetCursorPos(winptr, &xpos_1_R, &ypos_1_R);
		xpos_1_R /= N_dim; ypos_1_R /= N_dim;

		xpos_1_R = solver_utils::fit(xpos_1_R, 0.0f, 1.0f, -1.0f, 1.0f);
		ypos_1_R = solver_utils::fit(ypos_1_R, 0.0f, 1.0f, -1.0f, 1.0f);

		// Window Bounds - 
		if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_1_R = 0.0f, ypos_1_R = 0.0f;
	}
	else if (step_id == step::STEP_PREV)
	{
		glfwGetCursorPos(winptr, &xpos_0_R, &ypos_0_R);
		xpos_0_R /= N_dim; ypos_0_R /= N_dim;

		xpos_0_R = solver_utils::fit(xpos_1_R, 0.0f, 1.0f, -1.0f, 1.0f);
		ypos_0_R = solver_utils::fit(ypos_1_R, 0.0f, 1.0f, -1.0f, 1.0f);

		// Window Bounds - 
		if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_0_R = 0.0f, ypos_0_R = 0.0f;
	}
}

// Calc Mouse Velocity - 
void fluidsolver_3::updt_mousevel()
{
	// Prev and Cur Mouse Postions. 
	vec2<float> p_0(xpos_0_N, ypos_0_N);
	vec2<float> p_1(xpos_1_N, ypos_1_N);

	// Calc Vel (P Derivative) Diff over dt. 
	vec2<float> dp = p_1 - p_0; 
	dp /= dt;

	mouse_vel = dp; 

	//if (verbose) 
	std::cout << std::fixed << "DEBUG::Mouse Velocity = [" << mouse_vel.x << "," << mouse_vel.y << "] \n";
}

// Dissipate Grid Quanities By Multiplier - 
void fluidsolver_3::dissipate(grid3_scalar<float> *grid, float disp_mult, float dt)
{
	//disp_mult = std::max(0.0f, std::min(disp_mult, 1.0f)); // Enforce 0-1 Mult. 
	disp_mult = solver_utils::clamp(disp_mult, 0.0f, 1.0f); 

	#pragma omp parallel for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				// Don't Mult By dt for now. 
				float cur_scl = grid->getdata(i, j, k);
				cur_scl *= disp_mult;
				grid->setdata(cur_scl, i, j, k);
			}
		}
	}

}

void fluidsolver_3::dissipate(grid3_vector<vec3<float>> *grid, float disp_mult, float dt)
{
	disp_mult = solver_utils::clamp(disp_mult, 0.0f, 1.0f);
	#pragma omp parallel for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				// Don't Mult By dt for now. 
				vec3<float> cur_vel = grid->getdata(i, j, k);
				cur_vel *= disp_mult;
				grid->setdata(cur_vel, i, j, k);
			}
		}
	}
}



/*	====================================================
Solver Utils Implementation (Lambda Utility Funcs) -
==================================================== */

// Implement Static std::function stored Lambdas.

// CLAMP: Between Min and Max Range Lambda (Stored in std::function). 
std::function<float(float,float,float)> solver_utils::clamp = [&](float v, float min, float max) -> float
{
	return std::max(min, std::min(v, max));
};

// FIT: Value From Some Range A to Some New Range B (Eg 0 to 1, to -1 to 1 0.5 becomes 0) Lambda (Stored in std::function). 
std::function<float(float, float, float, float, float)> solver_utils::fit = [&](float val, float a_min, float a_max, float b_min, float b_max) -> float
{
	return b_min + (val - a_min)*(b_max - b_min) / (a_max - a_min);
};

// LERP: From One Float Value to Another By 0-1 Bias Value. 
std::function<float(float, float, float)> solver_utils::lerp = [&](float val_0, float val_1, float bias) -> float
{
	return (1.0f - bias) * val_0 + bias * val_1;
};

// 1D CosInterp. 
std::function<float(float, float, float)> solver_utils::cosinterp = [&](float val_0, float val_1, float bias) -> float
{
	float mu = (1.0f - std::cos(bias*PI)) / 2.0f;
	return (float) (val_0*(1.0f - mu) + val_1 * mu);
};


/*	====================================================
	DEBUG - Member FUnctions 
	==================================================== */

void fluidsolver_3::fill_test(int x)
{
	if (x >= NE_dim) x = NE_dim; 

	for (int k = 0; k < x; k++)
	{
		for (int j = 0; j < x; j++)
		{
			for (int i = 0; i < x; i++)
			{
				f3obj->dens->setdata(1.0f, i, j, k);
			}
		}
	}

}

// Test Implementation of Interactive Sphere Radius - Can Cause Pressure Solver Crashes. 
void fluidsolver_3::sphere_rad_test()
{
	if (glfwGetKey(winptr, GLFW_KEY_PAGE_UP) == GLFW_PRESS)
	{
		spherebound_radius += 0.001f; 
	}

	if (glfwGetKey(winptr, GLFW_KEY_PAGE_DOWN) == GLFW_PRESS)
	{
		spherebound_radius -= 0.0001f;
	}

	if (spherebound_radius <= 0.0001) spherebound_radius = 0.0001;
}
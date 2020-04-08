// Implementation of fluidsolver_3
#include "fluidsolver3d.h"

// Std Headers - 
#include <memory>
#include <chrono>
#include <thread>
#include <tuple>
#include <sstream>

#include <cassert>

// External Headers - 
#include <omp.h> 

#define dospherebound 0 
#define doedgebound 1 
#define DO_SPHEREBOUNDS_MT 1

#define RENDER_GL 1
#define RENDER_IMG 0

vec3<float> spherebound_offset(0.0f, 0.0f, 0.0f);  
const float spherebound_coliso = 0.0025f; 
float spherebound_radius = 0.005f; 
float impsource_radius = 0.008f;

extern int win_size_xy;
extern short verbose;


/* ====================================================
	fluidsolver_3 CONSTRUCTOR/DESTRUCTOR'S
   ==================================================== */

fluidsolver_3::fluidsolver_3(fluidobj_3d *f3dptr, float dtt) 
	: dt(dtt), f3obj(f3dptr) 
{
	// Get Dimensions,Spacing from FluidObject Ensure Values Match. 
	x_s = f3obj->x_s, y_s = f3obj->y_s, z_s = f3obj->z_s, e_s = f3obj->e_s;
	total_size = f3obj->t_s;

	N_dim =  f3dptr->x_s; // Single Dim Size, No Edge. 
	NE_dim = f3dptr->x_s + f3dptr->e_s; // Single Dim Size, With Edge Size. 

	// Init Solver (non temp) grids.
	spherebounds_sdf = new grid3_scalar<float>(N_dim, N_dim, N_dim, 2);

}

// fluidsolver_3 DESTRUCTOR - Toned down SmartPtr use for now. 
fluidsolver_3::~fluidsolver_3()
{
	delete spherebounds_sdf; spherebounds_sdf = nullptr; 

	del_pressure();  del_divergence();

	if (vort || vort != nullptr)
	{
		delete vort; vort = nullptr; 
	}

	delete render_obj; render_obj = nullptr; 

	// FluidObj Responsible for own deletion. 
}

/* ====================================================
	FIELD DELTETION
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


// ** EDGE_BOUNDS_SCALAR-FIELD IMPLEMENTATION ** \\ - 
// Assumes Grid is Cubed (N = X = Y = Z)
void fluidsolver_3::edge_bounds(grid3_scalar<float> *grid)
{
	// 6 Edge Faces of Grid (X-|X+|Y-|Y+|Z-|Z+)

	#pragma omp parallel for 
	for (int j = 1; j <= N_dim; j++)
	{
		#pragma omp parallel for 
		for (int i = 1; i <= N_dim; i++)
		{
			// X- Face Boundary 
			float x0i = grid->getdata(1, i, j); // [0,i,j] Boundary Values from Edge [1,i,j] Values.
			grid->setdata(x0i, 0, i, j);

			// X+ Edge Boundary
			float xNi = grid->getdata(N_dim, i, j); // [N+1,i,j] Boundary Values from Edge [N,i,j] Values.
			grid->setdata(xNi, (N_dim + 1), i, j);

			// Y- Edge Boundary
			float y0i = grid->getdata(i, 1, j); // [i, 0, j] Boundary Values from Edge [i, 1, j] Values.
			grid->setdata(y0i, i, 0, j);

			// Y+ Edge Boundary
			float yNi = grid->getdata(i, N_dim, j); // [i, N+1, j] Boundary Values from Edge [i,N,j] Values.
			grid->setdata(yNi, i, N_dim + 1, j);

			// Z- Edge Boundary
			float z0i = grid->getdata(i, j, 1); // [i, j, 0] Boundary Values from Edge [i,j,1] Values.
			grid->setdata(z0i, i, j, 0);

			// Z+ Edge Boundary
			float zNi = grid->getdata(i, j, N_dim); // [i, j, N+1] Boundary Values from Edge [i, j, N] Values.
			grid->setdata(zNi, i, j, N_dim + 1);

		}
	}


	/* 8 Corner Cells, ScalarGrid Edge Bounds Corner Adjacent Cell Neighbour Averages -
	   Self + or - XYZ (0(+1) or N+1(-1(N))) (0 and N+1 Axis Ghost/Edge Cells).
	   If At 0 For Coord Offset +1, if At N+1 For Coord Offset -1. Offset Axis, Keep Others Const. Like a PDE. */

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

// Set Grid Edge Bounds Faces, to reflect the perpendicular velocity component for cells within each six faces of grid. 
// ThreadSafe.
void fluidsolver_3::edge_bounds(grid3_vector<vec3<float>> *grid)
{
	// X +/- Grid Face Cells = Y Compoent Reflection
	// Y +/- Grid Face Cells = X Component Relfection
	// Z +/- Grid Face Cells = Z Component Relfection (Yes Z reflects Z)
	
	#pragma omp parallel for
	for (int j = 1; j <= N_dim; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= N_dim; i++)
		{
			// X -/+ Face Cells, Reflect Y Velocity Component -

			// X- Edge Boundary 
			float x0i = grid->getdata_y(1, i, j); // [0,i,j] Boundary Values from Edge [1,i,j] Values.
			x0i *= -1.0f;
			grid->setdata_y(x0i, 0, i, j);

			// X+ Edge Boundary
			float xNi = grid->getdata_y(N_dim, i, j); // [N+1,i,j] Boundary Values from Edge [N,i,j] Values.
			xNi *= -1.0f;
			grid->setdata_y(xNi, (N_dim + 1), i, j);

			// Y -/+ Face Cells, Reflect X Velocity Component -

			// Y- Edge Boundary
			float y0i = grid->getdata_x(i, 1, i); // [i, 0, i] Boundary Values from Edge [i, 1, j] Values.
			y0i *= -1.0f;
			grid->setdata_x(y0i, i, 0, j);

			// Y+ Edge Boundary
			float yNi = grid->getdata_x(i, N_dim, j); // [i, N+1, j] Boundary Values from Edge [i,N,j] Values.
			yNi *= -1.0f;
			grid->setdata_y(yNi, i, N_dim + 1, j);

			// Z -/+ Face Cells, Reflect Z Velocity Component -

			// Z- Edge Boundary
			float z0i = grid->getdata_z(i, j, 1); // [i, j, 0] Boundary Values from Edge [i,j,1] Values.
			z0i *= -1.0f;
			grid->setdata_z(z0i, i, j, 0);

			// Z+ Edge Boundary
			float zNi = grid->getdata_z(i, j, N_dim); // [i, j, N+1] Boundary Values from Edge [i, j, N] Values.
			zNi *= -1.0f;
			grid->setdata_z(zNi, i, j, N_dim + 1);
		}
	}

	// Corner Cell Interoplation/Averaging of Adjacent Neighbours for Ux,Vy,Wz Velocity Components -  
	// Do Each Corner Cell For each Component. Inline Averaging.

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

	#pragma omp parallel for
	for (int k = 1; k < N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				// Implicit Sphere/sphere Function: x^2 + y^2 + z^2 - r . >= 0 <= Thresh == Surface. > Thresh = Exterior. < 0 = Interior. Cells. 
				vec3<float> cell_gridSpace(float(i * h) - offset.x, float(j * h) - offset.y, float(k*h) - offset.z); // Index to Grid Space 0-1N. 
				float sphere_func = ((cell_gridSpace.x * cell_gridSpace.x) + (cell_gridSpace.y * cell_gridSpace.y) + (cell_gridSpace.z * cell_gridSpace.z)) - radius; 

				spherebounds_sdf->setdata(0.0f, i, j, k); 
				spherebounds_sdf->setdata(sphere_func, i, j, k); 
			}
		}
	}

}

// ** SPHERE_BOUNDS_EVAL - SCALAR-FIELD OVERRIDE ** \\ - 
// Eval SphereBounds - On Scalar Field. Unused (Surface and Exterior Condtions disabled for perf). 

void fluidsolver_3::sphere_bounds_eval(grid3_scalar<float> *grid, float col_iso)
{
	float h = 1.0f / N_dim; // Grid Spacing, Recoprical of One Dim Size (N). 
	
	//!MT #pragma omp parallel for num_threads(omp_get_max_threads())	

	#pragma omp parallel for
	for (int k = 1; k < N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
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

void fluidsolver_3::sphere_bounds_eval(grid3_vector<vec3<float>> *grid, float col_iso)
{
	float h = 1.0f / N_dim; // Grid Spacing, Recoprical of One Dim Size (N). 

	//!MT #pragma omp parallel for num_threads(omp_get_max_threads())
	#pragma omp parallel for
	for (int k = 1; k < N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				// Lookup Sphere/sphere SDF At CurCell i,j. 
				float sphere_func = spherebounds_sdf->getdata(i, j, k);

				/*
				// Surface
				if (sphere_func >= 0.0f && sphere_func <= col_iso) // func-radis > 0.0f but <= col_iso Cell is on "surface".
				{
					// Surface Operations -
				}
				*/

				// Inside - 
				if (sphere_func < 0.0f) // func-radius < 0.0f Cell is inside. 
				{
					// Inside Operations - 

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

				/* Mouse Vel Additon
				if (sphere_func <= col_iso) 
				{
					float mouse_velMult = 1.0f;
					// Get Stored Current Mouse Vel + to Current Grid Cell Velocity. 
					vec2<float> cur_vel = grid->getdata(i, j);
					grid->setdata(cur_vel + (mouse_vel *= mouse_velMult), i, j);
					grid->setdata(grid->getdata(i, j, k) + (mouse_vel *= mouse_velMult), i, j, k);
				}
				*/
			}
		}
	}

}


/* ====================================================
	DIFFUSION
   ==================================================== */
// Gauss-Seidel Based Diffusion, not Thread Safe Possible Race Condtion of Neighbour Cells own diffusion.

// ** DIFFUSION-SCALAR-FIELD-LINSOLVE IMPLEMENTATION ** \\ 

void fluidsolver_3::diffuse(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1, float diff, ushort iter)
{
	// Implicit Matrix
	float a = dt * diff * powf(N_dim, 2.0f);

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
		#if doedgebound == 1
		edge_bounds(grid_1); 
		#endif

		#if dospherebound == 1
		sphere_bounds_eval(grid_1, spherebound_coliso);
		#endif
	}
}

// ** DIFFUSION-VECTOR-FIELD-LINSOLVE IMPLEMENTATION ** \\ 

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
		#if doedgebound == 1
		edge_bounds(grid_1); // Generic Func Call, Pass Grid_1 (n+1). 
		#endif

		#if dospherebound == 1
		sphere_bounds_eval(grid_1, spherebound_coliso);
		#endif
	}
}

/* ====================================================
	ADVECTION
   ==================================================== */

// SEMI LAGRANGIAN - SINGLE STEP (EULER) ADVECTION \\

// Semi Lagrangian Advection - Single Step (Euler) - Scalar Overload - 
void fluidsolver_3::advect_sl(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1)
{
	// Tri-Cosine Interoplation of Grid Cell Indices ijk(0|+1) with Coefficents (t,s,r). 
	auto interpCos_scalar = [&](int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) -> float
	{
		// Tri-Cosine Interoplation .
		float C_000_001_t = cosinterp(grid_0->getdata(i0, j0, k0), grid_0->getdata(i0, j0, k1), t1);
		float C_010_011_t = cosinterp(grid_0->getdata(i0, j1, k0), grid_0->getdata(i0, j1, k1), t1);
		float C_100_101_t = cosinterp(grid_0->getdata(i1, j0, k0), grid_0->getdata(i1, j0, k1), t1);
		float C_110_111_t = cosinterp(grid_0->getdata(i1, j1, k0), grid_0->getdata(i1, j1, k1), t1);
		float C_A = cosinterp(C_000_001_t, C_010_011_t, s1);
		float C_B = cosinterp(C_100_101_t, C_110_111_t, s1);
		float C_F = cosinterp(C_A, C_B, r1);

		return C_F;
	};
	// Tri-Linear Interoplation of Grid Cell Indices ijk(0|+1) with Coefficents (t,s,r). 
	auto interpLin_scalar = [&](int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) -> float
	{
		float L_000_001_t = lerp(grid_0->getdata(i0, j0, k0), grid_0->getdata(i0, j0, k1), t1);
		float L_010_011_t = lerp(grid_0->getdata(i0, j1, k0), grid_0->getdata(i0, j1, k1), t1);
		float L_100_101_s = lerp(grid_0->getdata(i1, j0, k0), grid_0->getdata(i1, j0, k1), t1);
		float L_110_111_t = lerp(grid_0->getdata(i1, j1, k0), grid_0->getdata(i1, j1, k1), t1);
		float L_A = lerp(L_000_001_t, L_010_011_t, s1);
		float L_B = lerp(L_100_101_s, L_110_111_t, s1);
		float L_F = lerp(L_A, L_B, r1);

		return L_F;
	};

	float csize = (1.0f / (float) N_dim) * 0.5f; // Spacing / Size of Single Cell in 0.0-1.0 GS.

	for (int k = 1; k <= N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// CURRENT phi0
				float u = f3obj->vel->getdata_x(i, j, k); float v = f3obj->vel->getdata_y(i, j, k); float w = f3obj->vel->getdata_z(i, j, k);
				vec3<float> phi0_gs = idx_indexToGrid(i, j, k, N_dim);

				// Do SL Advection In Grid Space (BackTrace along vel) 
				float xf = phi0_gs.x - dt * u;
				xf = std::max(csize, std::min(xf, 1.0f + csize)); // Gridspace Clamp. 
				float yf = phi0_gs.y - dt * v;
				yf = std::max(csize, std::min(yf, 1.0f + csize)); 
				float zf = phi0_gs.z - dt * w;
				zf = std::max(csize, std::min(zf, 1.0f + csize)); 

				// Get BackTraced Cell (0|+1) Indices and resulting Coefficents - 
				vec3<float> phi1_idx = idx_gridToIndex(xf, yf, zf, N_dim); // IS.
				// Traced Cell and +1 Indices (Floored) -
				int i0f = (int) phi1_idx.x; int i1f = i0f + 1;
				int j0f = (int) phi1_idx.y; int j1f = j0f + 1;
				int k0f = (int) phi1_idx.z; int k1f = k0f + 1;
				// Traced Interoplation Coefficents - 
				float r1f = phi1_idx.x - (float) i0f; 
				float s1f = phi1_idx.y - (float) j0f;
				float t1f = phi1_idx.z - (float) k0f;

				// Trilinear or TriCosine Interpolation of BackTrace Location. 
				float phi1; 
				if (Parms.p_InteroplationType == Parms.Interoplation_Linear)
				{
					phi1 = interpLin_scalar(i0f, i1f, j0f, j1f, k0f, k1f, t1f, s1f, r1f);
				}
				else if (Parms.p_InteroplationType == Parms.Interoplation_Cosine)
				{
					phi1 = interpCos_scalar(i0f, i1f, j0f, j1f, k0f, k1f, t1f, s1f, r1f);
				}

				grid_1->setdata(phi1, i, j, k);
			}
		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	#if doedgebound == 1
	edge_bounds(grid_1); // Generic Func Call, Pass Grid_1 (n+1). 
	#endif

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}

// Semi Lagrangian Advection - Single Step (Euler) - Vector Overload - 
void fluidsolver_3::advect_sl(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1)
{
	// Interoplate Vector Components Quanitity At Traced Cell Indices (0,1) with Coefficents (t,s,r) - 
	auto cosvec = [&](int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) -> vec3<float>
	{
		// Interoplate Neighbours - for Velocity comp (U/x). 
		float U_000_001_t = cosinterp(grid_0->getdata_x(i0, j0, k0), grid_0->getdata_x(i0, j0, k1), t1);
		float U_010_011_t = cosinterp(grid_0->getdata_x(i0, j1, k0), grid_0->getdata_x(i0, j1, k1), t1);
		float U_100_101_s = cosinterp(grid_0->getdata_x(i1, j0, k0), grid_0->getdata_x(i1, j0, k1), t1);
		float U_110_111_t = cosinterp(grid_0->getdata_x(i1, j1, k0), grid_0->getdata_x(i1, j1, k1), t1);
		float U_A = cosinterp(U_000_001_t, U_010_011_t, s1);
		float U_B = cosinterp(U_100_101_s, U_110_111_t, s1);
		float U_F = cosinterp(U_A, U_B, r1);
		// Interoplate Neighbours - for Velocity comp (V/y). 
		float V_000_001_t = cosinterp(grid_0->getdata_y(i0, j0, k0), grid_0->getdata_y(i0, j0, k1), t1);
		float V_010_011_t = cosinterp(grid_0->getdata_y(i0, j1, k0), grid_0->getdata_y(i0, j1, k1), t1);
		float V_100_101_s = cosinterp(grid_0->getdata_y(i1, j0, k0), grid_0->getdata_y(i1, j0, k1), t1);
		float V_110_111_t = cosinterp(grid_0->getdata_y(i1, j1, k0), grid_0->getdata_y(i1, j1, k1), t1);
		float V_A = cosinterp(V_000_001_t, V_010_011_t, s1);
		float V_B = cosinterp(V_100_101_s, V_110_111_t, s1);
		float V_F = cosinterp(V_A, V_B, r1);
		// Interoplate Neighbours - for Velocity comp (W/z). 
		float W_000_001_t = cosinterp(grid_0->getdata_z(i0, j0, k0), grid_0->getdata_z(i0, j0, k1), t1);
		float W_010_011_t = cosinterp(grid_0->getdata_z(i0, j1, k0), grid_0->getdata_z(i0, j1, k1), t1);
		float W_100_101_s = cosinterp(grid_0->getdata_z(i1, j0, k0), grid_0->getdata_z(i1, j0, k1), t1);
		float W_110_111_t = cosinterp(grid_0->getdata_z(i1, j1, k0), grid_0->getdata_z(i1, j1, k1), t1);
		float W_A = cosinterp(W_000_001_t, W_010_011_t, s1);
		float W_B = cosinterp(W_100_101_s, W_110_111_t, s1);
		float W_F = cosinterp(W_A, W_B, r1);

		return vec3<float>(U_F, V_F, W_F);
	};

	// Return Min/Max Vector of Cells at 3D (0|+1) Interoplation Indices (For Limiter). Needs Optimizing (SIMD). 
	auto minmax = [&](int i0, int i1, int j0, int j1, int k0, int k1) -> std::tuple<vec3<float>, vec3<float>>
	{
		// Cell Indices ijk(0|+1) to get Min/Max of (For MacCormack Limiter of Quanitiy (Scalar))
		float C_000_x = grid_0->getdata_x(i0, j0, k0); float C_001_x = grid_0->getdata_x(i0, j0, k1);
		float C_010_x = grid_0->getdata_x(i0, j1, k1); float C_011_x = grid_0->getdata_x(i0, j1, k0);
		float C_100_x = grid_0->getdata_x(i1, j0, k0); float C_101_x = grid_0->getdata_x(i1, j0, k1);
		float C_110_x = grid_0->getdata_x(i1, j1, k1); float C_111_x = grid_0->getdata_x(i1, j1, k0);
		float C_000_y = grid_0->getdata_y(i0, j0, k0); float C_001_y = grid_0->getdata_y(i0, j0, k1);
		float C_010_y = grid_0->getdata_y(i0, j1, k1); float C_011_y = grid_0->getdata_y(i0, j1, k0);
		float C_100_y = grid_0->getdata_y(i1, j0, k0); float C_101_y = grid_0->getdata_y(i1, j0, k1);
		float C_110_y = grid_0->getdata_y(i1, j1, k1); float C_111_y = grid_0->getdata_y(i1, j1, k0);
		float C_000_z = grid_0->getdata_z(i0, j0, k0); float C_001_z = grid_0->getdata_z(i0, j0, k1);
		float C_010_z = grid_0->getdata_z(i0, j1, k1); float C_011_z = grid_0->getdata_z(i0, j1, k0);
		float C_100_z = grid_0->getdata_z(i1, j0, k0); float C_101_z = grid_0->getdata_z(i1, j0, k1);
		float C_110_z = grid_0->getdata_z(i1, j1, k1); float C_111_z = grid_0->getdata_z(i1, j1, k0);

		float f_min_x = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_x, C_001_x), C_010_x), C_011_x), C_100_x), C_101_x), C_110_x), C_111_x);
		float f_max_x = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_x, C_001_x), C_010_x), C_011_x), C_100_x), C_101_x), C_110_x), C_111_x);
		float f_min_y = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_y, C_001_y), C_010_y), C_011_y), C_100_y), C_101_y), C_110_y), C_111_y);
		float f_max_y = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_y, C_001_y), C_010_y), C_011_y), C_100_y), C_101_y), C_110_y), C_111_y);
		float f_min_z = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_z, C_001_z), C_010_z), C_011_z), C_100_z), C_101_z), C_110_z), C_111_z);
		float f_max_z = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_z, C_001_z), C_010_z), C_011_z), C_100_z), C_101_z), C_110_z), C_111_z);

		vec3<float> v_min(f_min_x, f_min_y, f_min_z);
		vec3<float> v_max(f_max_x, f_max_y, f_max_z);
		return std::make_tuple(v_min, v_max);
	};

	float csize = (1.0f / (float)N_dim) * 0.5f; // Spacing / Size of Single Cell in 0.0-1.0 GS.

	#pragma omp parallel for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp for
			for (int i = 1; i <= N_dim; i++)
			{
				// Grid0(prev) Current Cell (i,j,k) (phi0)
				float u = f3obj->vel->getdata_x(i, j, k); float v = f3obj->vel->getdata_y(i, j, k); float w = f3obj->vel->getdata_z(i, j, k);
				vec3<float> phi0_gs = idx_indexToGrid(i, j, k, N_dim);

				// Advection BackTrace - 
				float xf = phi0_gs.x - dt * u;
				xf = std::max(csize, std::min(xf, 1.0f + csize)); // GridSpace Clamp. 
				float yf = phi0_gs.y - dt * v;
				yf = std::max(csize, std::min(yf, 1.0f + csize));
				float zf = phi0_gs.z - dt * w;
				zf = std::max(csize, std::min(zf, 1.0f + csize));

				// Traced GridSpace Coords -> IDX Space (0,+1) Indices - 
				vec3<float> phi1prime_idx = idx_gridToIndex(xf, yf, zf, N_dim);
				int i0f = int(phi1prime_idx.x); int i1f = i0f + 1;
				int j0f = int(phi1prime_idx.y); int j1f = j0f + 1;
				int k0f = int(phi1prime_idx.z); int k1f = k0f + 1;

				// Traced Interoplation Coefficents Between Cell(0,+1) - 
				float r1f = phi1prime_idx.x - i0f;
				float s1f = phi1prime_idx.y - j0f;
				float t1f = phi1prime_idx.z - k0f;

				// Sample at BackTraced Location with Tri-Cosine Interoplation.
				vec3<float> phi1 = cosvec(i0f, i1f, j0f, j1f, k0f, k1f, t1f, s1f, r1f);

				grid_1->setdata(phi1, i, j, k);
			}
		}
	}

	// Call Boundary Condtions Post Advection (Vector) - 
	#if doedgebound == 1
	edge_bounds(grid_1); 
	#endif
	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}


// MACCORMACK - (EULER) ADVECTION \\

// MacCormack (Euler Per Trace (Forwards,Backwards)) Advection - SCALAR FIELD Overload - 
void fluidsolver_3::advect_mc(grid3_scalar<float> *grid_0, grid3_scalar<float> *grid_1)
{
	// Tri-Cosine Interoplate Scalar Quanitity At Traced Cell Indices (0,1) Of Grid0 (PrevFrame Grid) with Coefficents (t,s,r) - 
	auto interpCos_scalar = [&](int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) -> float
	{
		float C_000_001_t = cosinterp(grid_0->getdata(i0, j0, k0), grid_0->getdata(i0, j0, k1), t1);
		float C_010_011_t = cosinterp(grid_0->getdata(i0, j1, k0), grid_0->getdata(i0, j1, k1), t1);
		float C_100_101_t = cosinterp(grid_0->getdata(i1, j0, k0), grid_0->getdata(i1, j0, k1), t1);
		float C_110_111_t = cosinterp(grid_0->getdata(i1, j1, k0), grid_0->getdata(i1, j1, k1), t1);
		float C_A = cosinterp(C_000_001_t, C_010_011_t, s1);
		float C_B = cosinterp(C_100_101_t, C_110_111_t, s1);
		float C_F = cosinterp(C_A, C_B, r1);

		return C_F;
	};

	// Sample + Interoplate Vector Components Quanitity At Traced Cell Indices (0,1) Of Grid0 (PrevFrame Grid) with Coefficents (t,s,r) - 
	auto interpCos_vector = [&](int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) -> vec3<float>
	{
		// Interoplate Neighbours - for Velocity comp (U/x). 
		float U_000_001_t = cosinterp(f3obj->prev_vel->getdata_x(i0, j0, k0), f3obj->prev_vel->getdata_x(i0, j0, k1), t1);
		float U_010_011_t = cosinterp(f3obj->prev_vel->getdata_x(i0, j1, k0), f3obj->prev_vel->getdata_x(i0, j1, k1), t1);
		float U_100_101_s = cosinterp(f3obj->prev_vel->getdata_x(i1, j0, k0), f3obj->prev_vel->getdata_x(i1, j0, k1), t1);
		float U_110_111_t = cosinterp(f3obj->prev_vel->getdata_x(i1, j1, k0), f3obj->prev_vel->getdata_x(i1, j1, k1), t1);
		float U_A = cosinterp(U_000_001_t, U_010_011_t, s1);
		float U_B = cosinterp(U_100_101_s, U_110_111_t, s1);
		float U_F = cosinterp(U_A, U_B, r1);
		// Interoplate Neighbours - for Velocity comp (V/y). 
		float V_000_001_t = cosinterp(f3obj->prev_vel->getdata_y(i0, j0, k0), f3obj->prev_vel->getdata_y(i0, j0, k1), t1);
		float V_010_011_t = cosinterp(f3obj->prev_vel->getdata_y(i0, j1, k0), f3obj->prev_vel->getdata_y(i0, j1, k1), t1);
		float V_100_101_s = cosinterp(f3obj->prev_vel->getdata_y(i1, j0, k0), f3obj->prev_vel->getdata_y(i1, j0, k1), t1);
		float V_110_111_t = cosinterp(f3obj->prev_vel->getdata_y(i1, j1, k0), f3obj->prev_vel->getdata_y(i1, j1, k1), t1);
		float V_A = cosinterp(V_000_001_t, V_010_011_t, s1);
		float V_B = cosinterp(V_100_101_s, V_110_111_t, s1);
		float V_F = cosinterp(V_A, V_B, r1);
		// Interoplate Neighbours - for Velocity comp (W/z). 
		float W_000_001_t = cosinterp(f3obj->prev_vel->getdata_z(i0, j0, k0), f3obj->prev_vel->getdata_z(i0, j0, k1), t1);
		float W_010_011_t = cosinterp(f3obj->prev_vel->getdata_z(i0, j1, k0), f3obj->prev_vel->getdata_z(i0, j1, k1), t1);
		float W_100_101_s = cosinterp(f3obj->prev_vel->getdata_z(i1, j0, k0), f3obj->prev_vel->getdata_z(i1, j0, k1), t1);
		float W_110_111_t = cosinterp(f3obj->prev_vel->getdata_z(i1, j1, k0), f3obj->prev_vel->getdata_z(i1, j1, k1), t1);
		float W_A = cosinterp(W_000_001_t, W_010_011_t, s1);
		float W_B = cosinterp(W_100_101_s, W_110_111_t, s1);
		float W_F = cosinterp(W_A, W_B, r1);

		return vec3<float>(U_F, V_F, W_F);
	};

	// Return Min/Max (Scalar) of Cells at 3D (0|+1) Indices. 
	auto minmax = [&](int i0, int i1, int j0, int j1, int k0, int k1) -> std::tuple<float, float>
	{
		// Cell Indices ijk(0|+1) to get Min/Max of (For MacCormack Limiter of Quanitiy (Scalar))
		float C_000 = grid_0->getdata(i0, j0, k0); float C_001 = grid_0->getdata(i0, j0, k1);
		float C_010 = grid_0->getdata(i0, j1, k1); float C_011 = grid_0->getdata(i0, j1, k0);
		float C_100 = grid_0->getdata(i1, j0, k0); float C_101 = grid_0->getdata(i1, j0, k1);
		float C_110 = grid_0->getdata(i1, j1, k1); float C_111 = grid_0->getdata(i1, j1, k0);

		float f_min = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000, C_001), C_010), C_011), C_100), C_101), C_110), C_111);
		float f_max = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000, C_001), C_010), C_011), C_100), C_101), C_110), C_111);

		return std::make_tuple(f_min, f_max);
	};

	float csize = (1.0f / (float)N_dim) * 0.5f; // Spacing / Size of Single Cell in 0.0-1.0 GS.
	
	for (int k = 1; k <= N_dim; k++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				float phi0, phi0prime, phi1, phi1prime;

				// Grid0(prev) Current Cell (i,j,k) (phi0)
				float u = f3obj->vel->getdata_x(i, j, k); float v = f3obj->vel->getdata_y(i, j, k); float w = f3obj->vel->getdata_z(i, j, k);
				phi0 = grid_0->getdata(i, j, k);
				vec3<float> phi0_gs = idx_indexToGrid(i, j, k, N_dim); // To Grid Space. 

				// FORWARDS TRACE (phi1prime) 
				float dt_f = dt; // dt Time Forwards. 
				float xf = phi0_gs.x - dt_f * u;
				xf = std::max(csize, std::min(xf, 1.0f + csize)); // GridSpace Clamp. 
				float yf = phi0_gs.y - dt_f * v;
				yf = std::max(csize, std::min(yf, 1.0f + csize)); 
				float zf = phi0_gs.z - dt_f * w;
				zf = std::max(csize, std::min(zf, 1.0f + csize)); 

				// Forward Trace Grid Space to Index Space - 
				vec3<float> phi1prime_idx = idx_gridToIndex(xf, yf, zf, N_dim); // IS.
				// Forwad Trace Cell and +1 Indices (Floored). 
				int i0f = int(phi1prime_idx.x); int i1f = i0f + 1;
				int j0f = int(phi1prime_idx.y); int j1f = j0f + 1;
				int k0f = int(phi1prime_idx.z); int k1f = k0f + 1;
				// Forward Trace Interoplation Coefficents - 
				float r1f = phi1prime_idx.x - i0f;
				float s1f = phi1prime_idx.y - j0f;
				float t1f = phi1prime_idx.z - k0f;

				// Get Scalar Quanitiy at Forward Traced Location (phi1prime)
				phi1prime = interpCos_scalar(i0f, i1f, j0f, j1f, k0f, k1f, t1f, s1f, r1f);
				// Get Velocity (prev_vel) at Forward Traced Location (vel_f) 
				vec3<float> vel_f = interpCos_vector(i0f, i1f, j0f, j1f, k0f, k1f, t1f, s1f, r1f);

				// BACKWARDS TRACE (phi0prime) \\
				// (From Forward Traced Resulting "Particle" Location, along ForwardsTrace Sampled Vel, by -dt) 
				float dt_b = -dt; // dt Time Backwards. 
				float xb = xf - dt_b * vel_f.x;
				xb = std::max(csize, std::min(xb, 1.0f + csize)); // GridSpace Clamp. 
				float yb = yf - dt_b * vel_f.y;
				yb = std::max(csize, std::min(yb, 1.0f + csize));
				float zb = zf - dt_b * vel_f.z;
				zb = std::max(csize, std::min(zb, 1.0f + csize));

				// Backward Trace Grid Space to Index Space - 
				vec3<float> phi0prime_idx = idx_gridToIndex(xb, yb, zb, N_dim); // IS
				// Backward Trace Cell and +1 Indices. (Floored)
				int i0b = int(phi0prime_idx.x); int i1b = i0b + 1;
				int j0b = int(phi0prime_idx.y); int j1b = j0b + 1;
				int k0b = int(phi0prime_idx.z); int k1b = k0b + 1;
				// Backward Trace Interoplation Coefficents - 
				float r1b = phi0prime_idx.x - i0b;
				float s1b = phi0prime_idx.y - j0b;
				float t1b = phi0prime_idx.z - k0b;

				// Get Backwards Traced time Value. (phi0prime)
				phi0prime = interpCos_scalar(i0b, i1b, j0b, j1b, k0b, k1b, t1b, s1b, r1b);

				// Solve for error, UnClamped resulting McC value. 
				phi1prime = phi1prime + (phi0 - phi0prime) / 2.0f; 

				// Clamp Advected Quanity to min,max of neigbouring traced cells (limiter) -
				std::tuple<float, float> f_minmax = minmax(i0f, i1f, j0f, j1f, k0f, k1f); // Min/Max of ForwardTrace Loc Cells. 
				// Clamp Using Limiter (Min,Max)- 
				phi1prime = std::max(std::get<0>(f_minmax), std::min(phi1prime, std::get<1>(f_minmax)));

				grid_1->setdata(phi1prime, i, j, k);
			}
		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	#if doedgebound == 1
	edge_bounds(grid_1); // Generic Func Call, Pass Grid_1 (n+1). 
	#endif

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}

// MacCormack (Euler Per Trace (Forwards,Backwards)) Advection - VECTOR FIELD Overload - 
void fluidsolver_3::advect_mc(grid3_vector<vec3<float>> *grid_0, grid3_vector<vec3<float>> *grid_1)
{
	// Interoplate Vector Components Quanitity At Traced Cell Indices (0,1) with Coefficents (t,s,r) - 
	auto cosvec = [&](int i0, int i1, int j0, int j1, int k0, int k1, float t1, float s1, float r1) -> vec3<float>
	{
		// Interoplate Neighbours - for Velocity comp (U/x). 
		float U_000_001_t = cosinterp(grid_0->getdata_x(i0, j0, k0), grid_0->getdata_x(i0, j0, k1), t1);
		float U_010_011_t = cosinterp(grid_0->getdata_x(i0, j1, k0), grid_0->getdata_x(i0, j1, k1), t1);
		float U_100_101_s = cosinterp(grid_0->getdata_x(i1, j0, k0), grid_0->getdata_x(i1, j0, k1), t1);
		float U_110_111_t = cosinterp(grid_0->getdata_x(i1, j1, k0), grid_0->getdata_x(i1, j1, k1), t1);
		float U_A = cosinterp(U_000_001_t, U_010_011_t, s1);
		float U_B = cosinterp(U_100_101_s, U_110_111_t, s1);
		float U_F = cosinterp(U_A, U_B, r1);
		// Interoplate Neighbours - for Velocity comp (V/y). 
		float V_000_001_t = cosinterp(grid_0->getdata_y(i0, j0, k0), grid_0->getdata_y(i0, j0, k1), t1);
		float V_010_011_t = cosinterp(grid_0->getdata_y(i0, j1, k0), grid_0->getdata_y(i0, j1, k1), t1);
		float V_100_101_s = cosinterp(grid_0->getdata_y(i1, j0, k0), grid_0->getdata_y(i1, j0, k1), t1);
		float V_110_111_t = cosinterp(grid_0->getdata_y(i1, j1, k0), grid_0->getdata_y(i1, j1, k1), t1);
		float V_A = cosinterp(V_000_001_t, V_010_011_t, s1);
		float V_B = cosinterp(V_100_101_s, V_110_111_t, s1);
		float V_F = cosinterp(V_A, V_B, r1);
		// Interoplate Neighbours - for Velocity comp (W/z). 
		float W_000_001_t = cosinterp(grid_0->getdata_z(i0, j0, k0), grid_0->getdata_z(i0, j0, k1), t1);
		float W_010_011_t = cosinterp(grid_0->getdata_z(i0, j1, k0), grid_0->getdata_z(i0, j1, k1), t1);
		float W_100_101_s = cosinterp(grid_0->getdata_z(i1, j0, k0), grid_0->getdata_z(i1, j0, k1), t1);
		float W_110_111_t = cosinterp(grid_0->getdata_z(i1, j1, k0), grid_0->getdata_z(i1, j1, k1), t1);
		float W_A = cosinterp(W_000_001_t, W_010_011_t, s1);
		float W_B = cosinterp(W_100_101_s, W_110_111_t, s1);
		float W_F = cosinterp(W_A, W_B, r1);

		return vec3<float>(U_F, V_F, W_F);
	};

	// Return Min/Max Vector of Cells at 3D (0|+1) Interoplation Indices (For Limiter). Needs Optimizing (SIMD). 
	auto minmax = [&](int i0, int i1, int j0, int j1, int k0, int k1) -> std::tuple<vec3<float>, vec3<float>>
	{
		// Cell Indices ijk(0|+1) to get Min/Max of (For MacCormack Limiter of Quanitiy (Scalar))
		float C_000_x = grid_0->getdata_x(i0, j0, k0); float C_001_x = grid_0->getdata_x(i0, j0, k1);
		float C_010_x = grid_0->getdata_x(i0, j1, k1); float C_011_x = grid_0->getdata_x(i0, j1, k0);
		float C_100_x = grid_0->getdata_x(i1, j0, k0); float C_101_x = grid_0->getdata_x(i1, j0, k1);
		float C_110_x = grid_0->getdata_x(i1, j1, k1); float C_111_x = grid_0->getdata_x(i1, j1, k0);
		float C_000_y = grid_0->getdata_y(i0, j0, k0); float C_001_y = grid_0->getdata_y(i0, j0, k1);
		float C_010_y = grid_0->getdata_y(i0, j1, k1); float C_011_y = grid_0->getdata_y(i0, j1, k0);
		float C_100_y = grid_0->getdata_y(i1, j0, k0); float C_101_y = grid_0->getdata_y(i1, j0, k1);
		float C_110_y = grid_0->getdata_y(i1, j1, k1); float C_111_y = grid_0->getdata_y(i1, j1, k0);
		float C_000_z = grid_0->getdata_z(i0, j0, k0); float C_001_z = grid_0->getdata_z(i0, j0, k1);
		float C_010_z = grid_0->getdata_z(i0, j1, k1); float C_011_z = grid_0->getdata_z(i0, j1, k0);
		float C_100_z = grid_0->getdata_z(i1, j0, k0); float C_101_z = grid_0->getdata_z(i1, j0, k1);
		float C_110_z = grid_0->getdata_z(i1, j1, k1); float C_111_z = grid_0->getdata_z(i1, j1, k0);

		float f_min_x = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_x, C_001_x), C_010_x), C_011_x), C_100_x), C_101_x), C_110_x), C_111_x);
		float f_max_x = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_x, C_001_x), C_010_x), C_011_x), C_100_x), C_101_x), C_110_x), C_111_x);
		float f_min_y = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_y, C_001_y), C_010_y), C_011_y), C_100_y), C_101_y), C_110_y), C_111_y);
		float f_max_y = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_y, C_001_y), C_010_y), C_011_y), C_100_y), C_101_y), C_110_y), C_111_y);
		float f_min_z = std::min(std::min(std::min(std::min(std::min(std::min(std::min(C_000_z, C_001_z), C_010_z), C_011_z), C_100_z), C_101_z), C_110_z), C_111_z);
		float f_max_z = std::max(std::max(std::max(std::max(std::max(std::max(std::max(C_000_z, C_001_z), C_010_z), C_011_z), C_100_z), C_101_z), C_110_z), C_111_z);

		vec3<float> v_min(f_min_x, f_min_y, f_min_z);
		vec3<float> v_max(f_max_x, f_max_y, f_max_z);
		return std::make_tuple(v_min, v_max);
	};

	float csize = (1.0f / (float)N_dim) * 0.5f; // Spacing / Size of Single Cell in 0.0-1.0 GS.

	#pragma omp parallel for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp for
			for (int i = 1; i <= N_dim; i++)
			{
				vec3<float> phi0, phi0prime, phi1, phi1prime;

				// Grid0(prev) Current Cell (i,j,k) (phi0)
				float u = f3obj->vel->getdata_x(i, j, k); float v = f3obj->vel->getdata_y(i, j, k); float w = f3obj->vel->getdata_z(i, j, k);
				phi0 = grid_0->getdata(i, j, k);
				vec3<float> phi0_gs = idx_indexToGrid(i, j, k, N_dim);

				// FORWARDS TRACE (Along Postive dt, along curcell (phi0) sampled velocity (phi1prime) - 
				float xf = phi0_gs.x - dt * u;
				xf = std::max(csize, std::min(xf, 1.0f + csize)); // GridSpace Clamp. 
				float yf = phi0_gs.y - dt * v;
				yf = std::max(csize, std::min(yf, 1.0f + csize));
				float zf = phi0_gs.z - dt * w;
				zf = std::max(csize, std::min(zf, 1.0f + csize));

				// Forward-Traced GridSpace Coords -> IDX Space (0,+1). 
				vec3<float> phi1prime_idx = idx_gridToIndex(xf, yf, zf, N_dim);
				int i0f = int(phi1prime_idx.x); int i1f = i0f + 1;
				int j0f = int(phi1prime_idx.y); int j1f = j0f + 1;
				int k0f = int(phi1prime_idx.z); int k1f = k0f + 1;

				// Interoplation Coefficents ForwardTrace Between Cell(0,+1) - 
				float r1f = phi1prime_idx.x - i0f;
				float s1f = phi1prime_idx.y - j0f;
				float t1f = phi1prime_idx.z - k0f;

				// Get Forward Traced time Vector Value. (phi1prime)
				phi1prime = cosvec(i0f, i1f, j0f, j1f, k0f, k1f, t1f, s1f, r1f);

				// BACKWARDS TRACE (Along Negative dt, by Forward Trace Sampled Vel) (phi0prime) - 
				float xb = xf - (-dt) * phi1prime.x;
				xb = std::max(csize, std::min(xb, 1.0f + csize)); // GridSpace Clamp. 
				float yb = yf - (-dt) * phi1prime.y;
				yb = std::max(csize, std::min(yb, 1.0f + csize));
				float zb = zf - (-dt) * phi1prime.z;
				zb = std::max(csize, std::min(zb, 1.0f + csize));

				// Back-Traced GridSpace Coords -> IDX Space (0,+1) - 
				vec3<float> phi0prime_idx = idx_gridToIndex(xb, yb, zb, N_dim);
				int i0b = int(phi0prime_idx.x); int i1b = i0b + 1;
				int j0b = int(phi0prime_idx.y); int j1b = j0b + 1;
				int k0b = int(phi0prime_idx.z); int k1b = k0b + 1;
				// Interoplation Coefficents BackTrace Between Cell(0,+1) - 
				float r1b = phi0prime_idx.x - i0b;
				float s1b = phi0prime_idx.y - j0b;
				float t1b = phi0prime_idx.z - k0b;

				// Get Backwards Traced time Value (phi0prime).
				phi0prime = cosvec(i0b, i1b, j0b, j1b, k0b, k1b, t1b, s1b, r1b);

				// Solve for Error
				phi1prime = phi1prime + (phi0 - phi0prime) / 2.0f; // UnClamped McC Vector Value 

				// Limiter (Clamp to Forward Traced Min/Max Vector Values) - 
				vec3<float> phi1prime_uc, phi1prime_c;
				std::tuple<vec3<float>, vec3<float>> v_minmax = minmax(i0f, i1f, j0f, j1f, k0f, k1f);
				vec3<float> v_min = std::get<0>(v_minmax); vec3<float> v_max = std::get<1>(v_minmax);
				// Un-Clamped and Clamped Values for Limiter Strength. 
				phi1prime_uc = phi1prime;
				phi1prime_c = vec_clamp(std::get<0>(v_minmax), std::get<1>(v_minmax), phi1prime);
				// Lerp Limiter Strength 
				phi1prime = vec_lerp(phi1prime_uc, phi1prime_c, Parms.p_McC_LimiterStrength);

				// Set Final Cur_Grid Advected Vector Cell Value. 
				grid_1->setdata(phi1prime, i, j, k);
			}
		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	#if doedgebound == 1
	edge_bounds(grid_1); // Generic Func Call, Pass Grid_1 (n+1). 
	#endif

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}

/*	====================================================
	VORTICITY CONFINEMENT
==================================================== */





/* ====================================================
	VELOCITY PROJECTION / PRESSURE SOLVE 	
	==================================================== */

/* Velocity Projection / Pressure Solve 
   Calculate Pressure Field, whose Gradient when subtraced from velocity field, will result in Velocity Divergence ~0.
   LinearSystem of Matrixless Laplacian A to solve uknown Pressure Vector x by caluclated diverence vector b. */ 

// PROJECTION - GAUSS-SEIDEL + SUCESSIVE OVER RELAXATION (GS+SOR) - Single-Threaded Pressure Solve. 
void fluidsolver_3::project(int iter)
{
	float h = 1.0f/N_dim; // CellSize Assuming Cube'd Grid of Length (N_dim cells). 

	// Ensure previous Divergence and Pressure Temp Grids have been deleted before new grid ptr assignement/Alloc. 
	del_divergence(); del_pressure(); 

	// Init Solvers Own Divergence and Pressure Fields, Only needed for Scope of this Function.  
	divergence = new grid3_scalar<float>(x_s, y_s, z_s, e_s);
	pressure = new grid3_scalar<float>(x_s, y_s, z_s, e_s);

	// DIVERGENCE FIELD CALC \\ -

	// Thread Safe, Per Cell Operations. 
	#pragma omp for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp for
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
					+ f3obj->vel->getdata_z(i, j,k + 1) - f3obj->vel->getdata_z(i, j, k - 1));
 
				divergence->setdata(div, i, j, k);
				pressure->setdata(0.0f, i, j, k); // Init 0 guess. 
			}
		}
	}
	// Call Boundary Condtions on Divergence Field and Inital Pressure Field - 
	#if doedgebound == 1
	edge_bounds(divergence); edge_bounds(pressure);
	#endif
	#if dospherebound == 1
	sphere_bounds_eval(divergence, spherebound_coliso);
	sphere_bounds_eval(pressure, spherebound_coliso);
	#endif

	// PRESSURE FIELD LINEAR SOLVE \\ 

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
		#if doedgebound == 1
		edge_bounds(pressure);
		#endif
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

	// Thread Safe, Per Cell Operations. 
	#pragma omp parallel for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
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
	#if doedgebound == 1
	edge_bounds(f3obj->vel);
	#endif
	#if dospherebound == 1
	sphere_bounds_eval(f3obj->vel, spherebound_coliso);
	#endif
		
	// Delete Solver Temp Pressure and Divergence Fields - 
	del_divergence(); del_pressure(); // Faster than zeroing. 
}
// End of Velocity Projection Implementation (GS + SOR).


/* PROJECTION - JACOBI to Solve Pressure Poission Equation -
   Allows Multithreading as Cells Pressure Solve, Only lookup Previous Pressure Values, but results in slower Convergence. */

void fluidsolver_3::project_jacobi(int iter)
{
	float h = 1.0f/N_dim; // CellSize 1.0f/N (N_dim); 

	// Ensure prev D&P Temp Grids have been deleted.
	del_divergence(); del_pressure();
	
	// Alloc New DP Grids. 
	divergence = new grid3_scalar<float>(x_s, y_s, z_s, e_s);
	pressure = new grid3_scalar<float>(x_s, y_s, z_s, e_s);
	pressure_1 = new grid3_scalar<float>(x_s, y_s, z_s, e_s);

	// DIVERGENCE FIELD CALC \\ - 

	#pragma omp parallel for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				divergence->setdata(0.0f, i, j, k);

				// Compute Divergence Cell Value. 
				float div = -0.5 * h * (f3obj->vel->getdata_x(i + 1, j, k) - f3obj->vel->getdata_x(i - 1, j, k)
					+ f3obj->vel->getdata_y(i, j + 1, k) - f3obj->vel->getdata_y(i, j - 1, k)
					+ f3obj->vel->getdata_z(i, j, k + 1) - f3obj->vel->getdata_z(i, j, k - 1));

				divergence->setdata(div, i, j, k);
				pressure->setdata(0.0f, i, j, k); // Init 0 Guess. 
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

		edge_bounds(pressure_1);
		#if dospherebound == 1
		sphere_bounds_eval(pressure_1,spherebound_coliso); 
		#endif	

		// Swap Pressure Grid with Scratch Grid at end of k Iter After internal i,j,k MT'd Jacobi Projection Iteration is complete. ST
		pressure_1->swap(pressure);
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

	del_divergence(); del_pressure();
	
}
// End of Velocity Projection (Jacobi) Implementation.



/* ====================================================
	DESNITY SOLVE STEP - 
	==================================================== */
// Implementation of Simulation Solve step of DENSITY solve operations. 

void fluidsolver_3::density_step()
{
	// DIFFUSE Density - 
	// If Using Diffusion, If Not Need to Manually set Cur to Prev, rather than swapping.
	if (Parms.p_Do_Dens_Diff == true)
	{
		f3obj->dens->swap(f3obj->prev_dens); 
		// Gauss-Seidel Iterative Density (Scalar) Diffusion - 
		diffuse(f3obj->prev_dens, f3obj->dens, Parms.p_Dens_Diffuse_Str, Parms.p_Dens_Diff_iter);
		f3obj->dens->swap(f3obj->prev_dens); 
	}
	else 
	{
		// Use "SetCurToPrev" Funcs to Copy Cur Grid to Previous Grid at start of subsolve, Oppose to via Diffusion Swap - 
		f3obj->setcurtoprev(f3obj->prev_dens, f3obj->dens);
	}

	// ADVECT Density - 

	if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_Euler)
	{
		advect_sl(f3obj->prev_dens, f3obj->dens); 
	}
	else if (Parms.p_AdvectionType == Parms.Advect_MC_Euler)
	{
		advect_mc(f3obj->prev_dens, f3obj->dens);
	}

	// DISSIPATE Density - 

	if (Parms.p_Do_Dens_Disp)
	{
		dissipate(f3obj->dens, Parms.p_Dens_Disp_Mult, this->dt); 
	}

}
// End of Density Step Implemetnation.

/*	====================================================
	Velocity SOLVE STEP -
	==================================================== */

// ! Removed Pre Advection Projection Call, So I can use single Projection call post Advection, (with higher iter count). 

void fluidsolver_3::velocity_step()
{
	// If Using Diffusion, If Not Need to Manually set Cur to Prev, rather than swapping. 
	if (Parms.p_Do_Vel_Diff == true)
	{
		f3obj->vel->swap(f3obj->prev_vel); 
		// Gauss-Seidel Iterative Vector Diffusion - 
		diffuse(f3obj->prev_vel, f3obj->vel, Parms.p_Vel_Diffuse_Str, Parms.p_Vel_Diff_iter);
		f3obj->vel->swap(f3obj->prev_vel); 
	}
	else 
	{
		f3obj->setcurtoprev(f3obj->prev_vel, f3obj->vel);
	}

	// ADVECT VELOCITY FIELD (Self Advect) \\

	if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_Euler)
	{
		advect_sl(f3obj->prev_vel, f3obj->vel);
	}
	else if (Parms.p_AdvectionType == Parms.Advect_MC_Euler)
	{
		advect_mc(f3obj->prev_vel, f3obj->vel);
	}

	// DISSIPATE Velocity
	if (Parms.p_Do_Vel_Disp)
	{
		dissipate(f3obj->vel, Parms.p_Vel_Disp_Mult, this->dt); // Density Dissipation. 
	}

	// VORTEX CONFINEMENT 
	if (Parms.p_useVorticity)
	{
		// Apply Vortex Confinement Force to Velocity Field. 
	}

	// PROJECT VELOCITY FIELD \\  (Post Advect Only) 

	if (Parms.p_ProjectionType == Parms.Project_GaussSeidel || Parms.p_ProjectionType == Parms.Project_GaussSeidel_SOR)
	{
		// GS + SOR (ST) - (SOR May or May not be enabled by FluidSolver Member Paramter). 
		project(Parms.p_GS_Proj_iter);
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

void fluidsolver_3::solve_step(bool solve, int max_step)
{
	// Init Step/time vals. 
	int step_count = 0;
	float dt_acc = 0.0f; // Accumulated DeltaTime (SimTime)
	float total_elapsed = 0.0f; // Chrono Time

	// Render Object Creation/Setup - 
	int render_mode = 1; // 0 = Density, 1 = Vel. 
	render_obj = new renderobject_3D_OGL("OpenGL", 4, 2, vec2<int>(win_size_xy, win_size_xy), vec3<int>(x_s + e_s, y_s + e_s, z_s + e_s), this->winptr); 

	// SOLVE STEP AND RENDER - EXECUTION LOOP 
	while (solve == true && step_count <= max_step) // Infinite Solve Loop For Now. Will need to link this to drawing/waiting etc. 
	{
		// Time Start - 
		std::chrono::system_clock::time_point timestep_start = std::chrono::system_clock::now();
		std::stringstream log_out;

		// STEP INPUT OPERATIONS \\ ----------------. 
		// Interactive Implicit Source Scaling Basic. 
		sphere_rad_test(); 

		// Switch Advection (DEBUGGING) 
		if (glfwGetKey(winptr, GLFW_KEY_J) == GLFW_PRESS)
		{
			Parms.p_AdvectionType = Parms.Advect_SL_BackTrace_Euler;
			log_out << "DBG::Advection Type Changed = " << Parms.AdvType_Key[Parms.Advect_SL_BackTrace_Euler] << "\n";
		}
		else if (glfwGetKey(winptr, GLFW_KEY_K) == GLFW_PRESS)
		{
			Parms.p_AdvectionType = Parms.Advect_MC_Euler;
			log_out << "DBG::Advection Type Changed = " << Parms.AdvType_Key[Parms.Advect_MC_Euler] << "\n";
		}

		// Get CurFrame Mouse Pos And Calc Updated Mouse Vel. 
		std::this_thread::sleep_for(std::chrono::milliseconds(2));
		updt_mousepos(step::STEP_CUR); updt_mouseposNorm(step::STEP_CUR); // updt_mouseposRange(step::STEP_CUR);
		updt_mousevel(); 

		// PRESTEP OPERATIONS \\ ----------------
		// Eval and Pre-Cache SphereBounds_SDF - 
		#if dospherebound
		sphere_bounds_set(spherebound_radius, spherebound_coliso, spherebound_offset);
		#endif

		// STEP SOURCING OPERATIONS \\ ----------------
		// SinWave Emitter - 
		float offs = sin(((float)step_count / (float)max_step) * 500.0f) * 0.15f;
		float offb = cos(((float)step_count / (float)max_step) * 400.0f) * 1.0; 
		offb = (offb + 1) * 0.5f; 
		f3obj->implicit_sphere_source(0.25f, vec3<float>(0.0f, 1.0f, offb), vec3<float>(offs + 0.4f, 0.1f, 0.5f), impsource_radius); // 0.01f

		// Mouse Emitter -
		/*vec3<float> emit_vel = vec3<float>(0.0f, 1.0f, 0.0f) + (vec3<float>(mouse_vel.x, mouse_vel.y, 0.0f) *= 2.0f); 
		  f3obj->implicit_sphere_source(0.5f, emit_vel, vec3<float>(xpos_1_N, 1.0f-ypos_1_N, 0.5f), impsource_radius); */

		// Static Emitter -
		//f3obj->implicit_sphere_source(0.25f, vec3<float>(0.0f, 1.0f, 0.0f), vec3<float>(0.5f, 0.1f, 0.5f), impsource_radius);

		// STEP - SUB - SOLVER STEP OPERATIONS \\ -------------- 
		velocity_step();
		density_step();

		// STEP - RENDER CALLS \\ ------------------
		render_obj->get_input(vec2<float>(xpos_1_N, ypos_1_N)); 
		render_obj->et = step_count; 
		render_obj->shader_pipe(f3obj); 
		render_obj->render_loop(rend_state::RENDER_ACTIVE); 

		// STEP DEBUG CONSLE OPERATIONS \\ -------------------- 
		std::chrono::system_clock::time_point timestep_end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed = timestep_end - timestep_start;
		double elapsed_sec = elapsed.count();
		total_elapsed += elapsed_sec;

		// Store Mouse Pos as Prev Frame Pos.. 
		updt_mousepos(step::STEP_PREV); updt_mouseposNorm(step::STEP_PREV); updt_mouseposRange(step::STEP_PREV);

		log_out << "INFO::Solve_Step::SimTime Passed = " << dt_acc << " Seconds of " << dt << " timestep \n";
		log_out << "DEBUG::Solve_Step::Step_Count =  " << step_count << "\n";
		log_out << "DEBUG::Solve_Step::WallClock_Duration = " << elapsed_sec << " Seconds" << "\n";

		if (step_count == max_step)
		{
			log_out << "\n *** INFO::Solve_Step::TOTAL_DELTA-TIME = " << dt_acc << " Seconds at timestep of " << dt << "***\n";
			log_out << " \n *** DEBUG::Solve_Step::TOTAL_WALLCLOCK_CALC_DURATION = " << total_elapsed << " Seconds *** \n";
			log_out << " \n *** DEBUG::Solve_Step::TOTAL_FPS = ~" << double(max_step) / total_elapsed << " FPS *** \n \n";
		}

		// Print Log to stdout - 
		std::cout << log_out.str() << std::endl; 

		// STEP INCREMENT \\ ------------
		dt_acc += dt; 
		step_count++; 
	}
}
// End of Solve Step Implementation.


/*	====================================================
	Utility Member Functions - 
	==================================================== */

// Window Setter. Passed from Render Context Object.  
void fluidsolver_3::set_window(GLFWwindow *win)
{
	assert(win != nullptr); 
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

		xpos_1_N /= win_size_xy; ypos_1_N /= win_size_xy;
		xpos_1_N = clamp(xpos_1_N, 0.0f, 1.0f);
		ypos_1_N = clamp(ypos_1_N, 0.0f, 1.0f);

		// Window Bounds - 
	//	if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_1_N = 0.0f, ypos_1_N = 0.0f;
	}
	else if (step_id == step::STEP_PREV)
	{
		glfwGetCursorPos(winptr, &xpos_0_N, &ypos_0_N);

		xpos_0_N /= win_size_xy; ypos_0_N /= win_size_xy;
		xpos_0_N = clamp(xpos_0_N, 0.0f, 1.0f);
		ypos_0_N = clamp(ypos_0_N, 0.0f, 1.0f);

		// Window Bounds - 
	//	if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_0_N = 0.0f, ypos_0_N = 0.0f;
	}
}

// Call To Update Mouse Pos Norm Within Clamped Range
// Not Used Currently, may be useful for window remapping later. 
void fluidsolver_3::updt_mouseposRange(const step step_id)
{
	glfwPollEvents();
	if (step_id == step::STEP_CUR)
	{
		glfwGetCursorPos(winptr, &xpos_1_R, &ypos_1_R);
		xpos_1_R /= win_size_xy; ypos_1_R /= win_size_xy;

		xpos_1_R = fitRange(xpos_1_R, 0.0f, 1.0f, -1.0f, 1.0f);
		ypos_1_R = fitRange(ypos_1_R, 0.0f, 1.0f, -1.0f, 1.0f);

		// Window Bounds - 
		if (!glfwGetWindowAttrib(winptr, GLFW_HOVERED)) xpos_1_R = 0.0f, ypos_1_R = 0.0f;
	}
	else if (step_id == step::STEP_PREV)
	{
		glfwGetCursorPos(winptr, &xpos_0_R, &ypos_0_R);
		xpos_0_R /= win_size_xy; ypos_0_R /= win_size_xy;

		xpos_0_R = fitRange(xpos_1_R, 0.0f, 1.0f, -1.0f, 1.0f);
		ypos_0_R = fitRange(ypos_1_R, 0.0f, 1.0f, -1.0f, 1.0f);

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
	disp_mult = clamp(disp_mult, 0.0f, 1.0f); 

	#pragma omp parallel for
	for (int i = 0; i <= total_size; i++)
	{
		// Don't Mult By dt for now. 
		float cur_scl = grid->getdata(i);
		cur_scl *= disp_mult;
		grid->setdata(cur_scl, i);
	}

}

void fluidsolver_3::dissipate(grid3_vector<vec3<float>> *grid, float disp_mult, float dt)
{
	disp_mult = clamp(disp_mult, 0.0f, 1.0f);
	#pragma omp parallel for
	for (int i = 1; i <= N_dim; i++)
	{
		// Don't Mult By dt for now. 
		vec3<float> cur_vel = grid->getdata(i);
		cur_vel *= disp_mult;
		grid->setdata(cur_vel, i);
	}
}

void fluidsolver_3::add_velocity(const vec3<float> &vel)
{
	#pragma omp for
	for (int k = 1; k <= N_dim; k++)
	{
		#pragma omp for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				// Only Within Emission Sphere? Add alt ImpSource with MouseVel passed from solver..
				f3obj->add_velocity(vel, i, j, k);
			}
		}
	}
}

/*	====================================================
Utility Member Functions
==================================================== */






// Test Implementation of Interactive Sphere Radius
void fluidsolver_3::sphere_rad_test()
{
	if (glfwGetKey(winptr, GLFW_KEY_A) == GLFW_PRESS)
	{
		impsource_radius += 0.001f; 
	}

	if (glfwGetKey(winptr, GLFW_KEY_D) == GLFW_PRESS)
	{
		impsource_radius -= 0.001f;
	}
}

/*	====================================================
Static Utility Member Functions
==================================================== */

// LERP,Cosinterp,FitRange,Clamp moved to header (__forceinline)
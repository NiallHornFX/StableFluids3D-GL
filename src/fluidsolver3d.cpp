// Implementation of fluidsolver_2

//#include "fluidsolver_2.h"

// C++ Headers - 
#include <memory>
#include <chrono>
#include <thread>
#include <sstream>
#include <algorithm>
#include <cassert>

// Windows Specfic Headers - 
#include <omp.h> // OpenMP 2.0

#define cur2d "[" << i << "," << j << "]  " // For Printing.
#define dospherebound 1 // Enable Sphere Collisions
#define JACOBI_DO_SOR 1 // Do Sucessive Over Relaxiation Inside Jacobi Projection.
#define DO_SPHEREBOUNDS_MT 1

#define RENDER_GL 1
#define RENDER_IMG 0

// Replace with Parms ... 
// Sphere Bound Globals (For Sphere Bound MFunc Paramters) Oppose to Pass Per Nested MFunc Call. 
vec2<float> spherebound_offset(0.0f, 1.0f);  // NOW Set to x-h (y-k) So Sign IS Correct XY. 
const float spherebound_coliso = 0.0025f; // Surface Value Threshold (0-Collison Iso Value) of Sphere SDF Func. 
float spherebound_radius = 0.003f; // Much Smaller Values now no Sqrt on ImpFunc. 0.075f;

extern short verbose; // Get Verbose Global from main. 
/* 
	TO DO -
	Refer to trello. 
*/

// !!! OPTIMIZATION IN PROGRESS - GL RENDERING DISABLED + PROJECTION OF VEL FIELD !!! \\

/* ====================================================
	FLUIDSOLVER_2 CONSTRUCTOR/DESTRUCTOR'S
   ==================================================== */

// fluidsolver_2 CONSTRUCTOR (Explict Only) - 
// Added dt Param & HardCoded Spherebound Initalization List. 

fluidsolver_2::fluidsolver_2(fluidobj_2d *f2dptr, float dtt) 
	: dt(dtt), f2obj(f2dptr) //,spherebound_coliso(0.002f), spherebound_radius(0.01f), spherebound_offset(-0.5f, 0.75f)
{
	// Get Dimensions,Spacing from FluidObject Ensure Values Match. 
	x_s = f2obj->x_s, y_s = f2obj->y_s, e_s = f2obj->e_s; 
	total_size = (std::size_t) ((x_s + e_s) * (y_s + e_s));

	N_dim =  f2dptr->x_s; // Single Dim Size, No Edge. 
	NE_dim = f2dptr->x_s + f2dptr->e_s; // Single Dim Size, With Edge Size. 

	// Init Solver (non temp) grids.
	spherebounds_sdf = new grid2_scalar(N_dim, N_dim, 2, 8, 1.0f);

}

// fluidsolver_2 DESTRUCTOR 
fluidsolver_2::~fluidsolver_2()
{
	// Check if Pressure and Divergence Grids are not set to Null
	// (Incase theyve been deleted/and ptr Nulled already). 
	// Then ofc set ptrs to Null After Delete/Deallocation of Grids. 
	if (pressure != NULL && pressure != nullptr)
	{
		delete pressure;
		pressure = nullptr;
	}

	if (divergence != NULL && divergence != nullptr)
	{
		delete divergence;
		divergence = nullptr; 
	}

	// fluidsolver_2 has no other self alloacted data.
	// WE DO NOT WANT TO DELTE fluidobj_2d PTR HERE and DELTE FLUID OBJECT...
	// Let That Be Done in its own destructor/delete call. 

	// Delete FluidObject (NOT GLFW Window, as this is pointed to elsewhere !)
	//delete render_obj; render_obj = nullptr; DEBUGGING IF CAUSING ISSUES ...
}

/* ====================================================
	MANUAL USER FIELD DELTETION
   ==================================================== */

void fluidsolver_2::del_pressure()
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

void fluidsolver_2::del_divergence()
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

// FLUIDSOLVER_2::edge_bounds INFO - 
// Set Edge Cell Reflection of XY Velocity And Average Corner Cells of Vel and Dens Field. 

// ** EDGE_BOUNDS_SCALAR-FIELD IMPLEMENTATION ** \\ - 
void fluidsolver_2::edge_bounds(grid2_scalar *grid)
{
	// bool do_mt = Parms.p_MT_Global & Parms.p_MT_Global;

	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 1; i <= N_dim; i++)
	{
		// X- Edge Boundary 
		float pe_0i = grid->getdata(1, i); // [0,i] Boundary Values from Edge [1,i] Values.
		grid->setdata(pe_0i, 0, i);

		// X+ Edge Boundary
		float pe_Ni = grid->getdata(N_dim, i); // [N+1,i] Boundary Values from Edge [N,i] Values.
		grid->setdata(pe_Ni, (N_dim + 1), i);

		// Y+ Edge Boundary
		float pe_i0 = grid->getdata(i, 1); // [i, 0] Boundary Values from Edge [i, 1] Values.
		grid->setdata(pe_i0, i, 0);

		// Y- Edge Boundary
		float pe_iN = grid->getdata(i, N_dim); // [i, N+1] Boundary Values from Edge [i,N] Values.
		grid->setdata(pe_iN, i, N_dim + 1);
	}

	// ScalarGrid Edge Bounds Corner Cells - Neighbour Averages -
	float p_00 = 0.5 * (grid->getdata(1, 0) + grid->getdata(0, 1));
	grid->setdata(p_00, 0, 0);
	float p_01 = 0.5 * (grid->getdata(1, N_dim + 1) + grid->getdata(0, N_dim));
	grid->setdata(p_01, 0, N_dim + 1);
	float p_02 = 0.5 * (grid->getdata(N_dim, 0) + grid->getdata(N_dim + 1, 1));
	f2obj->dens->setdata(p_02, N_dim + 1, 0);
	float p_03 = 0.5 * (grid->getdata(N_dim, N_dim + 1) + grid->getdata(N_dim + 1, N_dim));
	grid->setdata(p_03, N_dim + 1, N_dim + 1);
}

// ** EDGE_BOUNDS_VECTOR-FIELD IMPLEMENTATION ** \\ - 

void fluidsolver_2::edge_bounds(grid2_vector *grid)
{
	// Loop Over Grid Not Including Edge Cells, 1D based Iteration. 1-N (N_dim). Set Edge Cells Layer Along 0 & N+1 On Each Axis +\-.  
	
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 1; i <= N_dim; i++)
	{
		/*
		Assuming X Velocity is ONLY Reflected by X Axis Edge Cells, and Y Velocity is ONLY Reflected
		by Y Axis Edge Cells. This is not true for NoSlip, which would solve both velocity axis/components for each
		edge cell layer axis. But for now basic Per Edge/Bound Axis, Per Velocity Axis Reflection is implemented.
		*/

		// X Velocity Reflection - X-Axis Bound/Edge Collum Cells 

		// Vel x (u) Collum (0,(i)) -

		// Get X Vel Of Next Cell In (1, (i))
		float next_x_1 = grid->getdata_x(1, i);
		next_x_1 *= -1.0f;
		// Set Inverse Vel To EdgeCells Along vel[0, i]
		grid->setdata_x(next_x_1, 0, i);

		// Vel x (u) Collum (N+1,(i)) -

		// Get X Vel Of Next Cell In (N, i)
		float next_x_N = grid->getdata_x(N_dim, i);
		next_x_N *= -1.0f;
		// Set Inverse Vel To EdgeCells Along vel[N+1, i]
		grid->setdata_x(next_x_N, N_dim + 1, i);

		// Y Velocity Reflect - Y-Axis Bound/Edge Row Cells 

		// Vel y (v) Row (i,0) -

		// Get Y Vel Of Next Cell In ((i), 1)
		float next_y_1 = grid->getdata_y(i, 1);
		next_y_1 *= -1.0f;
		// Set Inverse Vel To EdgeCells Along vel[i, 0]
		grid->setdata_y(next_y_1, i, 0);

		// Vel y (v) Collum ((i),N+1) -

		// Get Y Vel Of Next Cell In (i, N)
		float next_y_N = grid->getdata_y(i, N_dim);
		next_y_N *= -1.0f;
		// Set Inverse Vel To EdgeCells Along vel[i, N+1]
		grid->setdata_y(next_y_N, i, N_dim + 1);

	}

		// Set Single Corner Cells (Average Of Neighbours) - (Set Per Component Oppose to using vec2 opoverload add/avg Operations)

		// X Vel Neighbour Averages - 
		float x_00 = 0.5 * (grid->getdata_x(1, 0) + grid->getdata_x(0, 1));
		grid->setdata_x(x_00, 0, 0);
		float x_01 = 0.5 * (grid->getdata_x(1, N_dim + 1) + grid->getdata_x(0, N_dim));
		grid->setdata_x(x_01, 0, N_dim + 1);
		float x_02 = 0.5 * (grid->getdata_x(N_dim, 0) + grid->getdata_x(N_dim + 1, 1));
		grid->setdata_x(x_02, N_dim + 1, 0);
		float x_03 = 0.5 * (grid->getdata_x(N_dim, N_dim + 1) + grid->getdata_x(N_dim + 1, N_dim));
		grid->setdata_x(x_03, N_dim + 1, N_dim + 1);

		// Y Vel Neighbour Averages - 
		float y_00 = 0.5 * (grid->getdata_y(1, 0) + grid->getdata_y(0, 1));
		grid->setdata_y(y_00, 0, 0);
		float y_01 = 0.5 * (grid->getdata_y(1, N_dim + 1) + grid->getdata_y(0, N_dim));
		grid->setdata_y(y_01, 0, N_dim + 1);
		float y_02 = 0.5 * (grid->getdata_y(N_dim, 0) + grid->getdata_y(N_dim + 1, 1));
		grid->setdata_y(y_02, N_dim + 1, 0);
		float y_03 = 0.5 * (grid->getdata_y(N_dim, N_dim + 1) + grid->getdata_y(N_dim + 1, N_dim));
		grid->setdata_y(y_03, N_dim + 1, N_dim + 1);
}

// ** SPHERE_BOUNDS_SCALAR-FIELD IMPLEMENTATION ** \\ - 

// Set SphereBounds At Beginning Of Each Step into FluidSolver spherebounds_sdf Scalar Grid -

void fluidsolver_2::sphere_bounds_set(float radius, float col_iso, const vec2<float> &offset)
{
	float h = 1.0f / N_dim; // Grid Spacing, Recoprical of One Dim Size (N). 

	#pragma omp parallel for num_threads(omp_get_max_threads())	
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Implicit Sphere/Circle Function: x^2 + y^2 - r . >= 0 <= Thresh == Surface. > Thresh = Exterior. < 0 = Interior. Cells. 
			vec2<float> cell_gridSpace(float((i * h) - offset.x), float((j * h) - offset.y)); // Index to Grid Space 0-1N. 
			float sphere_func = ((cell_gridSpace.x * cell_gridSpace.x) + (cell_gridSpace.y * cell_gridSpace.y)) - radius; // NoMoreSqrt.

			// Write to  spherebounds_sdf grid for spherebounds_eval() to lookup. 
			spherebounds_sdf->setdata(0.0f, i, j); // Zero Out Prev SDF Grid Cell Values. 
			spherebounds_sdf->setdata(sphere_func, i, j); // Set New Cur Step SDF Cell Values. 
		}
	}

}

// ** SPHERE_BOUNDS_EVAL - SCALAR-FIELD OVERRIDE ** \\ - 
// Eval SphereBounds - On Scalar Field. Also Set "Col" Viz Grid For Cells Inside Sphere SDF For Render col viz. 
// Unused (Surface and Exterior Condtions disabled for perf). 

void fluidsolver_2::sphere_bounds_eval(grid2_scalar *grid, float col_iso)
{
	float h = 1.0f / N_dim; // Grid Spacing, Recoprical of One Dim Size (N). 
	
	#pragma omp parallel for num_threads(omp_get_max_threads())	
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Lookup Sphere/Circle SDF At CurCell i,j. 
			float sphere_func = spherebounds_sdf->getdata(i, j);
			f2obj->col->setdata(0.0f, i, j); // Zero Out Prev Collison Cell Values. 

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
					grid->setdata(0.0f, i, j); // Zero ANY Density inside Sphere Interior. 
				}
				
				// Inside Operations - Collision Scalar Grid
				// Set Collison Grid (HardCoded to pointed f2obj member) (Col Grid, probs should be Boolean Grid no? Evneutally MAC). 
				f2obj->col->setdata(1.0f, i, j);
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


// ** SPHERE_BOUNDS_EVAL - VECTOR-FIELD OVERRIDE ** \\ - 
// Eval SphereBounds - On Vector Field. Also add Mouse Velocity to Vel Grid, Within SphereBounds. 
// (This Function does not set col grid)
// Unused (Surface and Exterior Condtions disabled for perf). 

void fluidsolver_2::sphere_bounds_eval(grid2_vector *grid, float col_iso)
{
	float h = 1.0f / N_dim; // Grid Spacing, Recoprical of One Dim Size (N). 

	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Lookup Sphere/Circle SDF At CurCell i,j. 
			float sphere_func = spherebounds_sdf->getdata(i, j);

			/*
			// Surface
			if (sphere_func >= 0.0f && sphere_func <= col_iso) // func-radis > 0.0f but <= col_iso Cell is on "surface". 
			{
				// Surface Operations - 

				// Do Sphere Vel Reflection On Surface Cirlce/Sphere Cells - 
				// Reflect from offset vector (Sphere Center) to CellPos Diriection, oppose to origin^ - 

				//float input_spd = grid->getdata(i, j).length();
				//vec2 refl_dir = vec2(offset.x, offset.y) - vec2(float(i) * h, float(j) * h);
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
				vec2 refl_dir = vec2(float(i) * h, float(j) * h) - vec2(offset.x, offset.y) ;
				refl_dir.normalize();
				refl_dir *= input_spd * 1.0f;
				// Override Interior Col Vel With Reflect Dir Vel
				grid->setdata(refl_dir, i, j); // ISSUES WITH ADDED GRID VELOCITY ? */


				// Simply Zero Out Vector Grid Values In Interior Sphere Cells -   
				grid->setdata(vec2<float>(0.0f, 0.0f), i, j);
			}

			/*
			// Outside -
			if (sphere_func > col_iso) // func-radis > 0.0f+col_iso Cell is Outside.  
			{
				// Outside Operations.
			}
			*/

			/* ADD MOUSE VELOCITY - Also On Surface & Interior Cells - 
			Eventually add this as a Bool Param Option for Sphere_Bounds Vector, and expose Vel Multipler Param.*/
			if (sphere_func <= col_iso) // func-radius <= col_iso is Ethier Inside or on Surface. 
			{
				float mouse_velMult = 1.0f; 

				// Get Stored Current Mouse Vel + to Current Grid Cell Velocity. 
				vec2<float> cur_vel = grid->getdata(i, j);
				grid->setdata(cur_vel + (mouse_vel *= mouse_velMult), i, j);
			}
		}
	}
}


/* ====================================================
	DIFFUSION
   ==================================================== */

// ** DIFFUSION-SCALAR-FIELD-LINSOLVE IMPLEMENTATION ** \\ - 

void fluidsolver_2::diffuse(grid2_scalar *grid_0, grid2_scalar *grid_1, float diff, ushort iter)
{
	// Diffusion Value for "mat" A - 
	float a = dt * diff * powf(N_dim, 2.0f);

	// Scalar Field Diffusion - 
	for (int l = 0; l < iter; l++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				// Calc New Scalar Val via Diffusion Linear Solver. // x not b...
				float b_val = (grid_0->getdata(i, j) + a * (grid_1->getdata(i - 1, j)
				 + grid_1->getdata(i + 1, j) + grid_1->getdata(i, j - 1) + grid_1->getdata(i, j + 1)))
				 / (1 + 4 * a);

				// Set New Dens Val.
				grid_1->setdata(b_val, i, j);
			}
		}

		// Call (Re-Inforce Boundary Condtions) Boundary Calc MFs on each Relaxation Iter - 
		// Only Call BoundaryCondtions on Current Step Grids. 
		edge_bounds(grid_1); // Generic Func Call, Pass Grid_1. 

		#if dospherebound == 1
		//sphere_bounds(spherebound_radius, spherebound_coliso, spherebound_offset);
		sphere_bounds_eval(grid_1, spherebound_coliso);
		#endif
	}
}

// ** DIFFUSION-VECTOR-FIELD-LINSOLVE IMPLEMENTATION ** \\ - 

void fluidsolver_2::diffuse(grid2_vector *grid_0, grid2_vector *grid_1, float diff, ushort iter)
{
	// Diffusion Value for "mat" A - 
	float a = dt * diff * powf(N_dim, 2);

	for (int l = 0; l < iter; l++)
	{
		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{

				// Per Scalar (float) Vector component diffusion - 
				float u_val = (grid_0->getdata_x(i, j) + a * (grid_1->getdata_x(i - 1, j)
					+ grid_1->getdata_x(i + 1, j) + grid_1->getdata_x(i, j - 1) + grid_1->getdata_x(i, j + 1)))
					/ (1 + 4 * a);

				float v_val = (grid_0->getdata_y(i, j) + a * (grid_1->getdata_y(i - 1, j)
					+ grid_1->getdata_y(i + 1, j) + grid_1->getdata_y(i, j - 1) + grid_1->getdata_y(i, j + 1)))
					/ (1 + 4 * a);

				vec2<float> new_vel(u_val, v_val);

				// Set New Vector Val.
				grid_1->setdata(new_vel, i, j);

			}
		}

		// Call Boundary Calc MFs on each Relaxation Iter - 
		edge_bounds(grid_1);

		#if dospherebound == 1
		sphere_bounds_eval(grid_1, spherebound_coliso);
		#endif
	}
}

// ** DIFFUSION-SCALAR-FIELD-FDM IMPLEMENTATION ** \\ - 

void fluidsolver_2::diffuse_FDM(grid2_scalar *grid_0, grid2_scalar *grid_1, float diff)
{
	float a = dt * diff * powf(N_dim, 2.0f);
	a *= 0.001f; // Mult a down allows stablilty of diffusion, without fast propogation of density. But Hacky ofc. 
	
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			float dd = grid_0->getdata(i, j) + a*(grid_0->getdata(i - 1, j) + grid_0->getdata(i + 1, j)
				+ grid_0->getdata(i, j - 1) + grid_0->getdata(i, j + 1) - 4 * grid_0->getdata(i, j));

			grid_1->setdata(dd, i, j);
		}
	}

	// Call Boundary Calc MFs on each Relaxation Iter - 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}


// ** DIFFUSION-VECTOR-FIELD-FDM IMPLEMENTATION ** \\ - 
// Wont be used so leaving out for now. 


/* ====================================================
	ADVECTION
   ==================================================== */

// Semi Lagragain Advection - Overloads depend on Grid Type Passed. (Scalar,Vector Currently).
// Linear Backtrace along Velocity, Interoplate Neighbours at Backtraced Postion, to get our new advected value. 
// Single Step. MidPoint Coming soon. 

// ** ADVECT_SL(Semi-Lagragagin)_Scalar Overload ** \\ - 
// Assume ALL Advection is done using Main Vel Field. (Not using Input VelGrid Params). 
// Advection Done along Scaled Index Space - 

void fluidsolver_2::advect_sl(grid2_scalar *grid_0, grid2_scalar *grid_1)
{

	float dt0 = dt*N_dim; // Step Distance (Velocity Scaling) of DeltaTime * Grid Length, Hence GridSpace Dt Coeff (dt0). 

	// Density (Scalar Field Advection) - 

	#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Bind Vel Components Of CurCell to Vars to save multiple acess. 
			float u = f2obj->vel->getdata_x(i, j);
			float v = f2obj->vel->getdata_y(i, j);

			// X i 
			float x = i - dt0*u;
			if (x < 0.5) x = 0.5;
			if (x > N_dim + 0.5) x = N_dim + 0.5;

			// Interp Indices i 
			int i0 = int(x); int i1 = i0 + 1;

			// Y j 
			float y = j - dt0*v;
			if (y < 0.5) y = 0.5;
			if (y > N_dim + 0.5) y = N_dim + 0.5;

			// Interp Indices j 
			int j0 = int(y); int j1 = j0 + 1;

			float s1 = x - i0; float s0 = 1 - s1;
			float t1 = y - j0; float t0 = 1 - t1;

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

			// Set New Cur Density to Current Frame Density Grid cell value - 
			grid_1->setdata(new_dens, i, j);

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

void fluidsolver_2::advect_sl(grid2_vector *grid_0, grid2_vector *grid_1)
{
	float dt0 = dt*N_dim; // Step Distance (Velocity Scaling) of DeltaTime * Grid Length, Hence GridSpace Dt Coeff (dt0). 

	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Bind Vel Components Of CurCell to Vars to save multiple acess. 
			float u = f2obj->vel->getdata_x(i, j);
			float u0 = f2obj->prev_vel->getdata_x(i, j);
			float v = f2obj->vel->getdata_y(i, j);
			float v0 = f2obj->prev_vel->getdata_y(i, j);

			// Calc Advection Vars - (Not Grid Specfic) - 
			float x = i - dt0*u;

			if (x < 0.5f) x = 0.5f;
			if (x > N_dim + 0.5f) x = N_dim + 0.5f;
			int i0 = int(x); int i1 = i0 + 1;

			float y = j - dt0*v;

			if (y < 0.5f) y = 0.5f;
			if (y > N_dim + 0.5f) y = N_dim + 0.5f;
			int j0 = int(y); int j1 = j0 + 1;

			float s1 = x - i0; float s0 = 1 - s1;
			float t1 = y - j0; float t0 = 1 - t1;

			// Interoplate Velocity U,V Components.  
			float new_u, new_v; 
			if (Parms.p_InteroplationType == Parms.Interoplation_Linear)
			{
 
				/*new_u = s0 * (t0*grid_0->getdata_x(i0, j0) + t1 * grid_0->getdata_x(i0, j1))
					+ s1 * (t0 * grid_0->getdata_x(i1, j0) + t1 * grid_0->getdata_x(i1, j1));
				new_v = s0*(t0*grid_0->getdata_y(i0, j0) + t1 * grid_0->getdata_y(i0, j1))
					+ s1 * (t0 * grid_0->getdata_y(i1, j0) + t1 * grid_0->getdata_y(i1, j1));*/

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

			// Set New Cur vel vec2 to Velocity Grid cell value - 
			vec2<float> new_vec(new_u, new_v);
			grid_1->setdata(new_vec, i, j);

		}
	}

		// Call Boundary Condtions Post Advection (Density)- 
		edge_bounds(grid_1);

		#if dospherebound == 1
		sphere_bounds_eval(grid_1, spherebound_coliso);
		#endif
}

// MID POINT ADVECTION (RK2) WIP - 
// BackTrace to MidPoint, Sample (interoplate) Velocity at MidPoint, Then BackTrace from Cell using MidPoint Vel,
// to sample (interoplate) Final BackTraced Advection Quanitiy

// Advection in Grid Index Space. Dt Scaled to N_Dim size.
// Velocity Components split for Advection Sampling and Interoplation.

// ** ADVECT_SL_RK2(Semi-Lagragagin_MidPoint)_Scalar Overload ** \\ - 
void fluidsolver_2::advect_sl_mp(grid2_scalar *grid_0, grid2_scalar *grid_1)
{
	float dt0 = dt * N_dim; // Scale DeltaTime to Grid Dimension Size. 

	// Density (Scalar Field Advection) - 

	#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Vel At Cur Cell Postion. 
			float u_P = f2obj->vel->getdata_x(i, j);
			float v_P = f2obj->vel->getdata_y(i, j);

			// BackTrace Along Negative CurCell Vel - XG - dt0 * u(CurCell)
			// XG -> Midpoint = XG - dt0 * u(XG)
			float xxm = i - (0.5 * dt0) * u_P; float yym = j - (0.5 * dt0) * v_P; // BackTrace U,V to Midpoint.
			if (xxm < 0.5) xxm = 0.5; if (xxm > N_dim + 0.5) xxm = N_dim + 0.5; // Clamp
			if (yym < 0.5) yym = 0.5; if (yym > N_dim + 0.5) yym = N_dim + 0.5;
			
			// MidPoint - Mid Indices 
			int i_mid = int(xxm); int i_mid_1 = i_mid + 1; int j_mid = int(yym); int j_mid_1 = j_mid + 1;
			// MidPoint - Mid Interp Coefficents - 
			float sm1 = xxm - i_mid; float sm = 1 - sm1;
			float tm1 = yym - j_mid; float tm = 1 - tm1;

			/* Get Mid Point Velocity (Average Neighbours).
			float u_mid = (f2obj->vel->getdata_x(i_mid, j_mid) + f2obj->vel->getdata_x(i_mid_1, j_mid_1)) / 2.0f;
			float v_mid = (f2obj->vel->getdata_y(i_mid, j_mid) + f2obj->vel->getdata_y(i_mid_1, j_mid_1)) / 2.0f;
			*/

			// Get Mid Point Velocity (Bilinear) - 
			float u_mid = sm * (tm*f2obj->vel->getdata_x(i_mid, j_mid) + tm1 * f2obj->vel->getdata_x(i_mid, j_mid_1))
				+ sm1 * (tm * f2obj->vel->getdata_x(i_mid_1, j_mid) + tm1 * f2obj->vel->getdata_x(i_mid_1, j_mid_1));

			float v_mid = sm * (tm*f2obj->vel->getdata_y(i_mid, j_mid) + tm1 * f2obj->vel->getdata_y(i_mid, j_mid_1))
				+ sm1 * (tm * f2obj->vel->getdata_y(i_mid_1, j_mid) + tm1 * f2obj->vel->getdata_y(i_mid_1, j_mid_1));

			// BackTrace Along Negative Midpoint Vel - XG - dt0 * u(midpoint)
			float xxp = i - dt0 * u_mid; float yyp = j - dt0 * v_mid;
			if (xxp < 0.5) xxp = 0.5; if (xxp > N_dim + 0.5) xxp = N_dim + 0.5;
			if (yyp < 0.5) yyp = 0.5; if (yyp > N_dim + 0.5) yyp = N_dim + 0.5;

			// MidPoint - BackTrace Indices Test - 
			int i0 = int(xxp); int i1 = i0 + 1;
			int j0 = int(yyp); int j1 = j0 + 1;
			//
			// MidPoint - BackTrace Coefficents. 
			float s1 = xxp - i0; float s0 = 1 - s1;
			float t1 = yyp - j0; float t0 = 1 - t1;

			// Bilinearlly Sample Density at backtraced postion (via MidPoint vel). 
			float new_dens = s0 * (t0*grid_0->getdata(i0, j0) + t1 * grid_0->getdata(i0, j1))
				+ s1 * (t0 * grid_0->getdata(i1, j0) + t1 * grid_0->getdata(i1, j1));

			// Set New CurFrame Grid Density to Current Frame Density Grid cell value - 
			grid_1->setdata(new_dens, i, j);

		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif

}

// ** ADVECT_SL_RK2(Semi-Lagragagin_MidPoint)_Vector Overload ** \\ - 
void fluidsolver_2::advect_sl_mp(grid2_vector *grid_0, grid2_vector *grid_1)
{
	float dt0 = dt * N_dim; // Scale DeltaTime to Grid Dimension Size. 

	// Density (Scalar Field Advection) - 

	#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Vel At Cur Cell Postion. 
			float u_P = f2obj->vel->getdata_x(i, j);
			float v_P = f2obj->vel->getdata_y(i, j);

			// BackTrace Along Negative CurCell Vel - XG - dt0 * u(CurCell)
			// XG -> Midpoint = XG - dt0 * u(XG)
			float xxm = i - (0.5 * dt0) * u_P; float yym = j - (0.5 * dt0) * v_P; // BackTrace U,V to Midpoint.
			if (xxm < 0.5) xxm = 0.5; if (xxm > N_dim + 0.5) xxm = N_dim + 0.5; // Clamp
			if (yym < 0.5) yym = 0.5; if (yym > N_dim + 0.5) yym = N_dim + 0.5;

			// MidPoint - Mid Indices 
			int i_mid = int(xxm); int i_mid_1 = i_mid + 1; int j_mid = int(yym); int j_mid_1 = j_mid + 1;
			// MidPoint - Mid Interp Coefficents - 
			float sm1 = xxm - i_mid; float sm = 1 - sm1;
			float tm1 = yym - j_mid; float tm = 1 - tm1;

			/* Get Mid Point Velocity (Average Neighbours).
			float u_mid = (f2obj->vel->getdata_x(i_mid, j_mid) + f2obj->vel->getdata_x(i_mid_1, j_mid_1)) / 2.0f;
			float v_mid = (f2obj->vel->getdata_y(i_mid, j_mid) + f2obj->vel->getdata_y(i_mid_1, j_mid_1)) / 2.0f;
			*/

			// Get Mid Point Velocity (Bilinear) - 
			float u_mid = sm * (tm*f2obj->vel->getdata_x(i_mid, j_mid) + tm1 * f2obj->vel->getdata_x(i_mid, j_mid_1))
				+ sm1 * (tm * f2obj->vel->getdata_x(i_mid_1, j_mid) + tm1 * f2obj->vel->getdata_x(i_mid_1, j_mid_1));
			float v_mid = sm * (tm*f2obj->vel->getdata_y(i_mid, j_mid) + tm1 * f2obj->vel->getdata_y(i_mid, j_mid_1))
				+ sm1 * (tm * f2obj->vel->getdata_y(i_mid_1, j_mid) + tm1 * f2obj->vel->getdata_y(i_mid_1, j_mid_1));

			// BackTrace Along Negative Midpoint Vel - XG - dt0 * u(midpoint)
			float xxp = i - dt0 * u_mid; float yyp = j - dt0 * v_mid;
			if (xxp < 0.5) xxp = 0.5; if (xxp > N_dim + 0.5) xxp = N_dim + 0.5;
			if (yyp < 0.5) yyp = 0.5; if (yyp > N_dim + 0.5) yyp = N_dim + 0.5;

			// MidPoint - BackTrace Indices Test - 
			int i0 = int(xxp); int i1 = i0 + 1;
			int j0 = int(yyp); int j1 = j0 + 1;
			//
			// MidPoint - BackTrace Coefficents. 
			float s1 = xxp - i0; float s0 = 1 - s1;
			float t1 = yyp - j0; float t0 = 1 - t1;

			// Bilinearly Sample Velocity Components at Backtraced postion (via MidPoint Vel). 
			float new_u = s0 * (t0*grid_0->getdata_x(i0, j0) + t1 * grid_0->getdata_x(i0, j1))
				+ s1 * (t0 * grid_0->getdata_x(i1, j0) + t1 * grid_0->getdata_x(i1, j1));

			float new_v = s0 * (t0*grid_0->getdata_y(i0, j0) + t1 * grid_0->getdata_y(i0, j1))
				+ s1 * (t0 * grid_0->getdata_y(i1, j0) + t1 * grid_0->getdata_y(i1, j1));

			vec2<float> new_vec(new_u, new_v);

			// Set New Cur vel vec2 to Velocity Grid cell value - 
			grid_1->setdata(new_vec, i, j);

		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval (grid_1, spherebound_coliso);
	#endif

}

// ADVECTION - BackTracing in Grid Space - 

// ** ADVECT_SL_RK2(Semi-Lagragagin_MidPoint)_GRIDSPACE(GS) Scalar Overload ** \\ -
// Grid Space (oppose to Index Space) Backtrcing. 
void fluidsolver_2::advect_sl_mp_GS(grid2_scalar *grid_0, grid2_scalar *grid_1)
{
	// Index-Grid-Index Space Lambdas - 
	auto idx_indexToGrid = [&](int i, int j) -> vec2<float>
	{
		return vec2<float>((float)i / (float)N_dim, (float)j / (float)N_dim);
	};

	auto idx_gridToIndex = [&](auto x, auto y) -> vec2<int>
	{
		return vec2<int>((int) x * N_dim, (int) y * N_dim);
	};

	// Density (Scalar Field Advection) - 

	#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Cur Cell Index to Grid Space - 
			vec2<float> gsIDX = idx_indexToGrid(i, j);

			// Vel At Cur Cell Index. 
			float u_P = f2obj->vel->getdata_x(i, j);
			float v_P = f2obj->vel->getdata_y(i, j);

			// Back Trace in Grid Space. 
			// BackTrace Along Negative CurCell Vel - XG - dt0 * u(CurCell)
			// XG -> Midpoint = XG - dt0 * u(XG)
			float xxm_gs = gsIDX.x - (0.5 * dt) * u_P; float yym_gs = gsIDX.y - (0.5 * dt) * v_P; // BackTrace U,V to Midpoint Postion.

			// Clamp To 0-1 Grid Space to avoid out of bounds - 
			xxm_gs = solver_utils::clamp(xxm_gs, 0.0f, 1.0f); yym_gs = solver_utils::clamp(yym_gs, 0.0f, 1.0);

			// MidPoint Postion to Index Space - 
			vec2<int> mpIDX = idx_gridToIndex(xxm_gs, yym_gs);  // Ideally vec2 would return ints via Templated Types avoid float-int casting. 

			// Get Midpoint i,j indices
			int mpidx_i = (int) mpIDX.x; int mpidx_i1 = mpidx_i + 1;
			int mpidx_j = (int) mpIDX.y; int mpidx_j1 = mpidx_j + 1;

			// MidPoint - Mid Interp Coefficents - (Scalar (Grid->Index) - Integer (Grid->Index)). 
			float sm1 = mpIDX.x - mpidx_i; float sm = 1 - sm1;
			float tm1 = mpIDX.y - mpidx_j; float tm = 1 - tm1;

			// Sample Mid Point Velocity Components (Bilinear) - 
			float u_mid = sm * (tm*f2obj->vel->getdata_x(mpidx_i, mpidx_j) + tm1 * f2obj->vel->getdata_x(mpidx_i, mpidx_j1))
				+ sm1 * (tm * f2obj->vel->getdata_x(mpidx_i1, mpidx_j) + tm1 * f2obj->vel->getdata_x(mpidx_i1, mpidx_j1));

			float v_mid = sm * (tm*f2obj->vel->getdata_y(mpidx_i, mpidx_j) + tm1 * f2obj->vel->getdata_y(mpidx_i, mpidx_j1))
				+ sm1 * (tm * f2obj->vel->getdata_y(mpidx_i1, mpidx_j) + tm1 * f2obj->vel->getdata_y(mpidx_i1, mpidx_j1));

			// BackTrace Along Negative Midpoint Vel - XG - dt0 * u(midpoint)
			float xxp_gs = gsIDX.x - dt * u_mid; float yyp_gs = gsIDX.y - dt * v_mid;

			// Clamp To 0-1 Grid Space to avoid out of bounds - 
			xxp_gs = solver_utils::clamp(xxp_gs, 0.0f, 1.0f); yyp_gs = solver_utils::clamp(yyp_gs, 0.0f, 1.0);

			// BackTraced Postion to Index Space - 
			vec2<int> btIDX = idx_gridToIndex(xxp_gs, yyp_gs);
			// Get BackTraced i,j indices -
			int btidx_i = (int) btIDX.x; int btidx_i1 = btidx_i + 1;
			int btidx_j = (int) btIDX.y; int btidx_j1 = btidx_j + 1;

			// BackTrace Coefficents - (Scalar (Grid->Index) - Integer (Grid->Index)). 
			float s1 = btIDX.x - btidx_i; float s0 = 1 - s1;
			float t1 = btIDX.y - btidx_j; float t0 = 1 - t1;

			// Bilinearlly Sample Density at backtraced postion (via MidPoint vel). 
			float new_dens = s0 * (t0*grid_0->getdata(btidx_i, btidx_j) + t1 * grid_0->getdata(btidx_i, btidx_j1))
				+ s1 * (t0 * grid_0->getdata(btidx_i1, btidx_j) + t1 * grid_0->getdata(btidx_i1, btidx_j1));

			// Set New CurFrame Grid Density to Current Frame Density Grid cell value - 
			grid_1->setdata(new_dens, i, j);

		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}

// ** ADVECT_SL_RK2(Semi-Lagragagin_MidPoint)_GRIDSPACE(GS) Vector Overload ** \\ -
// Grid Space (oppose to Index Space) Backtrcing. 

void fluidsolver_2::advect_sl_mp_GS(grid2_vector *grid_0, grid2_vector *grid_1)
{
	// Index-Grid-Index Space Lambdas - 
	auto idx_indexToGrid = [&](int i, int j) -> vec2<float>
	{
		return vec2<float>((float)i / (float)N_dim, (float)j / (float)N_dim);
	};

	auto idx_gridToIndex = [&](auto x, auto y) -> vec2<int>
	{
		return vec2<int>(x * N_dim, y * N_dim);
	};

	// Density (Scalar Field Advection) - 

	#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Cur Cell Index to Grid Space - 
			vec2<float> gsIDX = idx_indexToGrid(i, j);

			// Vel At Cur Cell Index. 
			float u_P = f2obj->vel->getdata_x(i, j);
			float v_P = f2obj->vel->getdata_y(i, j);

			// Back Trace in Grid Space. 
			// BackTrace Along Negative CurCell Vel - XG - dt0 * u(CurCell)
			// XG -> Midpoint = XG - dt0 * u(XG)
			float xxm_gs = gsIDX.x - (0.5 * dt) * u_P; float yym_gs = gsIDX.y - (0.5 * dt) * v_P; // BackTrace U,V to Midpoint Postion.

			// Clamp To 0-1 Grid Space to avoid out of bounds - 
			xxm_gs = solver_utils::clamp(xxm_gs, 0.0f, 1.0f); yym_gs = solver_utils::clamp(yym_gs, 0.0f, 1.0);

			// MidPoint Postion to Index Space - 
			vec2<int> mpIDX = idx_gridToIndex(xxm_gs, yym_gs);  // Ideally vec2 would return ints via Templated Types avoid float-int casting. 
			// Get Midpoint i,j indices
			int mpidx_i = (int)mpIDX.x; int mpidx_i1 = mpidx_i + 1;
			int mpidx_j = (int)mpIDX.y; int mpidx_j1 = mpidx_j + 1;

			// MidPoint Interp Coefficents - (Scalar (Grid->Index) - Integer (Grid->Index)).  
			float sm1 = mpIDX.x - mpidx_i; float sm = 1 - sm1;
			float tm1 = mpIDX.y - mpidx_j; float tm = 1 - tm1;

			// Sample Mid Point Velocity Components (Bilinear) - 
			float u_mid = sm * (tm*f2obj->vel->getdata_x(mpidx_i, mpidx_j) + tm1 * f2obj->vel->getdata_x(mpidx_i, mpidx_j1))
				+ sm1 * (tm * f2obj->vel->getdata_x(mpidx_i1, mpidx_j) + tm1 * f2obj->vel->getdata_x(mpidx_i1, mpidx_j1));

			float v_mid = sm * (tm*f2obj->vel->getdata_y(mpidx_i, mpidx_j) + tm1 * f2obj->vel->getdata_y(mpidx_i, mpidx_j1))
				+ sm1 * (tm * f2obj->vel->getdata_y(mpidx_i1, mpidx_j) + tm1 * f2obj->vel->getdata_y(mpidx_i1, mpidx_j1));

			u_mid = f2obj->vel->getdata_x(mpidx_i, mpidx_j);
			v_mid = f2obj->vel->getdata_y(mpidx_i, mpidx_j);

			// BackTrace Along Negative Midpoint Vel - XG - dt0 * u(midpoint)
			float xxp_gs = gsIDX.x - dt * u_mid; float yyp_gs = gsIDX.y - dt * v_mid;

			// Clamp To 0-1 Grid Space to avoid out of bounds - 
			xxp_gs = solver_utils::clamp(xxp_gs, 0.0f, 1.0f); yyp_gs = solver_utils::clamp(yyp_gs, 0.0f, 1.0);

			// BackTraced Postion to Index Space - 
			vec2<int> btIDX = idx_gridToIndex(xxp_gs, yyp_gs);
			// Get BackTraced i,j indices -
			int btidx_i = (int)btIDX.x; int btidx_i1 = btidx_i + 1;
			int btidx_j = (int)btIDX.y; int btidx_j1 = btidx_j + 1;

			// BackTrace Coefficents - (Scalar (Grid->Index) - Integer (Grid->Index)). 
			float s1 = gsIDX.x - btidx_i; float s0 = 1 - s1;
			float t1 = gsIDX.y - btidx_j; float t0 = 1 - t1;

			// Bilinearly Sample Velocity Components at Backtraced postion (via MidPoint Vel). 
			float new_u = s0 * (t0*grid_0->getdata_x(btidx_i, btidx_j) + t1 * grid_0->getdata_x(btidx_i, btidx_j1))
				+ s1 * (t0 * grid_0->getdata_x(btidx_i1, btidx_j) + t1 * grid_0->getdata_x(btidx_i1, btidx_j1));

			float new_v = s0 * (t0*grid_0->getdata_y(btidx_i, btidx_j) + t1 * grid_0->getdata_y(btidx_i, btidx_j1))
				+ s1 * (t0 * grid_0->getdata_y(btidx_i1, btidx_j) + t1 * grid_0->getdata_y(btidx_i1, btidx_j1));

			new_u = grid_0->getdata_x(btidx_i, btidx_j);
			new_v = grid_0->getdata_y(btidx_i, btidx_j);

			vec2<float> new_vec(new_u, new_v);

			// Set New Cur vel vec2 to Velocity Grid cell value - 
			grid_1->setdata(new_vec, i, j);


		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	edge_bounds(grid_1);

	#if dospherebound == 1
	sphere_bounds_eval(grid_1, spherebound_coliso);
	#endif
}


/*
//  ADVECT_SL_RK2(Semi-Lagragagin_MidPoint)_Scalar Overload \\ - 
// WITH BICUBIC INTEROPLATION - WIP
void fluidsolver_2::advect_sl_mp_bc(grid2_scalar *grid_0, grid2_scalar *grid_1)
{
	float dt0 = dt * N_dim; // Scale DeltaTime to Grid Dimension Size. 

	// Density (Scalar Field Advection) - 
#pragma omp parallel for num_threads(omp_get_max_threads()) 
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Vel At Cur Cell Postion. 
			float u_P = f2obj->vel->getdata_x(i, j);
			float v_P = f2obj->vel->getdata_y(i, j);

			// BackTrace Along Negative CurCell Vel - XG - dt0 * u(CurCell)
			// XG -> Midpoint = XG - dt0 * u(XG)
			float xxm = i - (0.5 * dt0) * u_P; float yym = j - (0.5 * dt0) * v_P; // BackTrace U,V to Midpoint.
			if (xxm < 0.5) xxm = 0.5; if (xxm > N_dim + 0.5) xxm = N_dim + 0.5; // Clamp
			if (yym < 0.5) yym = 0.5; if (yym > N_dim + 0.5) yym = N_dim + 0.5;

			// MidPoint - Mid Indices 
			int i_mid = int(xxm); int i_mid_1 = i_mid + 1; int j_mid = int(yym); int j_mid_1 = j_mid + 1;
			// MidPoint - Mid Interp Coefficents - 
			float sm1 = xxm - i_mid; float sm = 1 - sm1;
			float tm1 = yym - j_mid; float tm = 1 - tm1;

			/* Get Mid Point Velocity (Average Neighbours).
			float u_mid = (f2obj->vel->getdata_x(i_mid, j_mid) + f2obj->vel->getdata_x(i_mid_1, j_mid_1)) / 2.0f;
			float v_mid = (f2obj->vel->getdata_y(i_mid, j_mid) + f2obj->vel->getdata_y(i_mid_1, j_mid_1)) / 2.0f;
			*/

/*
			// Get Mid Point Velocity (Bilinear) - 
			//float u_mid = sm * (tm*f2obj->vel->getdata_x(i_mid, j_mid) + tm1 * f2obj->vel->getdata_x(i_mid, j_mid_1))
				//+ sm1 * (tm * f2obj->vel->getdata_x(i_mid_1, j_mid) + tm1 * f2obj->vel->getdata_x(i_mid_1, j_mid_1));
			//float v_mid = sm * (tm*f2obj->vel->getdata_y(i_mid, j_mid) + tm1 * f2obj->vel->getdata_y(i_mid, j_mid_1))
			//	+ sm1 * (tm * f2obj->vel->getdata_y(i_mid_1, j_mid) + tm1 * f2obj->vel->getdata_y(i_mid_1, j_mid_1));

			float u_mid = solver_utils::bicubic_V(f2obj->vel, sm, sm1, tm, tm1, i, j, 0);
			float v_mid = solver_utils::bicubic_V(f2obj->vel, sm, sm1, tm, tm1, i, j, 1);

			// BackTrace Along Negative Midpoint Vel - XG - dt0 * u(midpoint)
			float xxp = i - dt0 * u_mid; float yyp = j - dt0 * v_mid;
			if (xxp < 0.5) xxp = 0.5; if (xxp > N_dim + 0.5) xxp = N_dim + 0.5;
			if (yyp < 0.5) yyp = 0.5; if (yyp > N_dim + 0.5) yyp = N_dim + 0.5;

			// MidPoint - BackTrace Indices Test - 
			int i0 = int(xxp); int i1 = i0 + 1;
			int j0 = int(yyp); int j1 = j0 + 1;
			//
			// MidPoint - BackTrace Coefficents. 
			float s1 = xxp - i0; float s0 = 1 - s1;
			float t1 = yyp - j0; float t0 = 1 - t1;

			// Bilinearlly Sample Density at backtraced postion (via MidPoint vel). 
			//float new_dens = s0 * (t0*grid_0->getdata(i0, j0) + t1 * grid_0->getdata(i0, j1))
				//+ s1 * (t0 * grid_0->getdata(i1, j0) + t1 * grid_0->getdata(i1, j1));

			float new_dens = solver_utils::bicubic_S(grid_0, sm, sm1, tm, tm1, i, j);

			// Set New CurFrame Grid Density to Current Frame Density Grid cell value - 
			grid_1->setdata(new_dens, i, j);

		}
	}

	// Call Boundary Condtions Post Advection (Scalar)- 
	edge_bounds(grid_1);

#if dospherebound == 1
	sphere_bounds(grid_1, spherebound_radius, spherebound_coliso, spherebound_offset);
#endif

}
*/

/* ====================================================
				FORCES AND MISC - 
	==================================================== */

// Loop Through Vel Grid, and Add passed force to velocity * dt. 
void fluidsolver_2::vel_force(vec2<float> ff, float dtt)
{
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			vec2<float> temp (ff.x, ff.y); // UN NEGATE AXIS for coord fix/dbg ?
			f2obj->integrate_force(temp, dtt, i, j);   // Call add_force MF in fluid_obj.
		}
	}
}

// Pass Custom Force Callback (Lambda) to run over each vector grid cell and integrate as force. 
void fluidsolver_2::custom_force(grid2_vector *grid, std::function <vec2<float>(vec2<int> idx)> &force)
{
	#pragma omp parallel for 
	for (int j = 1; j <= N_dim; j++)
	{
		#pragma omp parallel for 
		for (int i = 1; i <= N_dim; i++)
		{
			// Call Callback - Get Force Vector. 
			vec2<float> force_vec = force(vec2<int>(i, j));
			// Integrate Force -
			f2obj->integrate_force(force_vec, dt, i, j);
		}
	}
}

/*	====================================================
	VORTICITY CONFINEMENT
==================================================== */

// Compute and Integrate Vorticity Confinement Force, from/to Velocity Field.

// NOTES/Learnt - 
// CURL in 2D is a Scalar, not a vector. As theres no Cross Product for single XY Plane Vel.

/*
void fluidsolver_2::vorticty_confine(float strength)
{
	// VCForce = strength * cell size * (N cross W) ... Intergrate to Vel Field over Dt. 
	float dx = 1.0f / N_dim; // (h / cell size)

	// Allocate Vorticity Grids (Per Call)
	vort = new grid2_vector(X_dim, Y_dim, edge_size, 6, 1.0f); // Vorticity Field Grid
	deltavort = new grid2_vector(X_dim, Y_dim, edge_size, 7, 1.0f); // Vorticity Gradient Field Grid. 

	// Cap f2obj ptr by ref [&]. 
	std::function<float(int, int)> curl2d = [&](int i , int j) -> float 
	{
		return float(f2obj->vel->getdata_x(i, j + 1) - f2obj->vel->getdata_x(i, j - 1)) 
		+ (f2obj->vel->getdata_y(i + 1, j) - f2obj->vel->getdata_x(i, j - 1));
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
		curl_length = (f2obj->vel->getdata_x(i, j + 1) - f2obj->vel->getdata_x(i, j - 1)) - (f2obj->vel->getdata_y(i + 1, j) - f2obj->vel->getdata_x(i, j - 1));

	*/

	/*
	// Caluclate Vorticity (Curl w) From Velocity Grid, Write to Vorticity Grid. // Exclude Ghost/Edge Cells. 
	for (int j = 1; j < N_dim; j++)
	{
		for (int i = 1; i < N_dim; i++)
		{
		    // Central Diffrence To Calc Curl from XY Velocity - / i +/- (j) constant, then viceversa. Divide by h (dx) *2 Because Over 2 Cells +/- Distance h. 
		//	vec2 w_xy = (f2obj->vel->getdata(i + 1, j) - f2obj->vel->getdata(i - 1, j) /= dx*2) - (f2obj->vel->getdata(i, j + 1) - f2obj->vel->getdata(i, j - 1) /= dx*2);
			vec2 w_xy = (f2obj->vel->getdata(i + 1, j) - f2obj->vel->getdata(i - 1, j)) - (f2obj->vel->getdata(i, j + 1) - f2obj->vel->getdata(i, j - 1));
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
			vec2 grad(xx, yy); 
			//grad.normalize(); // Normalize Gradient N

			// Resulting Vorticity Confinement Vector W X N
			vec2 vc_vec = vort->getdata(i, j).cross(grad); 
			vec2 vc_force = vc_vec *= (strength * dx);

			// VCForce = Strength * dx * (W X N)
			// Vel += VCForce * dt; // Integrate VC Force to vel. 

			// Debuggin currently so messing with what I Intergrate to vel field... 

			// Get Current Cell Velocity From VelField - 
			vec2 vel = f2obj->vel->getdata(i, j);
			// vel = vel + (vc_force *= dt) ; // Integrate vc force. 
			vec2 vcurl = vort->getdata(i, j);
			//vcurl.normalize();
			vel = vel + (vcurl *= 0.002f); // Integrate vc force. 
			f2obj->vel->setdata(vel, i, j); // Set to new Vel Value. 
		//	f2obj->vel->setdata(vel + vort->getdata(i,j), i, j); // Set to new Vel Value. 
			

			//float xx = std::abs(curl2d(i, j+1)) - std::abs(curl2d(i, j-1));
			float xx = curl2d(i, j - 1) - curl2d(i, j + 1);
			float yy = curl2d(i+1, j) - curl2d(i-1, j);
			vec2 dir(xx, yy);

		//	float f = (10.0f / dir.length()) + 1e-05f;
		//	dir *= f; 
			vec2 old_vel = f2obj->vel->getdata(i, j);
		//	float dtc = dt * curl2d(i, j);
		//	vec2 vc = dir *= dtc; 
		//	vec2 new_vel = old_vel + vc; 
			f2obj->vel->setdata(old_vel + ((dir += 1e-05f) *= 5.0f * dt), i, j);


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

// Vorticity Confinement On The Fly (No Writing to Vorticity/Curl Grids, Curl, Curl Grad and Resulting VC Force
// calculated all within one loop, via Curl Lambda to calc  Curl (w) and Curl Gradient (w'). 
void fluidsolver_2::vorticty_confine_otf(float strength)
{
	// Calculate 2D Velocity Curl Scalar - Without Central Diff Offsets, for Calling within Curl gradient CentralDiff calc lines. 
	std::function<float(int, int)> curl2d = [&](int i, int j) -> float // Cap f2obj ptr by ref [&]. 
	{
		// i j offsets will be set on call. 
		return float(f2obj->vel->getdata_x(i, j) - f2obj->vel->getdata_x(i, j))
			+ (f2obj->vel->getdata_y(i, j) - f2obj->vel->getdata_x(i, j));
		// Should be / dx*2 ie central diff over h. 
	};

	// Curl With ij central diff offsets, for calling when eval curl only. 
	std::function<float(int, int)> curl2d_wo = [&](int i, int j) -> float // Cap f2obj ptr by ref [&]. 
	{
		// i j set with offset.
		return float(f2obj->vel->getdata_x(i, j+1) - f2obj->vel->getdata_x(i, j-1))
			+ (f2obj->vel->getdata_y(i+1, j) - f2obj->vel->getdata_x(i-1, j));
		// Should be / dx*2 ie central diff over h. 
	};

	for (int j = 1; j < N_dim - 1; j++)
	{
		for (int i = 1; i < N_dim - 1; i++)
		{
			// w' from w Calc both curl and gradient curl per axis using central diff. 
			float xx = curl2d(i, j-1) - curl2d(i, j+1);
			float yy = curl2d(i+1, j) - curl2d(i-1, j);
			vec2<float> dir(xx, yy); // resulting w' grad vector. 
			f2obj->vc->setdata(dir, i, j); // Pass To VC Grid for viz. Pre Normalization (will be clampped vals on GPU).
			dir += 1e-05; // Prevent /0 errors on normalization.
			dir.normalize();

			// Integrate VC to new vel value.
			vec2<float> old_vel = f2obj->vel->getdata(i, j);
			dir *= strength; // Mult By Strength Scalar Coeff. 
			//f2obj->vel->setdata(old_vel + (dir *= (dt * curl2d(i, j))), i, j);
			// v + strength*(Curl (w) and GradCurl (w')*dt) = v1 (v with VC). 
			f2obj->vel->setdata(old_vel + ((dir *= dt) *= curl2d_wo(i,j)), i, j);
		}
	}
}
// End of Vorticity_Confine_otf Member Func Implementation.

void fluidsolver_2::vorticity_confine_B(float strength)
{
	// Alloc Vorticity Grid - 
	vort = new grid2_scalar(x_s, y_s, e_s, 5, 1.0f);

	float dx = 1.0f / N_dim; // cellsize dx (h). 

	// Calculate 2D Velocity Curl Scalar - At Given Cell Indices- 
	std::function<float(int, int)> curl2d = [&](int i, int j) -> float // Cap f2obj ptr by ref [&]. 
	{
		return float(
			(f2obj->vel->getdata_x(i, j+1) - f2obj->vel->getdata_x(i, j-1) / (2 * dx))
			- 
			(f2obj->vel->getdata_y(i+1, j) - f2obj->vel->getdata_y(i-1, j) / (2 * dx))
			);
	};

	// Calculate 2D Vorticity Gradient Vector - At Given Cell Indices- - 
	std::function<vec2<float>(int, int)> grad2d = [&](int i, int j) -> vec2<float>
	{
		float xx = (vort->getdata(i + 1, j) - vort->getdata(i - 1, j)) / (2 * dx); // 2.0f * dx;
		float yy = (vort->getdata(i, j + 1) - vort->getdata(i, j - 1)) / (2 * dx); // 2.0f * dx;
		vec2<float> grad(xx, yy);
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
			vec2<float> N = grad2d(i, j);

			// Write W (as vec2(W,W)) OR N (VC Gradient) to FluidObj vc grid, for Viz/Rendering via RenderObject. 
		   // f2obj->vc->setdata(vec2(W, W), i, j); // W VIZ
			f2obj->vc->setdata(vec2<float>(N) *= 1.0f , i, j); // N VIZ

			vec2<float> VC = N; VC *= W; 
			//f2obj->vc->setdata(VC, i, j); // N VIZ
			//VC *= (strength * dx);
			VC *= strength; // Final Vort Con Vector

			// Intergrate - 
			vec2<float> vel_0 = f2obj->vel->getdata(i, j);
			//vec2 vel_1 = vel_0 + (VC *= dt);
			N *= strength; //10.0f; // strength; // Mult By Stength Coeff. 
			vec2<float> vel_1 = vel_0 + (N *= (curl2d(i, j) * dt));

			// Set Vel - 
			f2obj->vel->setdata(vel_1, i, j);

		}
	}

	// Delete/Dealloc Vort Grid. 
	delete vort; vort = nullptr; 
}

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

void fluidsolver_2::project(int iter)
{
	float h = 1.0f/N_dim; // CellSize 1.0f/N (N_dim); 

	// Ensure previous Divergence and Pressure Temp Grids have been deleted before new grid ptr assignement/Alloc. 
	del_divergence(); del_pressure(); 

	// Init Solvers Own Divergence and Pressure Fields, Only needed for Scope of this Function.  
	divergence = new grid2_scalar(x_s, y_s, e_s, 4, 1.0f);
	pressure = new grid2_scalar(x_s, y_s, e_s, 5, 1.0f);

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
			float div = -0.5 * h * (f2obj->vel->getdata_x(i + 1, j) - f2obj->vel->getdata_x(i - 1, j)
			+ f2obj->vel->getdata_y(i, j + 1) - f2obj->vel->getdata_y(i, j - 1));
			
			// Set Divergence Cell Value. 
			divergence->setdata(div, i, j);

			// Zero Out Pressure Grid, as Inital Value PreLinSolve. (Index Based oppose to calling grid2_scalar->clear()).
			pressure->setdata(0.0f, i, j);

			// Write Inital PreProjected VelField to Grid For dbg - 
			f2obj->preproj_vel->setdata(f2obj->vel->getdata(i, j), i, j);

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
	for (int k = 0; k < iter; k++)
	{
		double error = 0.0f; 
		double total_error = 0.0f; 

		for (int j = 1; j <= N_dim; j++)
		{
			for (int i = 1; i <= N_dim; i++)
			{
				float n0 = pressure->getdata(i, j); // n (Pressure (k(n)) for SOR)

				float pres = (divergence->getdata(i, j) + pressure->getdata(i - 1, j) + pressure->getdata(i + 1, j) +
				pressure->getdata(i, j - 1) + pressure->getdata(i, j + 1)) / 4.0f; 
				
				float n1 = pres; // n+1 (Pressure (k(n+1)) for SOR)

				// Use SOR or not - 
				if (Parms.p_ProjectionType == Parms.Project_GaussSeidel_SOR)
				{
					// SOR 
					float alpha = Parms.p_SOR_alpha;
					float sor_pres = alpha * n1 + (1 - alpha) * n0;
					pressure->setdata(sor_pres, i, j);
				}
				else if (Parms.p_ProjectionType == Parms.Project_GaussSeidel)
				{
					pressure->setdata(pres, i, j);
				}

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
	
	// SUBTRACT PRESSURE GRADEINT FROM VELOCITY FIELD \\  

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
			float new_vel_x = f2obj->vel->getdata_x(i, j) - grad_x;
			f2obj->vel->setdata_x(new_vel_x, i, j);
			float new_vel_y = f2obj->vel->getdata_y(i, j) - grad_y;
			f2obj->vel->setdata_y(new_vel_y, i, j);

		}
	}

	// Call and Enforce Boundary Condtions After Projection on Vel Field - 
	edge_bounds(f2obj->vel);
	#if dospherebound == 1
	sphere_bounds_eval(f2obj->vel, spherebound_coliso);
	#endif
		
	// Delete Solver Temp Pressure and Divergence Fields - 
	del_divergence(); del_pressure();

}
// End of Velocity Projection Implementation (GS + SOR).


/* PROJECTION - JACOBI to Solve Pressure Poission Equation -
	Allows Multithreading as Cells (inner i,j presure loop) Only lookup Previous Pressure Values, but results in slower Convergence. 
	Divergence and Pressure Gradient Subtraction loops are MT, these are thread safe.*/

void fluidsolver_2::project_jacobi(int iter)
{
	float h = 1.0f/N_dim; // CellSize 1.0f/N (N_dim); 

	// Ensure previous Divergence and Pressure Temp Grids have been deleted before new grid ptr assignement/Alloc. 
	del_divergence(); del_pressure();
	
	// Alloc New DP Grids. 
	divergence = new grid2_scalar(x_s, y_s, e_s, 4, 1.0f);
	pressure = new grid2_scalar(x_s, y_s, e_s, 5, 1.0f);
	pressure_1 = new grid2_scalar(x_s, y_s, e_s, 6, 1.0f);

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
			float div = -0.5 * h * (f2obj->vel->getdata_x(i + 1, j) - f2obj->vel->getdata_x(i - 1, j)
			+ f2obj->vel->getdata_y(i, j + 1) - f2obj->vel->getdata_y(i, j - 1));
			
			// Set Divergence Cell Value. 
			divergence->setdata(div, i, j);

			// Zero Out Pressure Grid, as Inital Value PreLinSolve. (Index Based oppose to calling grid2_scalar->clear()).
			pressure->setdata(0.0f, i, j);

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
	// Jacobi For Possion Pressure Solve. (Write to Seperate "Sratch Grid", Only Swap at end of kth Iteration.
	// pressure == main pressure grid (to read from), pressure_1 == scratch pressure grid (to write to and then swap). 

	for (int k = 0; k < iter; k++)
	{
		#pragma omp parallel for
		for (int j = 1; j <= N_dim; j++)
		{
			#pragma omp parallel for
			for (int i = 1; i <= N_dim; i++)
			{
				float pres_n0 = pressure->getdata(i, j);

				float pres = (divergence->getdata(i, j) + pressure->getdata(i - 1, j) + pressure->getdata(i + 1, j) +
				pressure->getdata(i, j - 1) + pressure->getdata(i, j + 1)) / 4.0f;

				pressure_1->setdata(pres, i, j); 
			}
		}
		// Call Boundary Condtion Functions On Pressure After Each Pressure Field Iteration.
		edge_bounds(pressure_1);
		#if dospherebound == 1
		//sphere_bounds_scalar(pressure_1, spherebound_radius, spherebound_coliso, spherebound_offset);
		sphere_bounds_eval(pressure_1,spherebound_coliso); // Optimized SphereBounds for pressure calc. 
		#endif	

		// Swap Pressure Grid with Scratch Grid at end of k Iter After internal i,j MT'd Jacobi Projection Iteration is complete.
		// (This is Single Threaded to ensure whole grid is swapped correctly together, and not within a Multithreaded Outer loop).
		pressure_1->swap(pressure);
		// Hence removal of Outer MT Kth Loop and OMP Crticial which was not correct, as Pressure grid was swapping atomically per x threads, not as a singlethread. 
	}
	
	// SUBTRACT PRESSURE GRADEINT FROM VELOCITY FIELD \\ -

	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Partial Derivatves for Each Pressure Gradient Components -
			float grad_x = 0.5 * (pressure->getdata(i + 1, j) - pressure->getdata(i - 1, j)) / h;
			float grad_y = 0.5 * (pressure->getdata(i, j + 1) - pressure->getdata(i, j - 1)) / h;

			// Subtract Gradient Components from Velocity Field Components and set to Velocity-
			float new_vel_x = f2obj->vel->getdata_x(i, j) - grad_x;
			f2obj->vel->setdata_x(new_vel_x, i, j);
			float new_vel_y = f2obj->vel->getdata_y(i, j) - grad_y;
			f2obj->vel->setdata_y(new_vel_y, i, j);
		}
	}
	// Call and Enforce Boundary Condtions After Projection on Vel Field - 
	edge_bounds(f2obj->vel);

	#if dospherebound == 1
	sphere_bounds_eval(f2obj->vel, spherebound_coliso);
	#endif

	// TEMP FIELD DELETION/DEALLOCATION - 
	del_divergence(); del_pressure();
	
}
// End of Velocity Projection (Jacobi) Implementation.

// PROJECT - Gauss-Seidel + SOR - SIMD TESTING - 

void fluidsolver_2::project_SIMD(int iter)
{
	float h = 1.0f / N_dim; // CellSize 1.0f/N (N_dim); 

	// Ensure previous Divergence and Pressure Temp Grids have been deleted before new grid ptr assignement/Alloc. 
	del_divergence(); del_pressure();

	// Init Solvers Own Divergence and Pressure Fields, Only needed for Scope of this Function.  
	divergence = new grid2_scalar(x_s, y_s, e_s, 4, 1.0f);
	pressure = new grid2_scalar(x_s, y_s, e_s, 5, 1.0f);

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
			float div = -0.5 * h * (f2obj->vel->getdata_x(i + 1, j) - f2obj->vel->getdata_x(i - 1, j)
				+ f2obj->vel->getdata_y(i, j + 1) - f2obj->vel->getdata_y(i, j - 1));

			// Set Divergence Cell Value. 
			divergence->setdata(div, i, j);

			// Zero Out Pressure Grid, as Inital Value PreLinSolve. (Index Based oppose to calling grid2_scalar->clear()).
			pressure->setdata(0.0f, i, j);

			// Write Inital PreProjected VelField to Grid For dbg - 
			f2obj->preproj_vel->setdata(f2obj->vel->getdata(i, j), i, j);

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
			float new_vel_x = f2obj->vel->getdata_x(i, j) - grad_x;
			f2obj->vel->setdata_x(new_vel_x, i, j);
			float new_vel_y = f2obj->vel->getdata_y(i, j) - grad_y;
			f2obj->vel->setdata_y(new_vel_y, i, j);

		}
	}

	// Call and Enforce Boundary Condtions After Projection on Vel Field - 
	edge_bounds(f2obj->vel);
	#if dospherebound == 1
	sphere_bounds_eval(f2obj->vel, spherebound_coliso);
	#endif

	// Delete Solver Temp Pressure and Divergence Fields - 
	del_divergence(); del_pressure();

}
// End of Velocity Projection Implementation (GS + SOR) SIMD.


/* ====================================================
	DESNITY SOLVE STEP - 
	==================================================== */
// Implementation of Simulation Solve step of Density solve operations. 

void fluidsolver_2::density_step(int diff_iter, float diffA, bool dodiff)
{
	// DIFFUSE Density - 
	// If Using Diffusion, If Not Need to Manually set Cur to Prev, rather than swapping.
	if (Parms.p_Do_Dens_Diff == true)
	{
		f2obj->dens->swap(f2obj->prev_dens); // Swap Density With Prev_Density. 
		// Gauss-Seidel Iterative Density (Scalar) Diffusion - 
		diffuse(f2obj->prev_dens, f2obj->dens, diffA, diff_iter);

		// Finite Diffrence Unstable Density (Scalar) Diffusion - 
		//diffuse_scalar_FDM(f2obj->prev_dens, f2obj->dens, diffA);
		f2obj->dens->swap(f2obj->prev_dens); // Re-Swap Density With Prev_Density. 
	}
	else if (Parms.p_Do_Dens_Diff == false)
	{
		// Use "SetCurToPrev" Funcs to Copy Grids, Oppose to via Diffusion - 
		f2obj->setcurtoprev_dens();

		if (Parms.p_useColour)
		{
			// SetCurtoPrev - Colour - 
			f2obj->setcurtoprev(f2obj->c_R_prev, f2obj->c_R);
			f2obj->setcurtoprev(f2obj->c_G_prev, f2obj->c_G);
			f2obj->setcurtoprev(f2obj->c_B_prev, f2obj->c_B);
		}

	}

	// ADVECT Density - 

	if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_Euler)
	{
		advect_sl(f2obj->prev_dens, f2obj->dens); 
	}
	else if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_RK2)
	{
		advect_sl_mp(f2obj->prev_dens, f2obj->dens);
		advect_sl_mp_GS(f2obj->prev_dens, f2obj->dens);
	}
	
	if (Parms.p_useColour)
	{
		// Advect Colour -
		if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_Euler)
		{
			advect_sl(f2obj->c_R_prev, f2obj->c_R);
			advect_sl(f2obj->c_G_prev, f2obj->c_G);
			advect_sl(f2obj->c_B_prev, f2obj->c_B);
		}
		else if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_RK2)
		{
			// RK2 (MidPoint) 
			advect_sl_mp(f2obj->c_R_prev, f2obj->c_R);
			advect_sl_mp(f2obj->c_G_prev, f2obj->c_G);
			advect_sl_mp(f2obj->c_B_prev, f2obj->c_B);
		}
	}

	if (Parms.p_Do_Dens_Disp)
	{
		// DISSIPATE Density by multipler - 
		dissipate(f2obj->dens, Parms.p_Dens_Disp_Mult, this->dt); // Density Dissipation. 
	}

}
// End of Density Step Implemetnation.

/*	====================================================
	Velocity SOLVE STEP -
	==================================================== */

// Implementation of Simulation Solve step of Velocity solve operations. 
// Excluding Sourcing which will be done via user in FluidObject.
// Removed Pre Advection Projection Calls, So I can use One Call Post Advection, With Higher Iter Counts.

void fluidsolver_2::velocity_step(int diff_iter, int proj_iter, float diffA, bool dodiff)
{
	// If Using Diffusion, If Not Need to Manually set Cur to Prev, rather than swapping. 
	if (dodiff == true)
	{
		f2obj->vel->swap(f2obj->prev_vel); // Swap Vel Field with Prev_VelField.
		diffuse(f2obj->prev_vel, f2obj->vel, diffA, diff_iter);
		f2obj->vel->swap(f2obj->prev_vel); // Re-Swap Vel With Prev_Vel. 
	}
	else if (dodiff == false)
	{
		f2obj->setcurtoprev_vel();
	}

	// ADVECT VELOCITY FIELD (Self Advect) \\

	if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_Euler)
	{
		advect_sl(f2obj->prev_vel, f2obj->vel);
	}
	else if (Parms.p_AdvectionType == Parms.Advect_SL_BackTrace_RK2)
	{
		advect_sl_mp(f2obj->prev_vel, f2obj->vel);
		//advect_sl_mp_GS(f2obj->prev_vel, f2obj->vel);
	}

	if (Parms.p_Do_Vel_Disp)
	{
		// DISSIPATE Density by multipler - 
		dissipate(f2obj->vel, Parms.p_Vel_Disp_Mult, this->dt); // Density Dissipation. 
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

// Implementation of A Single solve Step of fluidsolver_2. Once User Calls this
// Simulation will begin. And for now occur Untill program is closed. 

void fluidsolver_2::solve_step(bool solve, bool do_diffdens, bool do_diffvel, float dens_diff, float vel_diff, int proj_iter, int diff_iter, int max_step)
{
	// Init Step/time vals. 
	int step_count = 0;
	float dt_acc = 0.0f; // Accumulated DeltaTime (SimTime)
	float total_elapsed = 0.0f; // Chrono Time

	// Render Object Creation/Setup - 
	// Pass Render Context Window Ptr Now FluidSovler2Mem, to RenderObject. (Passed in To Solver in Main via set_window() MF.
	int render_mode = 1; // 0 = Density, 1 = Vel. // May be overriden by Input in SolveStep (RenderObj::ShaderPipe()).
	render_obj = new renderobject_2D_OGL("OpenGL", 4, 0, x_s + e_s, y_s + e_s, this->winptr, render_mode); 
	// Ideally Move render_obj setup to main/outside, and then pass into fluidsolver...?

	// Move this into RenderObject Initalization - 
	// Set Render_Object Texture Units to Sampler Uniforms Pre Solve Loop - 
	glUseProgram(render_obj->shader_prog); // Call Use (Shader) Program First. 

	// Set Tex Sampler Uniforms, to matching Texture Units.
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "d_tex"), 0); // Density = 0. 
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "v_u_tex"), 1); // Velocity u = 1. 
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "v_v_tex"), 2); // Velocity v = 2. 
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "c_tex"), 3); // Collision c = 3. 
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "vc_u_tex"), 4); // VortConfine u = 4. 
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "vc_v_tex"), 5); // VortConfine v = 5. 
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "img_rgb_tex"), 6); // Image RGB = 6.
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "ppv_u_tex"), 7); // PreProj Vel U = 7.
	glUniform1i(glGetUniformLocation(render_obj->shader_prog, "ppv_v_tex"), 8); // PreProj Vel V = 8.

	//Solve Step And Render - EXEC LOOP 
	while (solve == true && step_count <= max_step) // Infinite Solve Loop For Now. Will need to link this to drawing/waiting etc. 
	{
		// Time Start - 
		std::chrono::system_clock::time_point timestep_start = std::chrono::system_clock::now();
		std::stringstream log_out;

		// STEP INPUT OPERATIONS \\ ----------------
		// Interactive Sourcing ... (Implemented with GL WIP). Do Here vs Inside RenderObject (passed grids). Do they have scope withiN?
		// Interactive SphereBounds Collider Offset from Cursor Pos.
		// TESTING - DEC 19 
		// Check for Input to Modifiy SphereBounds Radius. 
		sphere_rad_test(); // Can Cause Pressure Solver Crashes. 
		if (glfwGetKey(winptr, GLFW_KEY_R) == GLFW_PRESS) // Messy inline for now. 
		{
			// Re-init Texture-ColourGrids - 
			f2obj->vel->clear(); // Also Clear Vel Field. 
			f2obj->RGB_imageLoad(f2obj->img_data.path);
		}
		// Interp Switching.
		//if (glfwGetKey(winptr, GLFW_KEY_I) == GLFW_PRESS) { Parms.p_InteroplationType == Parms.Interoplation_Linear; };
		//if (glfwGetKey(winptr, GLFW_KEY_C) == GLFW_PRESS) { Parms.p_InteroplationType == Parms.Interoplation_Cosine; };

		// Get CurFrame Mouse Pos And Update Mouse Vel. 
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
		//vec2 offset(-0.5, -0.85);
		//anim_offset.x = xpos; anim_offset.y = ypos; // Use Cursor Pos
		//vec2 anim_offset(0.4 + (float(sin(float(step_count) / float(max_step) * 50.0f))) * 0.1, 0.1); 
		// Fun Times Below - 
	//	vec2 anim_offset(0.4 + (float(sin(float(step_count) / float(max_step) * 50.0f))) * 0.1f, 0.4 + (float(cos(float(step_count) / float(max_step) * (float(step_count) / float(max_step) * 10.0f)))) * 0.2f);
	//	f2obj->implicit_source(0.1f, vec2(0.0f, 0.1f), anim_offset, 0.002f);
		//f2obj->implicit_source(0.25f, vec2(0.0f, 0.1f), vec2(0.5f, 0.1f), 0.025f); 

		f2obj->implicit_source(0.2f, vec2<float>(0.0f, 0.25f), vec2<float>(0.5f, 0.2f), 0.01f);

	//	f2obj->implicit_source(0.25f, vec2<float>(0.0f, 0.0f), vec2<float>(0.5f, 0.5f), 0.1f);

		// Forces- 
	//	if (step_count <= 20) f2obj->radial_force(vec2<float>(0.499f, 0.499f), 0.8f, this->dt);

		// STEP - SUB - SOLVER STEP OPERATIONS \\ -------------- 
		velocity_step(diff_iter, proj_iter, vel_diff, do_diffvel);
		density_step(diff_iter, dens_diff, do_diffdens);

		// STEP - RENDER CALLS \\ ------------------
		#if RENDER_GL == 1
		// Pass Cur Step - 
		render_obj->et = step_count; 
		// Uniform/Texture Set Calls (ShaderPipe) -
		render_obj->shader_pipe(f2obj); // Pass Grids Via Friend Acess to shader_pipe() MF
		// Render Operations Call - 
		render_obj->call_ren(rend_state::RENDER_ACTIVE); // Call Render Step Ops here within solver loop, ie NON Debug Mode (Pass RENDER_ACTIVE).
		#endif	

		#if RENDER_IMG == 1
		// STEP WRITE OUTPUT OPERATIONS - Non OGL Img output.
		//f2obj->writeto_img_all(step_count, pressure);
		f2obj->writeto_img(step_count); // Write Density Only, Debug.
	//	f2obj->writeto_img_vel(step_count); // Write Velocity Only, Debug. 
	//	f2obj->writeto_txt(step_count);
		#endif

		// STEP CONSLE DEBUG OPERATIONS \\ -------------------- 
		std::chrono::system_clock::time_point timestep_end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed = timestep_end - timestep_start;
		double elapsed_sec = elapsed.count();
		total_elapsed += elapsed_sec;

		// Store Mouse Pos as Prev Frame Pos.. 
		updt_mousepos(step::STEP_PREV); updt_mouseposNorm(step::STEP_PREV); updt_mouseposRange(step::STEP_PREV);

		// Differintate between SimTime Passed of accumlated DeltaTime, Vs Duration of SolveStep Caluclation !
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
			//f2obj->writeto_img(step_count);
			//std::cout << "\n"; f2obj->print_info();
			//f2obj->writeto_img_vel(step_count);
		}

		// Print Log_output - 
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
void fluidsolver_2::set_window(GLFWwindow *win)
{
	assert(win != nullptr); // Check for Passing a nullptr FluidSolver GLFW WinPtr. Assert this.
	this->winptr = win; 
}

// Call To Update Mouse Pos - (Pixel (Cell Index) Space)
void fluidsolver_2::updt_mousepos(const step step_id)
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
void fluidsolver_2::updt_mouseposNorm(const step step_id)
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

void fluidsolver_2::updt_mouseposRange(const step step_id)
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
void fluidsolver_2::updt_mousevel()
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
void fluidsolver_2::dissipate(grid2_scalar *grid, float disp_mult, float dt)
{
	//disp_mult = std::max(0.0f, std::min(disp_mult, 1.0f)); // Enforce 0-1 Mult. 
	disp_mult = solver_utils::clamp(disp_mult, 0.0f, 1.0f); 

	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Don't Mult By dt for now. 
			float cur_scl = grid->getdata(i, j);
			cur_scl *= disp_mult;
			grid->setdata(cur_scl, i, j);
		}
	}
}

void fluidsolver_2::dissipate(grid2_vector *grid, float disp_mult, float dt)
{
	disp_mult = solver_utils::clamp(disp_mult, 0.0f, 1.0f);
	for (int j = 1; j <= N_dim; j++)
	{
		for (int i = 1; i <= N_dim; i++)
		{
			// Don't Mult By dt for now. 
			vec2<float> cur_vel = grid->getdata(i, j);
			cur_vel *= disp_mult;
			grid->setdata(cur_vel, i, j);
		}
	}
}



/*	====================================================
Solver Utils Implementation (Lambda Utility Funcs) -
==================================================== */

// Implement Static Function Lambdas. Why not just use Member Funcs? Well I'm been fancy. 

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

// BILINEAR: 2D Linear Interoplation Of Cell i,j by Coefficents a,b. 
// Scalar Grid Input - (Cannot use Generics, due to inner Type diffs)
std::function<float(grid2_scalar*, float, float, int, int)> solver_utils::bilinear_S = [&](grid2_scalar *grid, float a, float b, int i, int j) -> float
{
	// Complement to get other coeffs. 
	float a1 = 1 - a; float b1 = 1 - b;

	// LERP X,Y(b,b1) (SameCoeffs b,b1 because sqaure grid) and between resuilting X,Y (by a,a1). 

	return a * (b * grid->getdata(i, j) + b1 * grid->getdata(i, j + 1)) +
		a1 * (b * grid->getdata(i + 1, j) + b1 * grid->getdata(i + 1, j + 1));
};
// Vector Grid Input - (short specfices x or y component) 
std::function<float(grid2_vector*, float, float, int, int, ushort)> solver_utils::bilinear_V = [&](grid2_vector *grid, float a, float b, int i, int j, ushort mode) -> float
{
	assert(mode <= 1); // 0 = X, 1 = Y.

	// Complement to get other coeffs. 
	float a1 = 1 - a; float b1 = 1 - b;

	// LERP X,Y(b,b1) (SameCoeffs b,b1 because sqaure grid) and between resuilting X,Y (by a,a1). 

	if (mode == 0)
	{
		return (float) a * (b * grid->getdata_x(i, j) + b1 * grid->getdata_x(i, j + 1)) +
			a1 * (b * grid->getdata_x(i + 1, j) + b1 * grid->getdata_x(i + 1, j + 1));
	}
	else if (mode == 1)
	{
		return (float) a * (b * grid->getdata_y(i, j) + b1 * grid->getdata_y(i, j + 1)) +
			a1 * (b * grid->getdata_y(i + 1, j) + b1 * grid->getdata_y(i + 1, j + 1));
	}

	return 0.0f; 
};

// BICUBIC
// Scalar - 
std::function<float(grid2_scalar*, float, float, float, float, int, int)> solver_utils::bicubic_S = [&]
	(grid2_scalar *grid,
	float wn1, float w0, float w1, float w2,
	int i, int j) -> float
{

	float qn1 = (wn1 * grid->getdata(i - 1, j - 1)) + (w0 * grid->getdata(i, j - 1)) + (w1 * grid->getdata(i + 1, j - 1)) + (w2 * grid->getdata(i - 2, j - 1));
	float q0 =  (wn1 * grid->getdata(i - 1, j)) + (w0 * grid->getdata(i, j)) + (w1 * grid->getdata(i + 1, j)) + (w2 * grid->getdata(i - 2, j));
	float q1 =  (wn1 * grid->getdata(i - 1, j + 1)) + (w0 * grid->getdata(i, j + 1)) + (w1 * grid->getdata(i + 1, j + 1)) + (w2 * grid->getdata(i - 2, j + 1));
	float q2 =  (wn1 * grid->getdata(i - 1, j + 2)) + (w0 * grid->getdata(i, j + 2)) + (w1 * grid->getdata(i + 1, j + 2)) + (w2 * grid->getdata(i - 2, j + 2));

	return (float) (wn1 * qn1) + (w0 * q0) + (w1 * q1) + (w2 * q2); 
};

// Vector - 
std::function<float(grid2_vector*, float, float, float, float, int, int, ushort)> solver_utils::bicubic_V = [&]
	(grid2_vector *grid,
	float wn1, float w0, float w1, float w2,
	int i, int j, ushort mode) -> float
{
	assert(mode <= 1); // 0 = X, 1 = Y.

	float qn1 = 0.0f, q0 = 0.0f, q1 = 0.0f, q2 = 0.0f; 

	if (mode == 0) // Vector.x Value
	{
		qn1 = (wn1 * grid->getdata_x(i - 1, j - 1)) + (w0 * grid->getdata_x(i, j - 1)) + (w1 * grid->getdata_x(i + 1, j - 1)) + (w2 * grid->getdata_x(i - 2, j - 1));
		q0 = (wn1 * grid->getdata_x(i - 1, j)) + (w0 * grid->getdata_x(i, j)) + (w1 * grid->getdata_x(i + 1, j)) + (w2 * grid->getdata_x(i - 2, j));
		q1 = (wn1 * grid->getdata_x(i - 1, j + 1)) + (w0 * grid->getdata_x(i, j + 1)) + (w1 * grid->getdata_x(i + 1, j + 1)) + (w2 * grid->getdata_x(i - 2, j + 1));
		q2 = (wn1 * grid->getdata_x(i - 1, j + 2)) + (w0 * grid->getdata_x(i, j + 2)) + (w1 * grid->getdata_x(i + 1, j + 2)) + (w2 * grid->getdata_x(i - 2, j + 2));
	}
	else if (mode == 1) // Vector.y Value
	{
		qn1 = (wn1 * grid->getdata_y(i - 1, j - 1)) + (w0 * grid->getdata_y(i, j - 1)) + (w1 * grid->getdata_y(i + 1, j - 1)) + (w2 * grid->getdata_y(i - 2, j - 1));
		q0 = (wn1 * grid->getdata_y(i - 1, j)) + (w0 * grid->getdata_y(i, j)) + (w1 * grid->getdata_y(i + 1, j)) + (w2 * grid->getdata_y(i - 2, j));
		q1 = (wn1 * grid->getdata_y(i - 1, j + 1)) + (w0 * grid->getdata_y(i, j + 1)) + (w1 * grid->getdata_y(i + 1, j + 1)) + (w2 * grid->getdata_y(i - 2, j + 1));
		q2 = (wn1 * grid->getdata_y(i - 1, j + 2)) + (w0 * grid->getdata_y(i, j + 2)) + (w1 * grid->getdata_y(i + 1, j + 2)) + (w2 * grid->getdata_y(i - 2, j + 2));
	}

	return (float)(wn1 * qn1) + (w0 * q0) + (w1 * q1) + (w2 * q2);
};

/*
qj-1 = W-1 i-1,j-1 + W0 i,j-1 + w1 i+1,j-1 + W2 i+2,j-1 // j Held Constant j-1 (Solve for x)
qj   = w-1 i-1,j   + w0 i,j   + w1 i+1,j   + w2 i+2,j   // j
qj+1 = w-1 i-1,j+1 + w0 i,j+1 + w1 i+1,j+1 + w2 i+2,j+1 // j+1
qj+2 = w-1 i-1,j+2 + w0 i,j+2 + w1 i+1,j+2 + w2 i+2,j+2 // j+2
q = w-1 qj-1 + w0 qj + w1 qj+1 + w2 qj+2
*/



/*	====================================================
	DEBUG - Member FUnctions 
	==================================================== */

void fluidsolver_2::fill_test(int x)
{
	if (x >= NE_dim) x = NE_dim; 

	for (int j = 0; j < x; j++)
	{
		for (int i = 0; i < x; i++)
		{
			f2obj->dens->setdata(1.0f, i, j);
		}
	}
}

// Test Implementation of Interactive Sphere Radius - Can Cause Pressure Solver Crashes. 
void fluidsolver_2::sphere_rad_test()
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
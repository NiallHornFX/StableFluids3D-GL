// Implementation of fluidobj
// std depends/includes via Header.

//#include "fluidobj.h"
#include <fstream> 
#include <string>
#include <sstream>
#include <chrono>
#include <cassert>

#include <omp.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define cur2d "[" << i << "," << j << "]  " // For Printing.

extern short verbose; // Get Verbose Global from main. 

// fluidobj_2d Implemenatation - 

// fluidobj_2d Constructor Main - 

fluidobj_2d::fluidobj_2d(int x_size, int y_size, int edge_size, int idd, float spc) 
	: x_s(x_size), y_s(y_size), e_s(edge_size), ID(idd), spacing(spc)
{
	// Heap Allocate Grids to Member Pointers, Manual ID Paramters. 
	dens = new grid2_scalar(x_size, y_size, edge_size, 1, spc);
	prev_dens = new grid2_scalar(x_size, y_size, edge_size, 0, spc);

	vel = new grid2_vector(x_size, y_size, edge_size, 3, spc);
	prev_vel = new grid2_vector(x_size, y_size, edge_size, 2, spc);
	col = new grid2_scalar(x_size, y_size, edge_size, 4, spc);

	curl = new grid2_scalar(x_size, y_size, edge_size, 5, spc);
	vc = new grid2_vector(x_size, y_size, edge_size, 6, spc);

	c_R = new grid2_scalar(x_size, y_size, edge_size, 10, spc);
	c_G = new grid2_scalar(x_size, y_size, edge_size, 11, spc);
	c_B = new grid2_scalar(x_size, y_size, edge_size, 12, spc);
	c_R_prev = new grid2_scalar(x_size, y_size, edge_size, 13, spc);
	c_G_prev = new grid2_scalar(x_size, y_size, edge_size, 14, spc);
	c_B_prev = new grid2_scalar(x_size, y_size, edge_size, 15, spc);

	// Debug Grids. 
	preproj_vel = new grid2_vector(x_size, y_size, edge_size, 16, spc);

}


fluidobj_2d::~fluidobj_2d()
{
	// Dealloc Grids Heap - 

	delete dens; dens = nullptr; 
	delete prev_dens; prev_dens = nullptr;
	delete vel; vel = nullptr;
	delete prev_vel; prev_vel = nullptr; 
	delete col; col = nullptr; 

	// If Vorticity Confinement Temp Grids Still Allocated - 
	if (curl) delete curl, curl = nullptr; 
	if (vc) delete vc, vc = nullptr; 

	// If Colour Grids were used/allocated. 
	if (c_R) 
	{
		delete c_R, c_R = nullptr; 
		delete c_G, c_G = nullptr; 
		delete c_B, c_B = nullptr; 
		delete c_R_prev, c_R_prev = nullptr; 
		delete c_G_prev, c_G_prev = nullptr; 
		delete c_B_prev, c_B_prev = nullptr; 
	}

	// If Debug Grids were used/allocated. 
	if (preproj_vel) delete preproj_vel, preproj_vel = nullptr; 

	// Free Image Buffer
	stbi_image_free(img_data.img_raw);


	// Check for Deletion. 
	/*
	if (!dens && !prev_dens && !vel && !prev_vel)
	{
		std::cout << "DEBUG::fluid_obj " << ID << " Destructor Called Sucessfully. \n";
	}
	*/
}

// Add Density/Velocity Member Functions - Pass Data to Grid_Data std::vector elements
// Using 2D Indices, via the grid2_*s setdata Member Funcs. 
// These Can be Done Per step, or just once on the inital frame, before the sim loop. 

// Pass In Float Density at Given Grid Location (2D Index) 
void fluidobj_2d::add_density(float dens, int i, int j)
{
	//this->dens->setdata(dens, i, j);
	this->dens->adddata(dens, i, j);
	// Oops use this because paramter is same as dens grid member name. 
}

// Pass In Vec2 Velocity at Given Grid Location (2D Index) 
void fluidobj_2d::add_velocity(const vec2<float> &vel, int i, int j)
{
	this->vel->adddata(vel, i, j);
	// Oops use this because paramter is same as dens grid member name. 
}

// Add Force To Velocity. Assume Mass/Density 1.0f. 
void fluidobj_2d::integrate_force(const vec2<float> &force, float dt, int i, int j)
{
	// Needs to be called per frame, to intergrate force to vel. 
	vec2<float> temp = force;
	temp *= dt; 
	vel->adddata(temp, i, j);
}

void fluidobj_2d::print_info()
{
	// Calc Size of Scalar and Vec2 Grids. Via (X+edge*Y+edge)
	//std::size_t sgrid_bytes = sizeof(float) * ((x_s + e_s) * (y_s + e_s));
	//std::size_t vgrid_bytes = sizeof(vec2) * ((x_s + e_s) * (y_s + e_s));

	// Just use sizeof (T) * GridData Size - 
	std::size_t sgrid_bytes = sizeof (float) * dens->griddataptr_getter()->size(); 
	std::size_t vgrid_bytes = sizeof (vec2<float>) * vel->griddataptr_getter()->size(); 

	// Correctly using 1024 B-KB-MB-GB Ranges. 
	std::size_t sgrid_KB = sgrid_bytes / 1024; // 1024 Btyes Per KB. 
	std::size_t sgrid_MB = sgrid_bytes / int(pow(1024, 2)); // 1048576 (1.048 Million) Btyes Per MB. 
	std::size_t sgrid_GB = sgrid_bytes / int(pow(1024, 3)); // 1073741824 (1.073 Billion) Btyes Per GB. 

	std::size_t vgrid_KB = vgrid_bytes / 1024; // 1024 Btyes Per KB. 
	std::size_t vgrid_MB = vgrid_bytes / int(pow(1024, 2)); // 1048576 (1.048 Million) Btyes Per MB. 
	std::size_t vgrid_GB = vgrid_bytes / int(pow(1024, 3)); // 1073741824 (1.073 Billion) Btyes Per GB. 

	// Print Memory of Fluid Object Density (Scalar) and Velocity (Vector) Grids Using FluidObject Dimension Vals -

	std::cout << "\nDEBUG::Fluid Object 2D " << ID << " INFO BEGIN - \n";
	// INFO Code - 
	std::cout << "DEBUG::Fluid Object 2D:: Scalar Grid Size = " << sgrid_bytes << " B \n";
	std::cout << "DEBUG::Fluid Object 2D:: Scalar Grid Size = " << sgrid_KB << " KB \n";
	std::cout << "DEBUG::Fluid Object 2D:: Scalar Grid Size = " << sgrid_MB << " MB \n";
	std::cout << "DEBUG::Fluid Object 2D:: Scalar Grid Size = " << sgrid_GB << " GB \n";
	std::cout << "DEBUG::Fluid Object 2D:: Scalar Cell Size = " << sizeof(float) << " B \n";
	std::cout << "============================================ \n";
	std::cout << "DEBUG::Fluid Object 2D:: vec2 Grid Size = " << vgrid_bytes << " B \n";
	std::cout << "DEBUG::Fluid Object 2D:: vec2 Grid Size = " << vgrid_KB << " KB \n";
	std::cout << "DEBUG::Fluid Object 2D:: vec2 Grid Size = " << vgrid_MB << " MB \n";
	std::cout << "DEBUG::Fluid Object 2D:: vec2 Grid Size = " << vgrid_GB << " GB \n";
	std::cout << "DEBUG::Fluid Object 2DLL vec2 Cell Size = " << sizeof(vec2<float>) << " B \n";
	std::cout << "DEBUG::Fluid Object 2D " << ID << " INFO END. \n \n";

	// Print Size of Total FluidObj Allocated grids memory. Based on above vec/scalar grid size. 
}

// Grid Member Delete MFs - 
// These probs will never be used, as these should only die on program/sim completion. 

void fluidobj_2d::del_dens()
{

	// Clear Grid Heap Data through mem ptrs. 
	delete dens; delete prev_dens;
	//std::cout << "DEBUG::fluid_obj " << ID << " Density Grids Dealloacted Sucesfully. \n";
	dens = nullptr; prev_dens = nullptr;

}

void fluidobj_2d::del_vel()
{
	// Clear Grid Heap Data through mem ptrs. 
	delete vel; delete prev_vel;
	//std::cout << "DEBUG::fluid_obj " << ID << " Velocity Grids Dealloacted Sucesfully. \n";
	vel = nullptr; prev_vel = nullptr;
}


// Input Handling GLFW Window Based MFs (Sourcing Coming Soon...) -
// Simple Sourcing For Now - 

// Implicit SDF Circle For Sourcing into Density and Velocity Grids. 

void fluidobj_2d::implicit_source(float dd, const vec2<float> &vv, const vec2<float> &offset, float rad)
{
	// Iterate Over Grid Cells, Test for Implicit Condtion. 
	// if (pow(@P.x, 2) * 10 <= -@P.y) where Px and Py = i,j.

	int N = x_s; // One Dim Size (Cell Count N). X,Y are Same ofc square grid.
	float h = 1.0f / N; // Recoprical of One Dim Size (N). 
	float surf_isothresh = 0.01f; 

	//#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int j = 1; j <= N; j++)
	{
		for (int i = 1; i <= N; i++)
		{
			// Implicit Sphere/Circle Function: x^2 + y^2 - r . >= 0 <= Thresh == Surface. > Thresh = Exterior. < 0 = Interior. Cells. 
			vec2<float> cell_gridSpace(float((i * h) - offset.x), float((j * h) - offset.y)); // Index to Grid Space 0-1N. 
			float sphere_func = ((cell_gridSpace.x * cell_gridSpace.x) + (cell_gridSpace.y * cell_gridSpace.y)) - rad; // NoMoreSqrt.

			// ALL Cells Sourcing - 

			// INSIDE Cells Sourcing - 
			if (sphere_func < 0.0f)
			{
				add_density(dd, i, j);
				add_velocity(vv, i, j);
			}

			// SURFACE Cells Sourcing -
			if (sphere_func >= 0.0f && sphere_func <= surf_isothresh)
			{
				// 
			}

			// EXTERIOR Cells Sourcing - 
			if (sphere_func > surf_isothresh)
			{
				//
			}

		}
	}
}

// Velocity Field Forces \\ 

// Radial Force - 
void fluidobj_2d::radial_force(const vec2<float> &orig_offset, float strength, float dt)
{
	// Thread Safe - 

	#pragma omp parallel for
	for (int j = 1; j <= y_s; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= x_s; i++)
		{
			// Index --> GridSpace Vector (x, y); 
			vec2<float> gs((float)i / x_s, (float)j / y_s);
			// Offset GridSpace Origin (x - h, y - k) (Offset (h,k) Passed in Cartesian GridSpace)
			gs -= orig_offset;

			// Set Radial Velocity, Normalize and Mult Speed Length. 
			vec2<float> rad_vel(gs.y - gs.x, -gs.y - gs.x);
			rad_vel.normalize() *= strength;

			// Integrate Force to Velocity - 
			integrate_force(rad_vel, dt, i, j);
		}
	}
}


// Inital Velocity Condtion (Velocity Override) Member Functions \\

// These Override Velocity, and Should Only be added Initally PreSolveStep call. Not within.

// Radial Velocity - 
void fluidobj_2d::radial_vel(const vec2<float> &orig_offset, float speed)
{
	for (int j = 1; j <= y_s; j++)
	{
		for (int i = 1; i <= x_s; i++)
		{
			// Index --> GridSpace Vector (x, y); 
			vec2<float> gs((float)i / x_s, (float)j / y_s);
			// Offset GridSpace Origin (x - h, y - k) (Offset (h,k) Passed in Cartesian GridSpace)
			gs -= orig_offset;

			// Set Radial Velocity, Normalize and Mult Speed Length. 
			vec2<float> rad_vel(gs.y - gs.x, -gs.y - gs.x);
			rad_vel.normalize() *= speed; 

			// Set Velocity In All Cells -
			vel->setdata(rad_vel, i, j);
		}
	}
}

// Sink Velocity (To Some GridSpace OriginOffset / Postion) - 

// (0.5f, 0.5f) Offset Vector Crashes? (x - 0.5, y - 0.5) (Offset Origin to Center) ? 
// (0.5f, 0.4999f) works... Must be rounding error. 
void fluidobj_2d::sink_vel(const vec2<float> &orig_offset, float speed)
{
	for (int j = 1; j <= y_s; j++)
	{
		for (int i = 1; i <= x_s; i++)
		{
			// Index --> GridSpace Vector (x, y); 
			vec2<float> gs((float)i / x_s, (float)j / y_s);

			// Offset GridSpace Origin (x - h, y - k) (Offset (h,k) Passed in Cartesian GridSpace)
			gs -= orig_offset;

			// Re-Normalize and Negate Grid Space Pos Vector to Face Offseted Origin. 
			gs.normalize(); gs.x *= -1, gs.y *= -1;

			// Set Velocity In All Cells -
			vel->setdata((gs *= speed), i, j);
		}
	}
}

// Write Grid Data To Image MFunctions \\ 

// Write to Image - (Per Frame/Step) 

// Add Code to set output res, and interoplate cell to upres image output res. 
int fluidobj_2d::writeto_img(int frame)
{
	// With Edge Cells - Edge Pixels.
	int w = x_s + e_s; 
	int h = y_s + e_s;

	// Doing Path via StringStream doesnt work for some reason. Will just use string. 
	//std::stringstream pathss;
	//pathss << "bin/img/outimg_" << frame << ".ppm" << std::endl;
	//std::string paths = pathss.str();
	//const char *pathc = paths.c_str();

	// Posiblly Encase in try/catch for expection safety output? With ofstream badbit/failbut except enabled. 

	// Out Path String (paths).
	// Append Time For Per Run Imgs. 
	auto time = std::chrono::system_clock::now();
	auto timet = std::chrono::system_clock::to_time_t(time);

	std::string paths = "build/img/seq/fluidimg_";
	//paths.append(std::to_string(timet) + "_");
	paths.append("_");
	paths.append(std::to_string(frame) + ".ppm");
	if (verbose) {
		std::cout << "DEBUG::WriteToImage_Output: " << paths << std::endl;
	}

	// Out Image Write OFilestream (img_out).
	std::ofstream img_out(paths); // .open is called on construction (OFStream). 
	img_out << "P3" << std::endl;
	img_out << w << " " << h << std::endl;
	img_out << "255" << std::endl;

	// 0 Based Iter, Output edge Cells/Pixels. 
	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			float value = dens->getdata(i, h-j); // FLIPPED Y Axis (hence h-j (N_Size - j) Backwards Y Iteration So UP is Postive Y/j. 
			//float uval = vel->getdata_x(i, j);
			//float vval = vel->getdata_y(i, j);
			float colv = col->getdata(i, j);

			int r = int(value * 255);
			int g = int(value * 255);
			int b = int(value * 255);

			//Visible Edges - 
			/*
			if (j == 0) { r = 20, g = 0, b = 0; }
			if (j == (h-1)) { r = 20, g = 0, b = 0; }
			if (i == 0) { r = 20, g = 0, b = 0; }
			if (i == (h-1)) { r = 20, g = 0, b = 0; }
			*/

			// Clamp Min/Max. 
			if (r <= 1) r = 1; if (r >= 255) r = 255;
			if (g <= 1) g = 1; if (g >= 255) g = 255;
			if (b <= 1) b = 1; if (b >= 255) b = 255;

			// Check for -Y Vel. 
			if (colv != 0.0f)
			{
				//g = 0; b = 0; // Color Cell->PixelValue Red. 
			}
			
			img_out << r << " " << g << " " << b << std::endl;
		}
	}

	return 1; // If Sucessful Write - WIP -  Need to encap try catch and return 0 if expects thrown from OFStream etc.. 
}

/*
// Write to Image B WIP - (Per Frame/Step) - WIP - Use Arbitary Width/Height Interpolate from grid Size...  NearestNeighbour Upscale...

// Add Code to set output res, and interoplate cell to upres image output res. 
void fluidobj_2d::writeto_img_b(int frame, int width, int height)
{
	// With Edge Cells - Edge Pixels.
	int w = x_s + e_s;
	int h = y_s + e_s;

	// Doing Path via StringStream doesnt work for some reason. Will just use string. 
	//std::stringstream pathss;
	//pathss << "bin/img/outimg_" << frame << ".ppm" << std::endl;
	//std::string paths = pathss.str();
	//const char *pathc = paths.c_str();

	// Posiblly Encase in try/catch for expection safety output? With ofstream badbit/failbut except enabled. 

	// Out Path String (paths).
	std::string paths = "bin/img/foo_";
	paths.append(std::to_string(frame) + ".ppm");
	std::cout << paths << std::endl;

	// Out Image Write OFilestream (img_out).
	std::ofstream img_out(paths);
	img_out << "P3" << std::endl;
	img_out << w << " " << h << std::endl;
	img_out << "255" << std::endl;

	// if i j are within grid range... else fillout with neighbour/0 values... WIP
	// 0 Based Iter, Output edge Cells/Pixels. 
	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			float value = dens->getdata(i, j);
			int r = int(value * 255);
			int g = int(value * 255);
			int b = int(value * 255);

			// Clamp Min/Max. 
			if (r <= 1) r = 1; if (r >= 255) r = 255;
			if (g <= 1) g = 1; if (g >= 255) g = 255;
			if (b <= 1) b = 1; if (b >= 255) b = 255;

			img_out << r << " " << g << " " << b << std::endl;
		}
	}

}
*/

// TEST - Wrtie To Img (Velocity Grid)
int fluidobj_2d::writeto_img_vel(int frame)
{
	// With Edge Cells - Edge Pixels.
	int w = x_s + e_s;
	int h = y_s + e_s;

	// Out Path String (paths).
	std::string paths = "build/img/seq/velimg_";
	paths.append(std::to_string(frame) + ".ppm");
	if (verbose)
	{
		std::cout << "DEBUG::WriteToImage_Velocity_Output: " << paths << std::endl;
	}
	
	// Out Image Write OFilestream (img_out).
	std::ofstream img_out(paths); // .open is called on construction (OFStream). 
	img_out << "P3" << std::endl;
	img_out << w << " " << h << std::endl;
	img_out << "255" << std::endl;

	// 0 Based Iter, Output edge Cells/Pixels. 
	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			// RGB Values Derive from Velocity XY Compoennts (B is 0 in 2D).
			float value_r = vel->getdata_x(i, j);
			float value_g = vel->getdata_y(i, j);
			float value_b = 0.0f; 

			int r = int(value_r * 255);
			int g = int(value_g * 255);
			int b = int(value_b * 255);

			if (r <= 1) r = 1; if (r >= 255) r = 255;
			if (g <= 1) g = 1; if (g >= 255) g = 255;
			if (b <= 1) b = 1; if (b >= 255) b = 255;

			img_out << r << " " << g << " " << b << std::endl;
		}
	}

	return 1; // If Sucessful Write - WIP -  Need to encap try catch and return 0 if expects thrown from OFStream etc.. 
}

// TEST - Write Pressure Field to Img, If FluidSolver has Kept it in memory.  Then Delete Pressure Field For Current Step.
// (To Avoid Keeping Each Steps Pressure Field (thus memory leak).
int fluidobj_2d::writeto_img_pressure(int frame, bool pressureexist, grid2_scalar *pressureptr)
{
	if (pressureexist == false)
	{
		return 0;
		// Early return, because Pressure Grid has been deleted by solver already. 
	}

	// With Edge Cells - Edge Pixels.
	int w = x_s + e_s;
	int h = y_s + e_s;

	// Out Path String (paths).
	std::string paths = "build/img/pressureimg_";
	paths.append(std::to_string(frame) + ".ppm");
	if (verbose)
	{
		std::cout << "DEBUG::WriteToImage_Pressure_Output: " << paths << std::endl;
	}
	
	// Out Image Write OFilestream (img_out).
	std::ofstream img_out(paths); // .open is called on construction (OFStream). 
	img_out << "P3" << std::endl;
	img_out << w << " " << h << std::endl;
	img_out << "255" << std::endl;

	// 0 Based Iter, Output edge Cells/Pixels. 
	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			// Acess Pressure Grid, through Passed In grid2_scalar ptr from FluidSolver obj.
			float value = pressureptr->getdata(i, j);
			value *= 2.5f; // Brightness Mult. 

			int r = int(value * 255);
			int g = int(value * 255);
			int b = int(value * 255);

			if (r <= 1) r = 1; if (r >= 255) r = 255;
			if (g <= 1) g = 1; if (g >= 255) g = 255;
			if (b <= 1) b = 1; if (b >= 255) b = 255;

			img_out << r << " " << g << " " << b << std::endl;
		}
	}

	/* This is now done by the FluidSolver itself, at the end of the timestep
	// Once Pressure is written to img, I DELETE/DEALLOC The Pressure Grid Here (so the next step pressure gird pointer
	//does not leave this steps pointer dangling. So I delete it in the scope of this func, via the passed pointer. 

	delete pressureptr;
	pressureptr = nullptr; // Set Passed Ptr to NULL So when this func is popped, Deallocation of passed pointer is not an issue. 
	*/ 

	return 1; // If Sucessful Write - WIP -  Need to encap try catch and return 0 if expects thrown from OFStream etc.. 
}
// End of fluidobj2::writeto_img_pressure implementation. 

// Write To Image All Implementation -

// Oppsoe to writing grid data out to file/ppm within the loop, pass all to a stringstream, that we then write out
// externally to avoid consant overhead of writing to ofstream per i*j iter. 
// Also ALL Grids are wrote out here, oppose to sepreate calls to writeto MFuncs, which is slow. 

int fluidobj_2d::writeto_img_all(int frame, grid2_scalar *pres_grid)
{
	// With Edge Cells - Edge Pixels.
	int w = x_s + e_s;
	int h = y_s + e_s;

	// Out Path String (paths).
	std::string path_o_dens = "build/img/img-dens_";
	path_o_dens.append(std::to_string(frame) + ".ppm");
	
	std::string path_o_vel = "build/img/img-vel_";
	path_o_vel.append(std::to_string(frame) + ".ppm");

	std::string path_o_pres = "build/img/img-pres_";
	path_o_pres.append(std::to_string(frame) + ".ppm");

	if (verbose)
	{
		std::cout << "DEBUG::WriteToImage_ALL_Output: " << path_o_dens << std::endl;
		std::cout << "DEBUG::WriteToImage_ALL_Output: " << path_o_vel << std::endl;
		std::cout << "DEBUG::WriteToImage_Output: " << path_o_pres << std::endl;
	}


	// Out Image Write OFilestream (img_out).
	std::ofstream img_out_dens(path_o_dens); // .open is called on construction (OFStream). 
	std::ofstream img_out_vel(path_o_vel);
	std::ofstream img_out_pres(path_o_pres);

	// Output PPM P3 (ASCII) Header.
	img_out_dens << "P3" << std::endl;
	img_out_dens << w << " " << h << std::endl;
	img_out_dens << "255" << std::endl;

	img_out_vel << "P3" << std::endl;
	img_out_vel << w << " " << h << std::endl;
	img_out_vel << "255" << std::endl;

	img_out_pres << "P3" << std::endl;
	img_out_pres << w << " " << h << std::endl;
	img_out_pres << "255" << std::endl;

	// String Streams to Write To. 
	std::stringstream sout_dens, sout_vel, sout_pres; 

	// 0 Based Iter, Output edge Cells/Pixels. 
	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			// Density Output - 
			{
				float value = dens->getdata(i, j);
				int r = int(value * 255);
				int g = int(value * 255);
				int b = int(value * 255);

				// Clamp Min/Max. 
				if (r < 1) r = 1; if (r > 255) r = 255;
				if (g < 1) g = 1; if (g > 255) g = 255;
				if (b < 1) b = 1; if (b > 255) b = 255;

				sout_dens << r << " " << g << " " << b << std::endl;
			}

			// Velocity Output - 
			{
				// Do Abs to account for inversed sign on y... 
				float value_r = std::fabsf(vel->getdata_x(i, j));
				float value_g = std::fabsf(vel->getdata_y(i, j));
				float value_b = 0.0f; 

				int r = int(value_r * 255);
				int g = int(value_g * 255);
				int b = int(value_b * 255);

				if (r > 255) r = 255;
				if (g > 255) g = 255;
				if (b > 255) b = 255;

				sout_vel << r << " " << g << " " << b << std::endl;
			}

			// Pressure Output - 
			{
				float value = pres_grid->getdata(i, j);
				value *= 2.5f; // Brightness Mult. 

				int r = int(value * 255);
				int g = int(value * 255);
				int b = int(value * 255);

				if (r <= 1) r = 1; if (r >= 255) r = 255;
				if (g <= 1) g = 1; if (g >= 255) g = 255;
				if (b <= 1) b = 1; if (b >= 255) b = 255;

				sout_pres << r << " " << g << " " << b << std::endl;
			}
		}
	}

	// Set ofstream exceptions. 
	img_out_dens.exceptions(std::ios::badbit | std::ios::failbit);
	img_out_vel.exceptions(std::ios::badbit | std::ios::failbit);
	img_out_pres.exceptions(std::ios::badbit | std::ios::failbit);

	// Write SStream Data to ofstream output, within exeception handler. 
	// (Exceptions on writing only, ofstreams should be open already on construction, unhandled). 
	// Indivudal Exception Handling, Per SStream->OStream output. 
	try
	{
		img_out_dens << sout_dens.rdbuf();
	}
	catch (std::ofstream::failure e)
	{
		std::cerr << "ERROR::WriteTo_Img_All::Exception Writing Density Grid To Image:: " << e.what() << "\n";
		return 0;
	}

	try
	{
		img_out_vel << sout_vel.rdbuf();
	}
	catch (std::ofstream::failure e)
	{
		std::cerr << "ERROR::WriteTo_Img_All::Exception Writing Velocity Grid To Image:: " << e.what() << "\n";
		return 0; 
	}

	try
	{
		img_out_pres << sout_pres.rdbuf();
	}
	catch (std::ofstream::failure e)
	{
		std::cerr << "ERROR::WriteTo_Img_All::Exception Writing Pressure Grid To Image:: " << e.what() << "\n";
		return 0;
	}

	// Close ofstreams - 
	img_out_dens.close(); img_out_vel.close(); img_out_pres.close();

	return 1;
}

// TEST - Wrtie To Text/ASCII (Velocity Grid)
void fluidobj_2d::writeto_txt(int frame)
{
	// With Edge Cells - Edge Pixels.
	int w = x_s + e_s;
	int h = y_s + e_s;

	// Out Path String (paths).
	std::string paths = "build/dump/veltxt_";
	paths.append(std::to_string(frame) + ".txt");

	// Out Image Write OFilestream (img_out).
	std::ofstream txt_out(paths); // .open is called on construction (OFStream). 

	// 0 Based Iter, Output edge Cells/Pixels. 
	for (int j = 0; j < h; j++)
	{
		txt_out << "\n";
		for (int i = 0; i < w; i++)
		{
			vec2<float> v = vel->getdata(i, j);
			txt_out << std::fixed << v.x << ",";
			if (v.y < 0) txt_out << "**" << v.y << "| "; else txt_out << v.y << "| ";
			
		}
	}

	txt_out.close(); 
}



// If Not using Diffusion MFunc in Solver to Pass Current(n) Values to Prev_Step Grids
// Need to pass them manaully before caluclating n+1 values in current step grids. 

// WHY Not also loop and set edge cells??

// Write Cur_Dens Field to Prev_Dens Field MF - 

void fluidobj_2d::setcurtoprev_dens()
{
	// Loop Main Grid (Non Edge Cells) Set PrevDens to Cur. 
	
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 1; i < int(x_s * y_s); i++)
	{
		// Get Cur Density at Cell i
		vec2<int> idx_2d = dens->idx_1Dto2D(i);
		float dens_i = dens->getdata(idx_2d.x, idx_2d.y);

		// Set Prev Density at Cell i
		prev_dens->setdata(dens_i, idx_2d.x, idx_2d.y);
	}
	

	/*
	
	// Is Memcpy 'n the Cur Grid to Prev Faster than ^ ? 

	// Ofc Assume Grids are same size in mem. Else somethings gone wrong. Grid Vecs should have same size and reserve. 
	assert((dens->dataptr_getter())->size() == (prev_dens->dataptr_getter())->size());

	float *curdens_begin = (dens->dataptr_getter())->data(); // Pointer to Interal vector Data Array begin Cur. 
	float *prevdens_begin = (prev_dens->dataptr_getter())->data(); // Pointer to Interal vector Data Array begin Prev. 
	memcpy(prevdens_begin, curdens_begin, (std::size_t) ((dens->dataptr_getter())->size()) * sizeof(float)); // Grid Data Vector Size * sizof(float)

	*/ 
}


// Write Cur_Vel Field to Prev_Vel Field MF - 

void fluidobj_2d::setcurtoprev_vel()
{
	// Loop Main Grid (W Edge Cells) Set PrevVel to Cur. 
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 1; i < int((x_s * y_s) + e_s); i++)
	{
		// Get Cur Vel at Cell i
		vec2<int> idx_2d = vel->idx_1Dto2D(i);
		vec2<float> vel_i = vel->getdata(idx_2d.x, idx_2d.y);

		// Set Prev Density at Cell i
		prev_vel->setdata(vel_i, idx_2d.x, idx_2d.y);
	}
}


void fluidobj_2d::setcurtoprev(grid2_scalar *grid0, grid2_scalar *grid1)
{
	// Loop Main Grid (W Edge Cells) Set PrevVel to Cur. 
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 1; i < int((x_s * y_s) + e_s); i++)
	{
		// Get Scalar at current cell i 
		vec2<int> idx_2d = grid1->idx_1Dto2D(i);
		float scl_i = grid1->getdata(idx_2d.x, idx_2d.y);

		// Set Scalar at current cell i
		grid0->setdata(scl_i, idx_2d.x, idx_2d.y);
	}
}

void fluidobj_2d::setcurtoprev(grid2_vector *grid0, grid2_vector *grid1)
{
	// Loop Main Grid (W Edge Cells) Set PrevVel to Cur. 
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 1; i < int((x_s * y_s) + e_s); i++)
	{
		// Get Cur Vel at Cell i
		vec2<int> idx_2d = grid1->idx_1Dto2D(i);
		vec2<float> vel_i = grid1->getdata(idx_2d.x, idx_2d.y);

		// Set Prev Density at Cell i
		grid0->setdata(vel_i, idx_2d.x, idx_2d.y);
	}
}


// Colour Grid - 
int fluidobj_2d::RGB_imageLoad(const char* path)
{
	// Check if path exists 
	std::ifstream validate (path);
	if (!validate.good()) 
	{
		std::cout << "ERROR::IMAGE NOT FOUND \n \n";
		return 0;
	}
	else validate.close(); 

	// Set CUr Img Path to Member If Diffrent Img Passed - 
	if (!this->img_data.path)
	{
		img_data.path = const_cast<char *>(path);
	}
	else if (strcmp(path, this->img_data.path) != 0)
	{
		img_data.path = const_cast<char *>(path); // Update Path Mem with new Path str. 
	}

	// Flip so read TopDown Y (Like Grids iter) GL will Flip back to BotUp Y When read as RGB Texture. 
	stbi_set_flip_vertically_on_load(true);
	// Load Image Into Buffer, Store Dimensions and NChannels. 
	img_data.img_raw = stbi_load(path, &img_data.width, &img_data.height, &img_data.nch, 0);

	// Assume Image Size matches Current Grid Size (and thus is sqaure) (Minus edge cells, Only fill in 1-N). 

	int N_pix = img_data.width * img_data.height;
	int N_dim = x_s;
	// N_dim should == sqrt(N_pix) assuming grid is square. 
	assert(((int)std::sqrt(N_pix)) == N_dim);

	// Byte Stride Offset for Each Pixel (BYTE(Uchar) * NChannels)
	int stride = (int) sizeof(BYTE) * img_data.nch;


	// Write Image Channels Data to R,G,B Colour Grids. 
	// Doesnt Account for Edge Cell offsets (Shouldnt matter because only iter to N_pix)? 
	// Solve will ofc, but ideally initalze colour values NOT in edge cells (0,N+1) ij.

	int ch_idx = 0;
	for (int i = 0; i < N_pix; i++)
	{
		ch_idx += stride; // Offset Channel Index By Stride. 

		vec2<int> idx_2d = c_R->idx_1Dto2D(i); // 1D Pixel ID -> 2D Cell Index. 

		// Set RGB Grid. 
		int pix_r = (int)img_data.img_raw[ch_idx];
		int pix_g = (int)img_data.img_raw[ch_idx+1];
		int pix_b = (int)img_data.img_raw[ch_idx+2];

		c_R->setdata((float)pix_r, idx_2d.x, idx_2d.y);
		c_G->setdata((float)pix_g, idx_2d.x, idx_2d.y);
		c_B->setdata((float)pix_b, idx_2d.x, idx_2d.y);

	}

	// Zero Out Edge Cells (0,N+1) ij (Oppose to Edge Cells get neighbour vals).
	// Actually Do -1 from Edge Cells to get clear edge inital boder. 
	for (int i = 0; i <= N_dim; i++)
	{
		// X- Edge Boundary
		c_R->setdata(0.0f, 1, i);
		c_G->setdata(0.0f, 1, i);
		c_B->setdata(0.0f, 1, i);

		// X+ Edge Boundary
		c_R->setdata(0.0f, N_dim-1, i);
		c_G->setdata(0.0f, N_dim-1, i);
		c_B->setdata(0.0f, N_dim-1, i);

		// Y+ Edge Boundary
		c_R->setdata(0.0f, i, 1);
		c_G->setdata(0.0f, i, 1);
		c_B->setdata(0.0f, i, 1);

		// Y- Edge Boundary
		c_R->setdata(0.0f, i, N_dim-1);
		c_G->setdata(0.0f, i, N_dim-1);
		c_B->setdata(0.0f, i, N_dim-1);
	}

}
// End of RGB_Image Load. 

// TEST - Wrtie To Img Color 
int fluidobj_2d::writeto_img_col()
{
	// With Edge Cells - Edge Pixels.
	int w = x_s + e_s;
	int h = y_s + e_s;

	// Out Path String (paths).
	std::string paths = "build/img/coltest.ppm";
	//paths.append(std::to_string(frame) + ".ppm");
	if (verbose)
	{
		std::cout << "DEBUG::WriteToImage_Velocity_Output: " << paths << std::endl;
	}

	// Out Image Write OFilestream (img_out).
	std::ofstream img_out(paths); // .open is called on construction (OFStream). 
	img_out << "P3" << std::endl;
	img_out << w << " " << h << "\n";
	img_out << "255" << "\n";

	
	// 0 Based Iter, Output edge Cells/Pixels. 
	for (int j = 0; j < h; j++)
	{
		for (int i = 0; i < w; i++)
		{
			// RGB Values Derive from Velocity XY Compoennts (B is 0 in 2D).
			float value_r = c_R->getdata(i, j);
			float value_g = c_G->getdata(i, j);
			float value_b = c_B->getdata(i, j);

			int r = int(value_r);
			int g = int(value_g);
			int b = int(value_b);

			if (r <= 0) r = 1; if (r >= 255) r = 255;
			if (g <= 0) g = 1; if (g >= 255) g = 255;
			if (b <= 0) b = 1; if (b >= 255) b = 255;

			img_out << r << " " << g << " " << b << "\n";
		}
	}
	
	return 1; // If Sucessful Write - WIP -  Need to encap try catch and return 0 if expects thrown from OFStream etc.. 
}

// Write out Raw Img Data to PPM for Debug Sake. (Check for Parsing, vs Loading Errors)> 
int fluidobj_2d::ppm_imgtest()
{
	// Out Path String (paths).
	std::string paths = "build/img/coltest_b.ppm";
	std::ofstream img_out(paths);

	img_out << "P3 \n" << img_data.width << " " << img_data.height << "\n" << "255 \n";

	for (int i = 0; i < (int)((img_data.width * img_data.height) * img_data.nch); i++)
	{
		if (i % img_data.nch == 0) img_out << "\n";
		img_out << (int) img_data.img_raw[i] << " ";
	}

	img_out.close(); 

	return 1; 
}
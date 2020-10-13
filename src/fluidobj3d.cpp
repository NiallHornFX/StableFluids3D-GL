// Implementation of FluidObject_3D
#include "fluidobj3d.h"

// Std Headers
#include <fstream> 
#include <string>
#include <sstream>
#include <chrono>
#include <cassert>

#include <omp.h>

extern short verbose; 

// fluidobj_3d Implemenatation - 

fluidobj_3d::fluidobj_3d(int x_size, int y_size, int z_size, int edge_size) 
	: x_s(x_size), y_s(y_size), z_s(z_size), e_s(edge_size)
{
	t_s = (x_s + e_s) * (y_s + e_s) * (z_s + e_s);

	dens = new grid3_scalar<float>(x_size, y_size, z_size, edge_size);
	prev_dens = new grid3_scalar<float>(x_size, y_size, z_size, edge_size);
	col = new grid3_scalar<float>(x_size, y_size, z_size, edge_size);
	vel = new grid3_vector<vec3<float>>(x_size, y_size, z_size, edge_size);
	prev_vel = new grid3_vector<vec3<float>>(x_size, y_size, z_size, edge_size);
	curl = new grid3_vector<vec3<float>>(x_size, y_size, z_size, edge_size);
	vc = new grid3_vector<vec3<float>>(x_size, y_size, z_size, edge_size);
	preproj_vel = new grid3_vector<vec3<float>>(x_size, y_size, z_size, edge_size);
}


fluidobj_3d::~fluidobj_3d()
{
	delete dens; dens = nullptr; 
	delete prev_dens; prev_dens = nullptr;
	delete vel; vel = nullptr;
	delete prev_vel; prev_vel = nullptr; 
	delete col; col = nullptr; 

	if (curl) delete curl, curl = nullptr; 
	if (vc) delete vc, vc = nullptr; 
	if (preproj_vel) delete preproj_vel, preproj_vel = nullptr; 
}

// Pass In float Density at Given Grid Location (3D Index) 
void fluidobj_3d::add_density(float d, int i, int j, int k)
{
	dens->adddata(d, i, j, k);
}

// Pass In Vec3<float> Velocity at Given Grid Location (3D Index) 
void fluidobj_3d::add_velocity(const vec3<float> &v, int i, int j, int k)
{
	vel->adddata(v, i, j, k);
}

// Add Force To Velocity. Assume Constant Mass/Density 1.0f for now. PerFrame.
void fluidobj_3d::integrate_force(const vec3<float> &force, float dt, int i, int j, int k)
{
	vec3<float> temp = force;
	temp *= dt; 
	vel->adddata(temp, i, j, k);
}

void fluidobj_3d::print_info()
{
	// Just use sizeof (T) * GridData Size - 
	std::size_t sgrid_bytes = sizeof (float) * dens->griddatavector_getter().size(); 
	std::size_t vgrid_bytes = sizeof (vec3<float>) * vel->griddatavector_getter().size(); 

	// Correctly using 1024 B-KB-MB-GB Ranges. 
	std::size_t sgrid_KB = sgrid_bytes / (std::size_t)1024; // 1024 Btyes Per KB. 
	std::size_t sgrid_MB = sgrid_bytes / std::size_t(pow(1024, 2)); // 1048576 (1.048 Million) Btyes Per MB. 
	std::size_t sgrid_GB = sgrid_bytes / std::size_t(pow(1024, 3)); // 1073741824 (1.073 Billion) Btyes Per GB. 

	std::size_t vgrid_KB = vgrid_bytes / (std::size_t)1024; // 1024 Btyes Per KB. 
	std::size_t vgrid_MB = vgrid_bytes / std::size_t(pow(1024, 2)); // 1048576 (1.048 Million) Btyes Per MB. 
	std::size_t vgrid_GB = vgrid_bytes / std::size_t(pow(1024, 3)); // 1073741824 (1.073 Billion) Btyes Per GB. 

	std::cout << "\nDEBUG::Fluid Object 3D " << " INFO BEGIN - \n";
	std::cout << "DEBUG::Fluid Object 3D:: Scalar Grid Size = " << sgrid_bytes << " B \n";
	std::cout << "DEBUG::Fluid Object 3D:: Scalar Grid Size = " << sgrid_KB << " KB \n";
	std::cout << "DEBUG::Fluid Object 3D:: Scalar Grid Size = " << sgrid_MB << " MB \n";
	std::cout << "DEBUG::Fluid Object 3D:: Scalar Grid Size = " << sgrid_GB << " GB \n";
	std::cout << "DEBUG::Fluid Object 3D:: Scalar Cell Size = " << sizeof(float) << " B \n";
	std::cout << "============================================ \n";
	std::cout << "DEBUG::Fluid Object 3D:: vec3 Grid Size = " << vgrid_bytes << " B \n";
	std::cout << "DEBUG::Fluid Object 3D:: vec3 Grid Size = " << vgrid_KB << " KB \n";
	std::cout << "DEBUG::Fluid Object 3D:: vec3 Grid Size = " << vgrid_MB << " MB \n";
	std::cout << "DEBUG::Fluid Object 3D:: vec3 Grid Size = " << vgrid_GB << " GB \n";
	std::cout << "DEBUG::Fluid Object 3D:: vec3 Cell Size = " << sizeof(vec3<float>) << " B \n";
	std::cout << "DEBUG::Fluid Object 3D " << " INFO END. \n \n";

}



// Implicit SDF sphere For Sourcing into Density and Velocity Grids. 
void fluidobj_3d::implicit_sphere_source(float dd, const vec3<float> &vv, const vec3<float> &offset, float rad)
{
	float h = 1.0f / x_s; 
	float surf_isothresh = 0.01f; 

	#pragma omp parallel for
    for (std::size_t k = 1; k <= z_s; k++)
	{
		#pragma omp parallel for
        for (std::size_t j = 1; j <= y_s; j++)
		{
			#pragma omp parallel for
            for (std::size_t i = 1; i <= x_s; i++)
			{
				// Implicit Sphere/sphere Function: (X-u)^2 + (Y-v)^2 + (Z-w)^2 - r . >= 0 <= Thresh == Surface. > Thresh = Exterior. < 0 = Interior. Cells. 
				vec3<float> cell_gridSpace( float((i * h) - offset.x), float((j * h) - offset.y), float((k * h) - offset.z)); // Index 0-N to Grid Space 0-1. 
				float sphere_func = ((cell_gridSpace.x * cell_gridSpace.x) + (cell_gridSpace.y * cell_gridSpace.y) + (cell_gridSpace.z * cell_gridSpace.z)) - rad;

				// INSIDE Cells Sourcing - 
                if (sphere_func < surf_isothresh)
				{
					add_density(dd, i, j, k);
					add_velocity(vv, i, j, k);
				}

				/*
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
				*/

			}
		}
	}

}


void fluidobj_3d::setcurtoprev(grid3_scalar<float> *grid0, grid3_scalar<float> *grid1)
{
	// Swap Cur with Prev, And then set Prev to Cur. Oppose to Looping through directly. Much Faster. 
	grid0->griddatavector_getter().swap(grid1->griddatavector_getter());
	grid1->griddatavector_getter() = grid0->griddatavector_getter(); // Reassign swapped grid1 grid_data (now in grid0) back to grid1 post swap. 
	

	/* SLOW ! 
	// Loop Main Grid (W Edge Cells) Set PrevVel to Cur. 
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 0; i < t_s; i++)
	{
		// Get Scalar at current cell i 
		vec3<int> idx_3d = grid1->idx_1Dto3D(i);

		// Set Scalar at current cell i
		grid0->setdata((grid1->getdata(idx_3d.x, idx_3d.y, idx_3d.z)), idx_3d.x, idx_3d.y, idx_3d.z);
	}
	*/
}

void fluidobj_3d::setcurtoprev(grid3_vector<vec3<float>> *grid0, grid3_vector<vec3<float>> *grid1)
{
	grid0->griddatavector_getter().swap(grid1->griddatavector_getter());
	grid1->griddatavector_getter() = grid0->griddatavector_getter(); // Reassign back to grid1 post swap. 

	/* SLOW !
	// Loop Main Grid (W Edge Cells) Set PrevVel to Cur. 
	#pragma omp parallel for num_threads(omp_get_max_threads())
	for (int i = 0; i < t_s; i++)
	{
		vec3<int> idx_3d = grid1->idx_1Dto3D(i);

		// Set Prev Density at Cell i
		grid0->setdata((grid1->getdata(idx_3d.x, idx_3d.y, idx_3d.z)), idx_3d.x, idx_3d.y, idx_3d.z);
	}
	*/
}

// TEST - Wrtie To Text/ASCII (Velocity Grid)
void fluidobj_3d::writeto_txt(int frame)
{
	// With Edge Cells - Edge Pixels.
	std::size_t w = x_s + e_s, h = y_s + e_s, d = z_s + e_s;

	// Out Path String (paths).
	std::string paths = "build/dump/veltxt_";
	paths.append(std::to_string(frame) + ".txt");

	std::ofstream txt_out(paths); 

	if (txt_out.is_open())
	{
		for (std::size_t k = 0; k < d; k++)
		{
			for (std::size_t j = 0; j < h; j++)
			{
				txt_out << "\n";
				for (std::size_t i = 0; i < w; i++)
				{
					txt_out << "\n";
					vec3<float> v = vel->getdata(i, j, k);
					txt_out << std::fixed << v.x << "," << v.y <<  "," << v.z << "| ";
				}
			}
		}
	}
	txt_out.close();
}

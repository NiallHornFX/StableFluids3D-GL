# Stable Fluids 3D - (CPU) Fluid Solver 

![Demo GIF - Updated Feb 12th](gif_demo.gif) ![Demo GIF - Updated Feb 14th](gif_demo_bb.gif)

#### Note this project is been re-factored into a CPU + GPU Based Realtime Eulerian Fluid Simulation Project, i'm working on (Fall 2020). Keeping this for memories sake. 

This is a 3D Implementation of Jos Stams famous paper *Stable Fluids Paper [Stable Fluids](https://d2f99xq7vri1nk.cloudfront.net/legacy_app_files/pdf/ns.pdf "StableFluidsPaper") [1]* on the CPU. Implemented in C++ in an OOP Framework style, with OpenGL currently the only implemented Graphics API, For Ray Marching the resulting 3D Density Grid. 
 
The Advection Method implemented currently uses Semi-Lagrangian (Single Step, Forward Euler) or MacCormack Advection based of [Selle et al. 2007](http://physbam.stanford.edu/~fedkiw/papers/stanford2006-09.pdf "MacCormackPaper") [2], with ethier (Tri) Linear or Cosine Interpolation. 
With optional use of Vortex Confinement force, based on [Fedkiw et al. 2001](https://web.stanford.edu/class/cs237d/smoke.pdf "VortConfinePaper") [3] to boost vortices lost due to numerical dissipation. The Pressure Projection step uses a matrix-less implementation of the Gauss-Seidel (Single-Threaded) with Successive Over Relaxation or Jacobi (Multi-Threaded) iterative linear solvers to compute the pressure grid. While the Jacobi method can be safeley Multi-Threaded increasing performance, its simultaneous displacements
result in much slower convergence than the Gauss-Seidel method, espeically when using Successive Over Relaxation which increases convergence simmilar to the
Conjugate Gradient Method, because of its successive displacements been applied within each solve iteration, thus its favoured as the default pressure solver.

I Chose to implement my own templated vector and matrix classes oppose to using GLM. These are been further refined. 
The Application/Solve and Render Loop is currently embedded as part of the `fluidsolver_3` class itself, within the `fluidsolver_3::solve_step()` member function.

### Task List - April 2020

- [x] MacCormack Advection 
- [x] Vorticity Confinement
- [x] OpenMP Multithreading
- [ ] RayMarcher Refactor

## Dependencies (See "thirdpaty_License.md")
* GLFW - OpenGL Window and Context Creation.
* GLEW - OpenGL Extensions Loading/Wangiling.

## Building ##
Currently uses shipped Visual Studio 2017 (vc141) Solution. Libaries of GLFW and GLEW are statically Linked to the application as of now. Tested on Windows 7 and Windows 10. 

## Use Case - Setup Example ##
````C++
	// Create Render Context instance, to handle Graphics API Context and Window Creation -
	render_context_OGL render_c (win_size_xy, win_size_xy, short(GLMajor), short(GLMinor)); 

	// Create Fluid Object, Containing Fluid Grids and Data - 
	fluidobj_3d test_fluidobj (gridres_x, gridres_y, gridres_z, edge_size); 

	// Create FluidSolver Instance,  Pass FluidObj Pointer to It - 
	fluidsolver_3 test_fluidsolver (&test_fluidobj, timestep);

	// Pre Solve Paramaters Initial Values Set -
	test_fluidsolver.Parms.p_Do_Dens_Diff = false; 
	test_fluidsolver.Parms.p_Do_Dens_Disp = true;  
	test_fluidsolver.Parms.p_ProjectionType = test_fluidsolver.Parms.Project_GaussSeidel_SOR; 
	test_fluidsolver.Parms.p_SOR_alpha = 1.9f;
	test_fluidsolver.Parms.p_GS_Proj_iter = 5; 
	test_fluidsolver.Parms.p_AdvectionType = test_fluidsolver.Parms.Advect_SL_BackTrace_Euler;

	// Pass RenderContext Window -
	test_fluidsolver.set_window(render_c.get_window());

	// Call FluidSolver Solve_Step Member Function to begin Simulation and Rendering - 
	test_fluidsolver.solve_step(true, solve_steps);
````

## Implemented Papers ##

1. - J. Stam. Stable fluids. In Proc. of SIGGRAPH 99, pages 121–128, 1999.
2. - S. Andrew, F. Ronald, Kim, Byungmoon, Liu, Yingjie, and Rossignac, Jarek. An unconditionally stable maccormack method. J. Sci. Comput., 35(2-3), June 2008
3. - R. Fedkiw, J. Stam, AND JENSEN, H. 2001. Visual simulation ofsmoke. In Proc. of ACM SIGGRAPH 2001, 15–22.

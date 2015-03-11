//a full-multigrid program 
//weighted red-black gauss-seidel solver for u''=f(x,y) over the domain 0<x<1, 0<y<1, with dirichlet bondary conditions at x=0 and x=1 and at y=0 and at y=1

//the no. of terms in the analytical solution
#define NOFTER 100

//some system parameters

//system boundaries along x
#define SYS_LEFT_END 0.0
#define SYS_RIGHT_END 1.0

//system boundaries along y
#define SYS_BOTTOM_END 0.0
#define SYS_UP_END 1.0

//no. of intervals along x
#define FINEST_NO_OF_INTERVALS_X 64 //this must be a power of 2

//no. of intervals along y
#define FINEST_NO_OF_INTERVALS_Y 64 //this must be a power of 2

//parameters related to the FULL MULTIGRID SOLVER    
#define NO_OF_LEVELS 4
#define MIN_NUM_ITER 4
#define TOLER 1e-6
#define TOLER_RATE 0.5
#define OMEGA 1.0

//the parameters related to the creation of random intial condition
#define INP_SEED 654321
#define NOISE_AMPLITUDE 1.0


//the choice for the right hand side (the other options are "ZERO_RHS" which would make it a Laplace equation)
#define POINT_SOURCE

//===========================================================================================//

//the location of the point source
#define LOCATE_POINT_SOURCE_X 0.67 //this is given as a fraction of the length of the entire domain along the x  
#define LOCATE_POINT_SOURCE_Y 0.67 //this is given as a fraction of the length of the entire domain along the y  	
#define POINT_SOURCE_MAGNITUDE 1000.0//the magnitude of the point source
//==========================================================================================//


//the choice for the boundary conditons (all are Dirichlet BCs)
//I am free to use non-homogeneous BCs as well
#define LEFT 0.0 //it is a no-flux BC along x on the left hand boundary 
#define RIGHT 0.0 //it is a dirichlet along x on the right hand boundary
#define UP 0.0 //it is a dirichlet along y on the upper hand boundary
#define BELOW 0.0//it is a no-flux along y on the lower boundary

//the choice for the inital guesses (the other option I have is RANDOM)
#define ZERO

//some header files we are going to use 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//============================================================================//
//some variables needed by the system

//declaring some structure types 
struct type_long_dim2 {
	long x;
	
	long y;
};

struct type_double_dim2 {
	double x;
	
	double y;
};	

struct type_long_dim2 nodes[NO_OF_LEVELS];//this will contain the information about the no. of nodes in the system  
struct type_double_dim2 d_sp[NO_OF_LEVELS];//this will contain the indormation on the grid spacing in the system 

struct type_long_dim2 nodes_sum; //the total no. of nodes in the system 

double g_s_prefac[NO_OF_LEVELS],sq_dx_times_sq_dy[NO_OF_LEVELS],sq_dx[NO_OF_LEVELS],sq_dy[NO_OF_LEVELS],corner_prefac[NO_OF_LEVELS];//some quantities that are required by the gausss-seidel solver

long SEED=INP_SEED;//this is required for generating the random nos., as the random no. generator doesn't recognize the pointers to macros 

long start[NO_OF_LEVELS],end[NO_OF_LEVELS];//they are going to hold the starting and finishing array_index for different grid sizes  
//'0' corresponds to the finest level and so on 

long total_no_of_nodes;//this variable stores the array lengths which are needed 
//=====================================================================//	


//====================================================================//
//some source function files
#include"ran2.c"
#include"initialise.c"
#include"set_dirichlet_boundary.c"
#include"set_flux_boundary.c"
#include"set_rhs.c"

#include"even_red_black_solver.c"
#include"odd_red_black_solver.c"
#include"relaxation.c"

#include"coarse_to_fine.c"

#include"compute_residual.c"
#include"fine_to_coarse.c"

#include"write_iterations.c"	
#include"multigrid_solver.c"

#include"write_profile.c"
#include"analytical_solution.c"
//===================================================================//


main()
{
	
	//==============================================================================//
	int i;//this will iterate through the no. of different meshes we have employed
	//==============================================================================//

	//==============================================================================//
	//computing the exponent (N) in 2^N=FINEST_NO_OF_INTERVALS
	struct type_long_dim2 finest_mesh_index;
	
	finest_mesh_index.x=(long)((log((double)(FINEST_NO_OF_INTERVALS_X))/log(2.0)));

	finest_mesh_index.y=(long)((log((double)(FINEST_NO_OF_INTERVALS_Y))/log(2.0)));

	printf("The finest mesh index along x is=%ld\t along y is=%ld\n",finest_mesh_index.x,finest_mesh_index.y);
	 	  
	if((int)(pow((double)2,(double)finest_mesh_index.x))!=FINEST_NO_OF_INTERVALS_X)
	{
		printf("The no. of intervals along x in the finest grid is not a power of 2\n");
		exit(1);
	}
	if((int)(pow((double)2,(double)finest_mesh_index.y))!=FINEST_NO_OF_INTERVALS_Y)
	{
		printf("The no. of intervals along y in the finest grid is not a power of 2\n");
		exit(1);
	}

	//==============================================================================//

	//==============================================================================//
	//finding total storage required considering all the grids andother important quantities
	//also computing the constants which are required in the computation of the derivatives
	
	nodes_sum.x=0;
	nodes_sum.y=0;

	struct type_long_dim2 mesh_index;

	struct type_long_dim2 mesh_intervals;	

	total_no_of_nodes=0;
 
	for(i=0;i<=NO_OF_LEVELS-1;i++)
	{
		//the mesh_index for the current level
		//for x
		mesh_index.x=finest_mesh_index.x-i;	
		//for y
		mesh_index.y=finest_mesh_index.y-i;	

		//the no. of mesh_intervals	
		//for x	
		mesh_intervals.x=(long)(pow((double)2,(double)mesh_index.x));
		//for y
		mesh_intervals.y=(long)(pow((double)2,(double)mesh_index.y));

		//the no. of nodes (including boundary points)
		//for x
		nodes[i].x=mesh_intervals.x+1;
		//for y
		nodes[i].y=mesh_intervals.y+1;

		//updating the sum of the nodes
		//for x
		nodes_sum.x+=nodes[i].x;
		//for y
		nodes_sum.y+=nodes[i].y;

		//computing the dx for the current level 
		//for x
		d_sp[i].x=(SYS_RIGHT_END-SYS_LEFT_END)/mesh_intervals.x;
		//for y
		d_sp[i].y=(SYS_UP_END-SYS_BOTTOM_END)/mesh_intervals.y;


		//using that arrays are going to be written in a flattened out fashion, i.e., first the entire 2D grid is  
		//computing the starting and ending array indices for every level such that each 2D grid is flattened out and stored one after the another
		if(i==0)//the finest level
		{
			start[i]=0;
			end[i]=nodes[i].x*nodes[i].y-1;

		}

		else //other coarser levels
		{
			
			start[i]=end[i-1]+1;//it starts just to the right of where the previous grid
			end[i]=start[i]+(nodes[i].x*nodes[i].y-1);	

			
		}		
			

		printf("The no. of nodes along x at the %dth level with dx=%lf is=%ld\n",i,d_sp[i].x,nodes[i].x);
		printf("The no. of nodes along y at the %dth level with dy=%lf is=%ld\n",i,d_sp[i].y,nodes[i].y);

		printf("it starts at=%ld and finishes at=%ld\n",start[i],end[i]); 

		

		//computing the total no of nodes		
		total_no_of_nodes+=nodes[i].x*nodes[i].y;



		//computing the prefactors
		g_s_prefac[i]=1.0/(2.0*((d_sp[i].x*d_sp[i].x)+(d_sp[i].y*d_sp[i].y)));

		corner_prefac[i]=1.0/((d_sp[i].x*d_sp[i].x)+(d_sp[i].y*d_sp[i].y));
	
		sq_dx_times_sq_dy[i]=d_sp[i].x*d_sp[i].x*d_sp[i].y*d_sp[i].y;

		sq_dx[i]=d_sp[i].x*d_sp[i].x;

		sq_dy[i]=d_sp[i].y*d_sp[i].y;	

		printf("At the level=%d\n",i);

		printf("The different prefactors are g_s_prefac=%lf and sq_dx_times_sq_dy=%lf\n",g_s_prefac[i],sq_dx_times_sq_dy[i]);




	}

	printf("the total no. of nodes along x=%ld\n",nodes_sum.x);

	printf("the total no. of nodes along y=%ld\n",nodes_sum.y);

	printf("the total no. of nodes in the system=%ld\n",total_no_of_nodes);


	getchar();

	//==============================================================================//
	
	//============================================================================//	
	//the arrays we are going to be needing
		
	double *y,*ana_y; //the dependent variable arrays
	double *f,*g; //the RHS array
	

	
	//allocating memory to the arrays
	
	if((y=(double *)malloc(total_no_of_nodes*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the y array\n");
		exit(1);
	}

	if((ana_y=(double *)malloc(nodes[0].x*nodes[0].y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the ana_y array\n");
		exit(1);
	}
		 	
	if((f=(double *)malloc(total_no_of_nodes*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the f array\n");
		exit(1);
	}

	if((g=(double *)malloc(nodes[0].x*nodes[0].y*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the g array\n");
		exit(1);
	}
	//============================================================================//	

	//============================================================================//
	//initial guess (see "initialise.c") for the unknown points
	init(y);
	//============================================================================//

	//============================================================================//
	//setting up the RHS (see "set_rhs.c")
	rhs(f,g);	
	//============================================================================//

	//============================================================================//
	//setting up the dirichlet boundary conditions (see "set_dirichlet_boundary.c")
	dirich_boundary(y);
	//============================================================================//	

	//============================================================================//
	//setting up the flux boundary conditions (see "set_flux_boundary.c")
	flux_boundary(f);
	//============================================================================//
	
	//===========================================================================//
	//solving the equations on multiple grids (see "multigrid_solver.c")
	mult_solve(y,f);
	//===========================================================================//	
	
	//===========================================================================//
	//writing the final files
	//see "write_profile.c"
	write(y);
	//===========================================================================//

	
	//computing the analytical solutions (see "analytical_solution.c")
	ana(ana_y,g,y);	

	free(ana_y);
	free(y);
	free(g);
	free(f);	
}
					

	
		
	

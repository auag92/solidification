//a full-multigrid program 
//weighted red-black gauss-seidel solver for y''=f(x) over the domain 0<x<1 with dirichlet bondary condition at x=0 and x=1

//some system parameters
#define SYS_LEFT_END 1.0
#define SYS_RIGHT_END 2.0
#define FINEST_NO_OF_INTERVALS 64 //this must be a power of 2   
#define NO_OF_LEVELS 4
#define MIN_NUM_ITER 4
#define TOLER 1e-6
#define TOLER_RATE 0.5
#define OMEGA 1.0
#define INP_SEED 654321
#define NOISE_AMPLITUDE 1.0


//the choice for the right hand side
#define X_RHS//the  choices can be "ZERO_RHS", 'X_RHS'

//the choice for the boundary conditons 
#define LEFT 0.0 //it denotes the value of the dependent variable at the left boundary  
#define RIGHT 0.0//it  denotes the value of the dependent variable at the right boundary  

//the choice for the inital guesses
#define ZERO //the options can be "RANDOM" OR "ZERO"

//the choice for the function 'a(x)'
#define A_X_BY_100 //this will reduce it to the usual formulation (the options are: 'A_ONE' and 'A_SIX_TIMES_X_SQ') 	


//some header files we are going to use 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//some variables needed by the system
long nodes_x[NO_OF_LEVELS],mesh_intervals[NO_OF_LEVELS];//the no. of nodes in the system
double dx[NO_OF_LEVELS];//the different grid spacings
long nodes_x_sum,mesh_intervals_x_sum; //the total no. of nodes in the system; and the total no. of mesh intervals in the system
long start[NO_OF_LEVELS],end[NO_OF_LEVELS];//they are going to hold the starting and finishing array_index for different grid sizes  
long start_mesh_inter[NO_OF_LEVELS],end_mesh_inter[NO_OF_LEVELS];//starting and ending of the mesh intervals 
long SEED=INP_SEED;
//'0' corresponds to the finest level and so on 


//some source function files
#include"ran2.c"
#include"initialise.c"
#include"set_dirichlet_boundary.c"
#include"creating_a.c"
#include"set_rhs.c"

#include"red_black_solver.c"
#include"relaxation.c"

#include"coarse_to_fine.c"
#include"fine_to_coarse.c"

#include"write_iterations.c"	
#include"multigrid_solver.c"

#include"write_profile.c"


main()
{
	
	//==============================================================================//
	int i;
	//==============================================================================//

	//==============================================================================//
	//computing the exponent (N) in 2^N=FINEST_NO_OF_INTERVALS
	int finest_mesh_index;
	
	finest_mesh_index=(int)((log((double)(FINEST_NO_OF_INTERVALS))/log(2.0)));

	printf("The finest mesh index is=%d\n",finest_mesh_index);
	 	  
	if((int)(pow((double)2,(double)finest_mesh_index))!=FINEST_NO_OF_INTERVALS)
	{
		printf("The no. of intervals in the finest grid is not a power of 2\n");
		exit(1);
	}
	//==============================================================================//

	//==============================================================================//
	//finding total storage required considering all the grids andother important quantities
	
	nodes_x_sum=0;

	mesh_intervals_x_sum=0;

	int mesh_index;

	
 
	for(i=0;i<=NO_OF_LEVELS-1;i++)
	{
		//the mesh_index for the next level
		mesh_index=finest_mesh_index-i;	

		//the no. of mesh_intervals in the x-direction		
		mesh_intervals[i]=(long)(pow((double)2,(double)mesh_index));

		//the no. of nodes in the x-direction (including boundary points)
		nodes_x[i]=mesh_intervals[i]+1;

		//updating the sum of the nodes_x
		nodes_x_sum+=nodes_x[i];

		//updating the sum of the mesh_intervals
		mesh_intervals_x_sum+=mesh_intervals[i];

		//computing the dx for the current level 
		dx[i]=(SYS_RIGHT_END-SYS_LEFT_END)/mesh_intervals[i];

		//computing the starting and ending array indices for every level 
		if(i==0)//the finest level
		{
			start[i]=0;
			end[i]=nodes_x[i]-1;

			start_mesh_inter[i]=0;
			end_mesh_inter[i]=mesh_intervals[i]-1;
		}

		else //other coarser levels
		{
			start[i]=end[i-1]+1;//it starts just to the right of where the previous grid
			end[i]=start[i]+(nodes_x[i]-1);	


			start_mesh_inter[i]=end_mesh_inter[i-1]+1;//it starts just to the right of where the previous grid
			end_mesh_inter[i]=start_mesh_inter[i]+(mesh_intervals[i]-1);	
		}		
			

		printf("The no. of nodes at the %dth level with dx=%lf is=%ld\n",i,dx[i],nodes_x[i]);

		printf("it starts at=%ld and finishes at=%ld\n",start[i],end[i]); 

		
	}

	printf("the total no. of nodes_x=%ld\n",nodes_x_sum);

	printf("the total no. of mesh_intervals=%ld\n",mesh_intervals_x_sum);


	getchar();

	//==============================================================================//
	

	
	//============================================================================//	
	//the arrays we are going to be needing
		
	double *y; //the dependent variable arrays
	double *f; //the RHS array
	double *a;//the variable coefficient function
	

	
	//allocating memory to the arrays
	
	if((y=(double *)malloc(nodes_x_sum*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the y array\n");
		exit(1);
	}
		 	
	if((f=(double *)malloc(nodes_x_sum*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the f array\n");
		exit(1);
	}

	if((a=(double *)malloc(mesh_intervals_x_sum*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the a array\n");
		exit(1);
	}

	
	//============================================================================//	

	//============================================================================//
	//initial guess (see "initialise.c") for the unknown points
	init(y);
	//============================================================================//

	//============================================================================//
	//setting up the RHS (see "set_rhs.c")
	rhs(f);	
	//============================================================================//

	//============================================================================//
	//creating the 'a' array (see "creating_a.c")
	cr_a(a);
	//============================================================================//

	//============================================================================//
	//setting up the boundary conditions for the finest grid only(see "set_dirichlet_boundary.c")
	//the coarser grids will have zero boundary conditions
	dirichlet_boundary(y);
	//============================================================================//	
	 		
	
	//===========================================================================//
	//solving the equations on multiple grids (see "multigrid_solver.c")
	mult_solve(y,f,a);
	//===========================================================================//	
	
	//===========================================================================//
	//writing the final files
	//see "write_profile.c"
	write(y);
	//===========================================================================//
}
					

	
		
	

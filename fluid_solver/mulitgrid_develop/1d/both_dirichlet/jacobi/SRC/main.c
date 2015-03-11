//weighted jacobi solver for y''=f(x) over the domain 0<x<1 with dirichlet bondary condition at x=0 and x=1

//some system parameters
#define sys_size_x 1.0
#define dx 0.01
#define omega 0.6667
#define SEED 654321
#define noise_amplitude 1.0
#define toler 1e-6
#define output_interval 5000 

//the choice for the right hand side
#define X_RHS

//the choice for the boundary conditons
#define LEFT 0.0
#define RIGHT 0.0

//the choice for the inital guesses
#define ZERO	

//some header files we are going to use 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//some variables needed by the system
int nodes_x;//the no. of nodes in the system

//some source function files
#include"ran2.c"
#include"initialise.c"
#include"set_boundary.c"
#include"set_rhs.c"
#include"weighted_jacobi_solver.c"
#include"error_compute.c"
#include"copying_solutions.c"
#include"write_profile.c"



main()
{
	
	//==============================================================================//
	int iter=0;

	double error=0.0;
	//==============================================================================//

	//============================================================================//
	nodes_x=(int)(sys_size_x/dx)+1;

	printf("the no. of nodes in the system=%d\n",nodes_x);
	//============================================================================//


	//============================================================================//	
	double *y_old,*y; //the dependent variable arrays
	double *f; //the RHS array

	//allocating memory to the arrays
	if((y_old=(double *)malloc(nodes_x*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the y_old array\n");
		exit(1);
	}
	
	if((y=(double *)malloc(nodes_x*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the y array\n");
		exit(1);
	}
		 	
	if((f=(double *)malloc(nodes_x*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the f array\n");
		exit(1);
	}
	
	//============================================================================//	

	//============================================================================//
	//initial guess (see "initialise.c") for the unknown points
	init(y_old);
	//============================================================================//

	//============================================================================//
	//setting up the boundary conditions (see "set_boundary.c")
	boundary(y_old);
	//============================================================================//	

	//============================================================================//
	//setting up the RHS (see "set_rhs.c")
	rhs(f);	
	//============================================================================//

	//strating the iterations
	while(iter<2 || error>toler)	
	{

		//============================================================================//
		//solving the equation (see "weighted_jacobi_solver.c")
		jacobi(y_old,y,f);
		//===========================================================================//	

		//===========================================================================//
		//computing the errors (see "error_compute.c")
		error=comp_error(y_old,y);	
		//===========================================================================//
	
		//===========================================================================//
		iter++;
		printf("iter=%d\terror=%lf\n",iter,error);
		//===========================================================================//

		//checking whether files have to be written
		if(iter%output_interval==0)
		{
			printf("Files are going to be written\n");
		
			//see "write_profile.c"
			write(y,iter);
		}

		//transferring y to y_old (see "copying_solutions.c")
		copy_soln(y,y_old);		

	}//closing the iteration loop

	//again_writing the final files
	//see "write_profile.c"
	write(y,iter);	
}
	
		
				

						

			

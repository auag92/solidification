//weighted red-black gauss-seidel solver for y''=f(x) over the domain 0<x<1 with dirichlet bondary condition at x=0 and x=1

//some system parameters
#define sys_size_x 1.0
#define dx 0.015625
#define omega 1.0
#define SEED 654321
#define noise_amplitude 1.0
#define toler 1e-6
#define output_interval 1000

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
#include"red_black_solver.c"
#include"write_profile.c"



main()
{
	
	//==============================================================================//
	int iter=0;

	double error=0.0;

	int start;
	//==============================================================================//

	//============================================================================//
	nodes_x=(int)(sys_size_x/dx)+1;

	printf("the no. of nodes in the system=%d\n",nodes_x);
	//============================================================================//


	//============================================================================//	
	double *y; //the dependent variable arrays
	double *f; //the RHS array

	//allocating memory to the arrays
	
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
	init(y);
	//============================================================================//

	//============================================================================//
	//setting up the boundary conditions (see "set_boundary.c")
	boundary(y);
	//============================================================================//	

	//============================================================================//
	//setting up the RHS (see "set_rhs.c")
	rhs(f);	
	//============================================================================//

	//strating the iterations
	while(iter<2 || error>toler)	
	{

		error=0.0;

		//============================================================================//
		//solving the equation for the even update only(see "red_black_solver.c")
		start=2;
		rb_solve(y,f,start,&error);
		//===========================================================================//	

		
		//============================================================================//
		//solving the equation for the odd update only(see "red_black_solver.c")
		start=1;
		rb_solve(y,f,start,&error);
		//===========================================================================//	
	
		//===========================================================================//
		error=sqrt(error*dx);
		//===========================================================================//
	
		//===========================================================================//
		iter++;
		printf("iter=%d\terror=%lf\n",iter,error);

	//	getchar();
		//===========================================================================//

		//checking whether files have to be written
		if(iter%output_interval==0)
		{
			printf("Files are going to be written\n");
		
			//see "write_profile.c"
			write(y,iter);
		}

		

	}//closing the iteration loop

	//again_writing the final files
	//see "write_profile.c"
	write(y,iter);	
}
	
		
				

						

			
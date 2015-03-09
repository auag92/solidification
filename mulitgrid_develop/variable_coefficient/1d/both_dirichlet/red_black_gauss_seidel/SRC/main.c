//weighted red-black gauss-seidel solver for y''=f(x) over the domain 0<x<1 with Dirichlet at x=0 and at x=1
//but the exact poisson equation we are going to solve here is one with variable coefficients:
// i.e., \nabla \cdot (a(x)\nabla y)=f(x); where 'a(x)' is a kown function of x


//some system parameters

//=====================================//
//for setting up the domains
#define SYS_LEFT_END 1.0
#define SYS_RIGHT_END 2.0
#define NO_OF_INTERVALS 64
#define OMEGA 1.0
#define INP_SEED 654321
#define NOISE_AMPLITUDE 1.0
#define TOLER 1e-10
#define OUTPUT_INTERVAL 1000
//====================================//

//the choice for the right hand side
#define X_RHS//the  choices can be "ZERO_RHS", 'X_RHS'

//the choice for the boundary conditons 
#define LEFT 0.0 //it denotes the value of the dependent variable at the left boundary  
#define RIGHT 0.0//it  denotes the value of the dependent variable at the right boundary  

//the choice for the inital guesses
#define ZERO //the options can be "RANDOM" OR "ZERO"

//the choice for the function 'a(x)'
#define A_SIX_TIMES_X_SQ //this will reduce it to the usual formulation (the options are: 'A_ONE' and 'A_SIX_TIMES_X_SQ') 	

//some header files we are going to use 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//some variables needed by the system
int nodes_x;//the no. of nodes in the system
double dx;

//some source function files
#include"ran2.c"
#include"initialise.c"
#include"creating_a.c" 
#include"set_dirichlet_boundary.c"
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

	//==============================================================================//
	//computing the exponent (N) in 2^N=NO_OF_INTERVALS
	int mesh_index;
	
	mesh_index=(int)((log((double)(NO_OF_INTERVALS))/log(2.0)));

	printf("The mesh index is=%d\n",mesh_index);
	 	  
	if((int)(pow((double)2,(double)mesh_index))!=NO_OF_INTERVALS)
	{
		printf("The no. of intervals in the grid is not a power of 2\n");
		exit(1);
	}
	//==============================================================================//

	//==============================================================================//
	//computing the dx for the current level 
	dx=(SYS_RIGHT_END-SYS_LEFT_END)/NO_OF_INTERVALS;

	nodes_x=NO_OF_INTERVALS+1;

	printf("the no. of nodes in the system=%d\n",nodes_x);
	//============================================================================//


	//============================================================================//	
	double *y; //the dependent variable arrays
	double *f; //the RHS array
	double *a;//the varaible coefficient function

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

	if((a=(double *)malloc(NO_OF_INTERVALS*sizeof(double)))==NULL) //I am going to create an 'a' array which contains values at the centre of the intervals
	{
		printf("Space creation failed for the 'a' array\n");
		exit(1);
	}
	
	//============================================================================//	

	//============================================================================//
	//initial guess (see "initialise.c") for the unknown points
	init(y);
	//============================================================================//

	//============================================================================//
	//creating the 'a' array (see "creating_a.c")
	cr_a(a);
	//============================================================================//
	
	//============================================================================//
	//setting up the RHS (see "set_rhs.c")
	rhs(f);	
	//============================================================================//

	//============================================================================//
	//setting up the boundary conditions (see "set_dirichlet_boundary.c")
	dirichlet_boundary(y);
	//============================================================================//

	
	//strating the iterations
	while(iter<2 || error>TOLER)	
	{

		error=0.0;

		

		//============================================================================//
		//solving the equation for the even update only(see "red_black_solver.c")
		start=2;
		rb_solve(y,f,a,start,&error);
		//===========================================================================//	

		
		//============================================================================//
		//solving the equation for the odd update only(see "red_black_solver.c")
		start=1;
		rb_solve(y,f,a,start,&error);
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
		if(iter%OUTPUT_INTERVAL==0)
		{
			printf("Files are going to be written\n");
		
			//see "write_profile.c"
			write(y,iter);
		}

		

	}//closing the iteration loop

	//again_writing the final files
	//see "write_profile.c"
	write(y,iter);

	free(y);
	free(f);	
	free(a);
}
	
		
				

						

			
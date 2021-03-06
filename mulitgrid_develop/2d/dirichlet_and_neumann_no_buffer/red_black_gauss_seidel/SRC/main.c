//weighted red-black gauss-seidel solver for u''=f(x,y) over the domain 0<x<1,0<y<1 with no-flux on boundaries adjoining the origin and dirichlet and set to '0' on the opposite boundaries
//in the code we should keep a provision for specifying the non-zero values 

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
#define NO_OF_INTERVALS_X 64 //this must be a power of 2 

//no. of intervals along y
#define NO_OF_INTERVALS_Y 64 //this must be a power of 2

//some other parameters
#define OMEGA 1.0 //this is for the weighted gauss-seidel

//this is for the random intial guess 
#define INP_SEED 654321
#define NOISE_AMPLITUDE 1.0

//solver specifications
#define TOLER 1e-6
#define OUTPUT_INTERVAL 1000

//the choice for the right hand side (the other options are "ZERO_RHS" which would make it a Laplace equation and "FUNCTION_RHS")
#define POINT_SOURCE

//==============================================================================================//
//they are related to the point source

//the location of the point source
#define LOCATE_POINT_SOURCE_X 0.67 //this is given as a fraction of the length of the entire domain along the x  
#define LOCATE_POINT_SOURCE_Y 0.67 //this is given as a fraction of the length of the entire domain along the y  	
#define POINT_SOURCE_MAGNITUDE 1000.0//the magnitude of the point source
//==============================================================================================//


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

//=====================================//
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

struct type_long_dim2 nodes;//this will contain the information about the no. of nodes in the system  
struct type_double_dim2 d_sp;//this will contain the indormation on the grid spacing in the system 

long tot_array_length;

double g_s_prefac,sq_dx_times_sq_dy,sq_dx,sq_dy,corner_prefac;

long SEED=INP_SEED;

//=====================================//	


//=====================================//
//some source function files
#include"ran2.c"
#include"initialise.c"
#include"set_dirichlet_boundary.c"
#include"set_flux_boundary.c"
#include"set_rhs.c"
#include"even_red_black_solver.c"
#include"odd_red_black_solver.c"
#include"write_profile.c"
#include"analytical_solution.c"
//====================================//



main()
{
	
	//==============================================================================//
	int iter=0;

	double error=0.0;

	//==============================================================================//
	//computing the no. of nodes in the system 
	//the no. of nodes in the x-direction (including boundary points)
	nodes.x=NO_OF_INTERVALS_X+1;

	//the no. of nodes in the y-direction (including boundary points)
	nodes.y=NO_OF_INTERVALS_Y+1;

	printf("the no. of nodes in the system along x=%ld\n",nodes.x);
	
	printf("the no. of nodes in the system along y=%ld\n",nodes.y);

	//computing the length of the array
	tot_array_length=nodes.x*nodes.y;
	
	printf("The total array length=%ld\n",tot_array_length);

	//============================================================================//

	//============================================================================//
	//computing the d_sp for the grid we have employed 

	//along x	
	d_sp.x=(SYS_RIGHT_END-SYS_LEFT_END)/NO_OF_INTERVALS_X;
	//along y
	d_sp.y=(SYS_UP_END-SYS_BOTTOM_END)/NO_OF_INTERVALS_Y;


	printf("The grid spacing along x=%lf\t and along y=%lf\n",d_sp.x,d_sp.y);
	//=============================================================================// 
		


	//============================================================================//	
	double *y; //the dependent variable arrays
	double *f,*g; //the RHS array

	double *ana_y;

	//allocating memory to the arrays
	
	if((y=(double *)malloc(tot_array_length*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the y array\n");
		exit(1);
	}

	if((ana_y=(double *)malloc(tot_array_length*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the ana_y array\n");
		exit(1);
	}
		 	
	if((f=(double *)malloc(tot_array_length*sizeof(double)))==NULL)
	{
		printf("Space creation failed for the f array\n");
		exit(1);
	}

	if((g=(double *)malloc(tot_array_length*sizeof(double)))==NULL)
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
	//setting up the RHS (see "set_rhs.c")
	//I am not going to modify my rhs to take into account the flux BC's in this routine 
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
	

	//============================================================================//
	//computing the prefactors
	g_s_prefac=1.0/(2.0*((d_sp.x*d_sp.x)+(d_sp.y*d_sp.y)));
	
	sq_dx_times_sq_dy=d_sp.x*d_sp.x*d_sp.y*d_sp.y;

	sq_dx=d_sp.x*d_sp.x;

	sq_dy=d_sp.y*d_sp.y;	

	corner_prefac=1.0/((d_sp.x*d_sp.x)+(d_sp.y*d_sp.y));	

	printf("The different prefactors are g_s_prefac=%lf and sq_dx_times_sq_dy=%lf\tcorner_prefac=%lf\n",g_s_prefac,sq_dx_times_sq_dy,corner_prefac);

	//============================================================================//

	//strating the iterations
	while(iter<2 || error>TOLER)	
	{

		error=0.0;

		//============================================================================//
		//solving the equation for the even update only(see "even_red_black_solver.c")
		
		even_rb_solve(y,f,&error);
		//===========================================================================//	

		
		//============================================================================//
		//solving the equation for the odd update only(see "odd_red_black_solver.c")
		
		odd_rb_solve(y,f,&error);
		//===========================================================================//	
	
		//===========================================================================//
		error=sqrt(error*d_sp.x*d_sp.y);
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

	//computing the analytical solutions (see "analytical_solution.c")
	ana(ana_y,g,y);	

	free(ana_y);
	free(y);
	free(g);
	free(f);	
	
}
	
		
				

						

			

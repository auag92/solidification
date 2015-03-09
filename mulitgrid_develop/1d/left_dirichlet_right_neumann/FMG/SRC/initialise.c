void init(double *y)
{
	long i;

	//in principle I could have initialised anything for the finest grid and kept the other grids at zero

	//initialising the array in the finest grid 

	//we already know that the level zero is the finest grid

	for(i=start[0]+1;i<=end[0];i++) //I have initialised everything excluding Dirichlet boundary points for the finest grid; the boundary conditions will bw set correctly when I call the next function from main
	//the boundary conditions are tackled by modifiying the rhs in this case 
	//as we are going to solve for a Neumann boundary condition at the right end and a Dirichlet at the left end 
	{
		#ifdef ZERO //I am going to fill it with zeroes 
 		y[i]=0.0;
		#endif

		#ifdef RANDOM //I am going to fill it with random nos. which can be useful when I solve y''=0 
		y[i]=(2.0*ran2(&SEED)-1.0)*NOISE_AMPLITUDE;
		#endif
	}

	//initialising the other portions of the solution array; they must be initialised to zero regardless of the exact problem I am trying to solve

	//they are just solving the residual equation so the initial guess should be zero

	//the boundary conditons for the coarser grids are going to be implicitly taken care of by the setting the correct rhs

	for(i=start[1];i<=end[NO_OF_LEVELS-1];i++) //we are starting at the grid which is second finest and we go till the end of the array
	{
		y[i]=0.0;
	} 
	 
	  	

}

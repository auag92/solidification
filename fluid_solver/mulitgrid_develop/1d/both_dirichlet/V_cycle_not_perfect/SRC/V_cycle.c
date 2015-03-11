void solve_V_cycle(double *y,double *f)
{
	int count;

	//============================================================================================//
	//going from the finest to the coarsest 		
	for(count=0;count<no_of_levels-1;count++)//the loop which will iterate over the no. of levels except for the coarsest grid 
	{
		//relaxing on the grid using red-black gauss-seidel (see "relaxation.c")
		relax(y,f,count);

		//computing the residual to the above calculation and restricting this to the next coarser grid (using the "full weighting" scheme)
		//see "compute_residual.c"
		comp_res(y,f,count);
	}
	//============================================================================================//

	printf("here\n");
	getchar();

	//============================================================================================//
	//going from the coarsest to the finest
	for(count=no_of_levels-1;count>0;count--)	 //the loop which will iterate over the no. of levels except for the finest grid
	{
		//solving the equation
		//relaxing on the grid using red-black gauss-seidel (see "relaxation.c")
		relax(y,f,count);

		//creating the intial guess for the next finer level
		//setting the dependent variables to zero for the current level
		//see "compute_initial_guess.c"
		init_guess(y,count);  
	}

	//solving the equation at the coarsest grids
	//relaxing on the grid using red-black gauss-seidel (see "relaxation.c")
	relax(y,f,count);
		
}
		
		  

				
			 	

	
					
		
										
		
					

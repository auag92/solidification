void relax(double *y,double *f,int count)
{
	double error;

	int start;

	int iter=0;

	printf("At level %d\n",count);

	//strating the iterations
	//while(iter<2 || error>toler)	
	while(iter<=num_iter-1)	
	{

		error=0.0;

		//============================================================================//
		//solving the equation for the even update only(see "red_black_solver.c")
		start=2;
		rb_solve(y,f,start,&error,count);
		//===========================================================================//	

		
		//============================================================================//
		//solving the equation for the odd update only(see "red_black_solver.c")
		start=1;
		rb_solve(y,f,start,&error,count);
		//===========================================================================//	
	
		//===========================================================================//
		error=sqrt(error*dx[count]);
		//===========================================================================//
	
		//===========================================================================//
		iter++;
		printf("iter=%d\terror=%lf\n",iter,error);
		//===========================================================================//		

	}
}

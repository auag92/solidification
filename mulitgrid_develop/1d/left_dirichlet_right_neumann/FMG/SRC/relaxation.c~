double relax(double *y,double *f,int mesh_tracker)
{
	double error=0.0;

	int begin;

	int iter=0;

	//printf("At level %d\n",mesh_tracker);
	
	//============================================================================//
	//solving the equation for the even update only(see "red_black_solver.c")
	begin=0;
	rb_solve(y,f,begin,&error,mesh_tracker);
	//===========================================================================//	

		
	//============================================================================//
	//solving the equation for the odd update only(see "red_black_solver.c")
	begin=1;
	rb_solve(y,f,begin,&error,mesh_tracker);
	//===========================================================================//	
	
	//===========================================================================//
	error=sqrt(error*dx[mesh_tracker]);
	//===========================================================================//
	
	return error;
	
}

void rhs(double *f)
{	
	int i;

	for(i=1;i<nodes_x;i++)//I am following the convention of "Briggs" where the equations at the boundary nodes at which Neumann bundary conditons are imposed are scaled by 2   
	{
		#ifdef ZERO_RHS
		f[i]=0.0;
		#endif

		#ifdef ONE_MIN_2_X_SQ_RHS
		f[i]=1.0-(2*(i*dx)*(i*dx));
		#endif
	}

	//scaling the end locations corresponding to the Neumann conidtion
	f[nodes_x-1]/=2.0; 
}						

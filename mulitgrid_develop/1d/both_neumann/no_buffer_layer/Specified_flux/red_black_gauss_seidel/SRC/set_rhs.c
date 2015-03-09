void rhs(double *f)
{	
	int i;

	for(i=0;i<nodes_x;i++)//I am initialising everywhere
	//I am following the convention of "Briggs" where the equations at the boundary nodes are scaled by 2   
	{
		#ifdef ZERO_RHS
		f[i]=0.0;
		#endif

		#ifdef ONE_MIN_2X_RHS
		f[i]=1.0-(2*i*dx);
		#endif
	}

	//scaling the end locations
	f[0]/=2.0;
	f[nodes_x-1]/=2.0; 
}						

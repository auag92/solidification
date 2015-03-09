void rhs(double *f)
{	
	int i;
	
	//I am going to initialise the rhs for the interior points on the finest mesh; indicated by mesh tracker value of 0 

	printf("rhs initialised for %ldth to %ldth\n",start[0],end[0]);
	
	for(i=start[0];i<=end[0];i++)//I am initialising everywhere
	//I am following the convention of "Briggs" where the equations at the boundary nodes are scaled by 2   
	{
		#ifdef ZERO_RHS //this will convert the poisson equation into a laplace equation  
		f[i]=0.0;
		#endif

		#ifdef ONE_MIN_2X_RHS
		f[i]=1.0-(2*(i-start[0])*dx[0]);
		#endif
	}
	//scaling the end locations
	f[start[0]]/=2.0;
	f[end[0]]/=2.0; 
}						

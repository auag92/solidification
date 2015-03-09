void rhs(double *f)
{	
	int i;
	
	//I am going to initialise the rhs for the interior points on the finest mesh 

	printf("rhs initialised for %ldth to %ldth\n",start[0]+1,end[0]-1);
	
	for(i=start[0]+1;i<=end[0]-1;i++)// I am initialising the 'f' array for the finest grid 
	//after leaving out the buffer points
	{
		#ifdef ZERO_RHS //this will convert the poisson equation into a laplace equation  
		f[i]=0.0;
		#endif

		#ifdef ONE_MIN_2X_RHS
		f[i]=1.0-(2*(i-(start[0]+1))*dx[0]);
		#endif
	}
}						

void rhs(double *f)
{	
	int i;

	printf("rhs initialised for %dth to %dth\n",left[0]+1,right[0]-1);
	
	for(i=left[0]+1;i<right[0];i++)// I am initialising the 'f' array for the finest grid 
	//after leaving out the boundaries
	{
		#ifdef ZERO_RHS
		f[i]=0.0;
		#endif

		#ifdef X_RHS
		f[i]=i*dx[0];
		#endif
	}
}						

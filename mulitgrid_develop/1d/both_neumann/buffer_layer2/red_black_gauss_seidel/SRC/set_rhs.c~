void rhs(double *f)
{	
	int i;

	for(i=1;i<nodes_x-1;i++)//I am not initialising the arrays on the buffer layer locations  
	{
		#ifdef ZERO_RHS
		f[i]=0.0;
		#endif

		#ifdef ONE_MIN_2X_RHS
		f[i]=1.0-(2*(i-1)*dx);
		#endif
	}
}						

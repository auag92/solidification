void rhs(double *f)
{	
	int i;

	for(i=0;i<nodes_x;i++)//though I am initialising the 'f' array at the boundary locations thwy won't be put to use
	{
		#ifdef ZERO_RHS
		f[i]=0.0;
		#endif

		#ifdef X_RHS
		f[i]=i*dx;
		#endif
	}
}						

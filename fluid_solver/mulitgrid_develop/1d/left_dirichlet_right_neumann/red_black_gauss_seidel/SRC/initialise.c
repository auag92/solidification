void init(double *y)
{
	int i;

	for(i=1;i<nodes_x;i++)//as we are going to solve for a Neumann boundary condition at the right end and a Dirichlet at the left end
	 
	{
		#ifdef ZERO
 		y[i]=0.0;
		#endif

		#ifdef RANDOM
		y[i]=(2.0*ran2(&SEED)-1.0)*noise_amplitude;
		#endif
	}
}

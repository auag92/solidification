void init(double *y)
{
	int i;

	for(i=1;i<nodes_x-1;i++)//leaving out the boundary points for the Dirichlet boundary condition 
	{
		#ifdef ZERO
 		y[i]=0.0;
		#endif

		#ifdef RANDOM
		y[i]=(2.0*ran2(&SEED)-1.0)*noise_amplitude;
		#endif
	}
}

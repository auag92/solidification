void init(double *y)
{
	int i;

	for(i=1;i<=nodes_x-2;i++)//a Dirichlet conditon at both ends
	//we will set the ends when we set the boundary condition 
	{
		#ifdef ZERO
 		y[i]=0.0;
		#endif

		#ifdef RANDOM
		y[i]=(2.0*ran2(&SEED)-1.0)*noise_amplitude;
		#endif
	}
}

void init(double *y)
{
	int i;

	//in principle I could have initialised anything for the finest grid and kept the other grids at zero

	for(i=0;i<nodes_x_sum;i++) //I have initialised everything including boundary points
	{
		#ifdef ZERO
 		y[i]=0.0;
		#endif

		#ifdef RANDOM
		y[i]=(2.0*ran2(&SEED)-1.0)*noise_amplitude;
		#endif
	}
}

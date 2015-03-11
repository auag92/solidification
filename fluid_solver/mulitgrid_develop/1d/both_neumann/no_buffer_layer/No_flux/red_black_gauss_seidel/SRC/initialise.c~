void init(double *y)
{
	int i;

	for(i=1;i<nodes_x-1;i++)//I am initialising everywhere in the physical domain
	//still the buffer layers have been left out
	{
		#ifdef ZERO
 		y[i]=0.0;
		#endif

		#ifdef RANDOM
		y[i]=(2.0*ran2(&SEED)-1.0)*noise_amplitude;
		#endif
	}
}

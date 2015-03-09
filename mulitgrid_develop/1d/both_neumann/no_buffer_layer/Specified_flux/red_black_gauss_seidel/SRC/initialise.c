void init(double *y)
{
	int i;

	for(i=0;i<nodes_x;i++)//I am initialising everywhere in the physical domain
	//as we are going to solve for a Neumann boundary condition, we need to solve the equation everywhere   
	{
		#ifdef ZERO
 		y[i]=0.0;
		#endif

		#ifdef RANDOM
		y[i]=(2.0*ran2(&SEED)-1.0)*noise_amplitude;
		#endif
	}
}

void init(double *y)
{
	int i,j;

	long array_index;	

	for(j=1;j<nodes.y-1;j++)//leaving out the boundary points for the Dirichlet boundary condition 
	{	
		for(i=1;i<nodes.x-1;i++)//leaving out the boundary points for the Dirichlet boundary condition 
		{
			array_index=i+nodes.x*j;
	
			#ifdef ZERO
 			y[array_index]=0.0;
			#endif

			#ifdef RANDOM
			y[array_index]=(2.0*ran2(&SEED)-1.0)*NOISE_AMPLITUDE;
			#endif
		}
	}
}

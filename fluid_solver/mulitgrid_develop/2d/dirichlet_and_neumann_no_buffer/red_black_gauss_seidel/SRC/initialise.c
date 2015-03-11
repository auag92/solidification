void init(double *y)
{
	int i,j;

	long array_index;	

	for(j=0;j<nodes.y-1;j++)//leaving out the boundary points for the Dirichlet boundary condition 
	{	
		for(i=0;i<nodes.x-1;i++)//leaving out the boundary points for the Dirichlet boundary condition 
		//but we have initialised where flux boundary conditions have to be applied
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

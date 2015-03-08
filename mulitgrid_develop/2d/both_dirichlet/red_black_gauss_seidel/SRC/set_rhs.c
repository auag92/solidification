void rhs(double *f)
{	
	int i,j;

	long array_index;

	#ifdef POINT_SOURCE

	long x_index= (long)(LOCATE_POINT_SOURCE_X*nodes.x);
		
	long y_index= (long)(LOCATE_POINT_SOURCE_Y*nodes.y);	

	printf("The point source is located at=(%ld,%ld)",x_index,y_index);

	#endif
			

	for(j=1;j<nodes.y-1;j++)//leaving out the boundary points for the Dirichlet boundary condition 
	{
		for(i=1;i<nodes.x-1;i++)//leaving out the boundary points for the Dirichlet boundary condition 
		{
			array_index=i+nodes.x*j;
		
			#ifdef ZERO_RHS
			f[array_index]=0.0;
			#endif

			#ifdef POINT_SOURCE
			if(i==x_index && j==y_index)
			{	
				f[array_index]=POINT_SOURCE_MAGNITUDE;
			}
			else
			{
				f[array_index]=0.0;
			}
			#endif

		}
	}
}						

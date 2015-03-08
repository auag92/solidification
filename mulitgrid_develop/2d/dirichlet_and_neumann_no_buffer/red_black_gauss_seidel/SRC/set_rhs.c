void rhs(double *f,double *g)
{	
	int i,j;

	long array_index;

	#ifdef POINT_SOURCE

	long x_index= (long)(LOCATE_POINT_SOURCE_X*nodes.x);
		
	long y_index= (long)(LOCATE_POINT_SOURCE_Y*nodes.y);	

	printf("The point source is located at=(%ld,%ld)",x_index,y_index);

	#endif
			

	for(j=0;j<=nodes.y-1;j++)
	{
		for(i=0;i<=nodes.x-1;i++)
		//we have not left out the locations where Dirichlet boundary conditions have to be applied
		//but we won't solve our equations there
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

			#ifdef FUNCTION_RHS
			f[array_index]=1.0;
			#endif
			
			g[array_index]=f[array_index];

		}
	}

												
}								

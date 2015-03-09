void dirich_boundary(double *y)
{
	int i,j;

	long array_index;

	//setting the dirichlet boundary condition
	//currently the dirichlet boundary conditions are at 'right' and 'up'
	
	//setting the one at 'right' 
	
	i=nodes.x-1;

	for(j=0;j<=nodes.y-1;j++)//it supersedes the flux BC at one corner dirichlet BC at another
	{
		array_index=i+nodes.x*j;
		
		y[array_index]=RIGHT;
	}

	//setting the one at "up"
	j=nodes.y-1;							
	
	for(i=0;i<=nodes.x-2;i++)//this supersedes the flux bc at the left corner and is secondary to the dirichelt bc at he right corner
	{
		array_index=i+nodes.x*j;
		
		y[array_index]=UP;
	}   
}		
  

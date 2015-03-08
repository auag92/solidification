void boundary(double *y,double *f)
{
	int i,j;

	long array_index;	
	
		
	//initialising the left boundary (a flux condition)
	for(j=0;j<nodes.y;j++)
	{
		i=0; //left

		array_index=i+nodes.x*j;		

		f[array_index]+=(LEFT/d_sp.x);
	}

	//initialising the right boundary (a dirichlet condition)
	for(j=0;j<nodes.y;j++)
	{
		i=nodes.x-1; //right 

		array_index=i+nodes.x*j;		

		y[array_index]=RIGHT;
	}

	//initialising the up boundary (a dirichlet condition)
	for(i=1;i<nodes.x-1;i++)//I am going to consider the corner points as a part of the vertical boundaries  
	{
		j=nodes.y-1; //up 

		array_index=i+nodes.x*j;		

		y[array_index]=UP;
	}
	
	//initialising the below boundary (a flux condition)
	for(i=1;i<nodes.x-1;i++) //I am going to consider the corner points as a part of the vertical boundaries  
	{
		j=0; //below 

		array_index=i+nodes.x*j;		

		f[array_index]+=(BELOW/d_sp.y);
	}
	
}

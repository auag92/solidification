void boundary(double *y)
{
	int i,j;

	long array_index;	
	
		
	//initialising the left boundary
	for(j=0;j<nodes.y;j++)
	{
		i=0; //left

		array_index=i+nodes.x*j;		

		y[array_index]=LEFT;
	}

	//initialising the right boundary
	for(j=0;j<nodes.y;j++)
	{
		i=nodes.x-1; //right 

		array_index=i+nodes.x*j;		

		y[array_index]=RIGHT;
	}

	//initialising the up boundary
	for(i=1;i<nodes.x-1;i++)//I am going to consider the corner points as a part of the vertical boundaries  
	{
		j=nodes.y-1; //up 

		array_index=i+nodes.x*j;		

		y[array_index]=UP;
	}
	
	//initialising the below boundary
	for(i=1;i<nodes.x-1;i++) //I am going to consider the corner points as a part of the vertical boundaries  
	{
		j=0; //below 

		array_index=i+nodes.x*j;		

		y[array_index]=BELOW;
	}
	
}

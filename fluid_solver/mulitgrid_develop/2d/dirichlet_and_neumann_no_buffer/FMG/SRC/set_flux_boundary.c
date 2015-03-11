void flux_boundary(double *f)
{
	int i,j;

	long array_index;

	//setting the different flux conditions for the finest grid only

	//setting the RHS to account for the different flux conditions


	//==========================================================//	
	//at the corner (0,0)
		
	i=0;
	
	j=0;

	array_index=(i+start[0])+nodes[0].x*j;
	
	f[array_index]/=2;

	f[array_index]+=(LEFT/d_sp[0].x)+(BELOW/d_sp[0].y);

	//===========================================================//

	//==========================================================//	
	//at the left end 
	
	i=0;

	for(j=1;j<=nodes[0].y-2;j++)//the corner (0,0) has already been tackled and the other corner has a dirichlet
	{
		array_index=(i+start[0])+nodes[0].x*j;	

		f[array_index]+=((2.0*LEFT)/d_sp[0].x);
	}
	//==========================================================//

	//==========================================================//	
	//at the below end 
	
	j=0;

	for(i=1;i<=nodes[0].x-2;i++)//the corner (0,0) has already been tackled and the other corner has a dirichlet
	{
		array_index=(i+start[0])+nodes[0].x*j;	

		f[array_index]+=((2.0*BELOW)/d_sp[0].y);
	}
	//==========================================================//
				
}

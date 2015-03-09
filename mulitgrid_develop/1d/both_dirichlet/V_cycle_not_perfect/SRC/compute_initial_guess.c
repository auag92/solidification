void init_guess(double *y,int count)
{
	int i;
	int array_index_2i,array_index_2i_plus_1;
	int array_index_i,array_index_i_plus_1;

	//interpolating the solutions
	for(i=0;i<nodes_x[count]-1;i++)//this is going over the coarser grid 
	//the range is selected by the algorithm 
	{
		//computing the different array indices
		//for the coarser grid
		array_index_i=i+left[count];
		array_index_i_plus_1=(i+1)+left[count];				

		//for the finer grid
		array_index_2i=left[count-1]+2*i;
		array_index_2i_plus_1=left[count-1]+2*i+1;		

		//interpolating to the finer grid
		y[array_index_2i]+=y[array_index_i]; 
		y[array_index_2i_plus_1]+=0.5*(y[array_index_i]+y[array_index_i_plus_1]);

		//setting the coarse grid to zero
		y[array_index_i]=0.0;
	}
} 
			

void ftoc(double *y,double *f,int mesh_tracker)
{
	double res_2i,res_2i_min_1,res_2i_plus_1;
	int array_index_2i,array_index_2i_min_1,array_index_2i_plus_1,array_index_2i_min_2,array_index_2i_plus_2;
	
	int mesh_inter_2i,mesh_inter_2i_min_1,mehs_inter_2i_plus_1,mesh_inter_2i_min_2,mesh_inter_2i_plus_2;
	
	int i;	

	int array_index_i;
	
	//I have to first initialise the coarse grid array to zero (because it might have a solution from a previous iteration
	for(i=0;i<=nodes_x[mesh_tracker+1]-1;i++)//I am iterating over the entire coarse grid array including the boundary locations
	{
		array_index_i=i+start[mesh_tracker+1];		
		 
		y[array_index_i]=0.0;
	}		
				
	//I am going to employ full weighting for the interpolation 
	for(i=1;i<=nodes_x[mesh_tracker+1]-2;i++)//the range of iterations is selected by the algorithm 
	//I am iterating over the points in the coarser grid	
	{
		//===========================================================================//
		//computing the residuals in the finer grid
		array_index_2i=start[mesh_tracker]+2*i;
		array_index_2i_min_1=start[mesh_tracker]+2*i-1;
		array_index_2i_plus_1=start[mesh_tracker]+2*i+1;
		
		array_index_2i_min_2=start[mesh_tracker]+2*i-2;
		array_index_2i_plus_2=start[mesh_tracker]+2*i+2;	

		res_2i=f[array_index_2i]-((y[array_index_2i_min_1]-2.0*y[array_index_2i]+y[array_index_2i_plus_1])/(dx[mesh_tracker]*dx[mesh_tracker]));

		res_2i_min_1=f[array_index_2i_min_1]-((y[array_index_2i_min_2]-2.0*y[array_index_2i_min_1]+y[array_index_2i])/(dx[mesh_tracker]*dx[mesh_tracker]));	

		res_2i_plus_1=f[array_index_2i_plus_1]-((y[array_index_2i]-2.0*y[array_index_2i_plus_1]+y[array_index_2i_plus_2])/(dx[mesh_tracker]*dx[mesh_tracker]));
		//===========================================================================//	


		//===========================================================================//		
		//moving the errors to the coarser grid using the full weighing operator  	
		array_index_i=i+start[mesh_tracker+1];
		f[array_index_i]=0.25*(res_2i_min_1+2.0*res_2i+res_2i_plus_1);
		//===========================================================================//
	
	}
}				
		

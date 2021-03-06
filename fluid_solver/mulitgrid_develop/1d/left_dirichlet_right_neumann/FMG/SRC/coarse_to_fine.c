void ctof(double *y,int mesh_tracker)
{
	int i;
	int array_index_2i,array_index_2i_plus_1;
	int array_index_i,array_index_i_plus_1;

	//interpolating the solutions

	//the first loop will map the entries over the entire coarse mesh to array-indices which are twice that of them in the fine mesh    
	for(i=0;i<nodes_x[mesh_tracker];i++) //I am including boundary points (in a domain having even no. of sub-intervals the boundaries are always at even locations)
	//in contrast to the dirichlet case as I am solving for the bounadry points also I need to update them from the the information from the coarser meshes     
	//for the Dirichlet they shouldn't get updated   
	{
		//for the coarser grid
		array_index_i=i+start[mesh_tracker];

		//for the finer grid
		array_index_2i=2*i+start[mesh_tracker-1];
		
		//interpolating		
		y[array_index_2i]+=y[array_index_i];
	} 
		
	//this second loop will just tackle odd array indices in the finer grid  
	for(i=0;i<nodes_x[mesh_tracker]-1;i++)//this is going over the coarser grid 
	//the end of the loop is decided by the fact that we use 'i+1' 
	{
		//computing the different array indices
	
		//for the coarser grid
		array_index_i=i+start[mesh_tracker];
		array_index_i_plus_1=(i+1)+start[mesh_tracker];				

		//for the finer grid
		array_index_2i_plus_1=(2*i+1)+start[mesh_tracker-1];		

		//interpolating to the finer grid
		y[array_index_2i_plus_1]+=0.5*(y[array_index_i]+y[array_index_i_plus_1]);

		
	}
} 
			

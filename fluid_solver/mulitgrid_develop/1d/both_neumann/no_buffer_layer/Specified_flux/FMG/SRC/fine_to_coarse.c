void ftoc(double *y,double *f,int mesh_tracker)
{
	double res_2i,res_2i_min_1,res_2i_plus_1;
	int array_index_2i,array_index_2i_min_1,array_index_2i_plus_1,array_index_2i_min_2,array_index_2i_plus_2;
	int i;	

	int array_index_i;
	
	//I have to first initialise the coarse grid array to zero (because it might have a solution from a previous iteration
	for(i=0;i<=nodes_x[mesh_tracker+1]-1;i++)//I am iterating over the entire coarse grid array including the boundary locations
	{
		array_index_i=i+start[mesh_tracker+1];		
		 
		y[array_index_i]=0.0;
	}		
				
	//I am going to employ full weighting for the interpolation 
	for(i=0;i<=nodes_x[mesh_tracker+1]-1;i++)//the range of iterations is selected by the algorithm 
	//I am iterating over the points in the coarser grid	
	{

		//===========================================================================//
		//computing the array indices in the finer grid
		//some of them might turn out to be negative but we won't use them
		array_index_2i=start[mesh_tracker]+2*i;
		array_index_2i_min_1=start[mesh_tracker]+2*i-1;
		array_index_2i_plus_1=start[mesh_tracker]+2*i+1;
			
		array_index_2i_min_2=start[mesh_tracker]+2*i-2;
		array_index_2i_plus_2=start[mesh_tracker]+2*i+2;
		//===========================================================================//

		if(i==0)
		{
			//===================================================================//
			//computing the residuals
		
			//=================================//	
			//the residual is computed using the equation at i=0, which is different in character from what used for the interior nodes  			
			res_2i=f[array_index_2i]-((-y[array_index_2i]+y[array_index_2i_plus_1])/(dx[mesh_tracker]*dx[mesh_tracker]));
			//================================//

			res_2i_plus_1=f[array_index_2i_plus_1]-((y[array_index_2i]-2.0*y[array_index_2i_plus_1]+y[array_index_2i_plus_2])/(dx[mesh_tracker]*dx[mesh_tracker]));

			
			//=========================================================================//	
			//moving the errors to the coarser grid using the full weighing operator  	
	
			//this equation is modified by dropping the term which corresponds to a negative index
			array_index_i=i+start[mesh_tracker+1];
			f[array_index_i]=0.25*(2.0*res_2i+res_2i_plus_1);
			//===========================================================================//
		}

		else if(i==nodes_x[mesh_tracker+1]-1)
		{
			//===================================================================//
			//computing the residuals
		
			//=================================//	
			//the residual is computed using the equation at i=nodes_x[mesh_tracker+1]-1, which is different in character from what used for the interior nodes  		
			res_2i=f[array_index_2i]-((y[array_index_2i_min_1]-y[array_index_2i])/(dx[mesh_tracker]*dx[mesh_tracker]));
			//================================//

			res_2i_min_1=f[array_index_2i_min_1]-((y[array_index_2i_min_2]-2.0*y[array_index_2i_min_1]+y[array_index_2i])/(dx[mesh_tracker]*dx[mesh_tracker]));	

			//=========================================================================//	
			//moving the errors to the coarser grid using the full weighing operator
			//this equation is modified by dropping the term which corresponds to a unphysical index  	
			array_index_i=i+start[mesh_tracker+1];
			f[array_index_i]=0.25*(res_2i_min_1+2.0*res_2i);
			//===========================================================================//
		}	


		else
		{
				
			res_2i=f[array_index_2i]-((y[array_index_2i_min_1]-2.0*y[array_index_2i]+y[array_index_2i_plus_1])/(dx[mesh_tracker]*dx[mesh_tracker]));

			res_2i_min_1=f[array_index_2i_min_1]-((y[array_index_2i_min_2]-2.0*y[array_index_2i_min_1]+y[array_index_2i])/(dx[mesh_tracker]*dx[mesh_tracker]));	

			res_2i_plus_1=f[array_index_2i_plus_1]-((y[array_index_2i]-2.0*y[array_index_2i_plus_1]+y[array_index_2i_plus_2])/(dx[mesh_tracker]*dx[mesh_tracker]));
			


			//=========================================================================//	
			//moving the errors to the coarser grid using the full weighing operator  	
			array_index_i=i+start[mesh_tracker+1];
			f[array_index_i]=0.25*(res_2i_min_1+2.0*res_2i+res_2i_plus_1);
			//===========================================================================//
		}
	
	}
}				
		

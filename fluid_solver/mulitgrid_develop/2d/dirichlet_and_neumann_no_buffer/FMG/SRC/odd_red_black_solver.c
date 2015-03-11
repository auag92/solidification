void odd_rb_solve(double *y,double *f,double *error,int mesh_tracker)
{
	int i,j;
	double y_star,new_val;
	
	long array_index,array_index_up,array_index_below,array_index_left,array_index_right;

	long odd_check;

	for(j=0;j<nodes[mesh_tracker].y-1;j++) //leaving out the boundary points for the Dirichlet boundary condition 
	{
		for(i=0;i<nodes[mesh_tracker].x-1;i++)//leaving out the boundary points for the Dirichlet boundary condition 
		{
			odd_check=i+j;

			if(odd_check%2!=0)//the sum of the indices are odd 	
			{
				//computing the different array indices
				//this has to be done considering that 'i' includes the starting indices 
				array_index=start[mesh_tracker]+i+nodes[mesh_tracker].x*j;
				array_index_left=(start[mesh_tracker]+i-1)+nodes[mesh_tracker].x*j;
				array_index_right=(start[mesh_tracker]+i+1)+nodes[mesh_tracker].x*j;
				array_index_up=start[mesh_tracker]+i+nodes[mesh_tracker].x*(j+1);
				array_index_below=start[mesh_tracker]+i+nodes[mesh_tracker].x*(j-1);
			
				if(i==0 && j==0)
					y_star=	corner_prefac[mesh_tracker]*(sq_dy[mesh_tracker]*y[array_index_right]+sq_dx[mesh_tracker]*y[array_index_up]-f[array_index]*sq_dx_times_sq_dy[mesh_tracker]);

				else if(i==0 && j!=0)				
					y_star=g_s_prefac[mesh_tracker]*(y[array_index_right]*2.0*sq_dy[mesh_tracker]+(y[array_index_up]+y[array_index_below])*sq_dx[mesh_tracker]-sq_dx_times_sq_dy[mesh_tracker]*f[array_index]);

				else if(i!=0 && j==0)				
					y_star=g_s_prefac[mesh_tracker]*(y[array_index_up]*2.0*sq_dx[mesh_tracker]+(y[array_index_left]+y[array_index_right])*sq_dy[mesh_tracker]-sq_dx_times_sq_dy[mesh_tracker]*f[array_index]);
			
				else
					y_star=g_s_prefac[mesh_tracker]*((y[array_index_left]+y[array_index_right])*sq_dy[mesh_tracker]+(y[array_index_up]+y[array_index_below])*sq_dx[mesh_tracker]-sq_dx_times_sq_dy[mesh_tracker]*f[array_index]);
				
				new_val=y[array_index]+OMEGA*(y_star-y[array_index]);	

				*error+=(new_val-y[array_index])*(new_val-y[array_index]);

				y[array_index]=new_val;
			}
		}
	}

	
}				
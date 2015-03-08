void rb_solve(double *y,double *f,int begin,double *error,int mesh_tracker)
{
	int i;
	double y_star,new_val;
	
	//'start' corresponds to the left boundary 
	//'end' corresponds to the right boundary
	//'begin' will have values either '1' or '2'
		

	for(i=begin+start[mesh_tracker];i<=end[mesh_tracker];i+=2)
	{
		if((i-start[mesh_tracker])==0)
			y_star=y[i+1]-dx[mesh_tracker]*dx[mesh_tracker]*f[i];

		else if((i-start[mesh_tracker])==nodes_x[mesh_tracker]-1)
			y_star=y[i-1]-dx[mesh_tracker]*dx[mesh_tracker]*f[i];	
		
		else
			y_star=0.5*(y[i-1]+y[i+1]-dx[mesh_tracker]*dx[mesh_tracker]*f[i]);
	
		new_val=y[i]+OMEGA*(y_star-y[i]);	

		*error+=(new_val-y[i])*(new_val-y[i]);

		y[i]=new_val;
	}

	
}				

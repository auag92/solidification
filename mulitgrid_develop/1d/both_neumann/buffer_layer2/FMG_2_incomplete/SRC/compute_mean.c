double comp_mean(double *y,int mesh_tracker)
{
	int i;
	
	double sum=0;
	
	for(i=1+start[mesh_tracker];i<=end[mesh_tracker]-1;i++)
	{
		sum+=y[i];
	}
	
	return(sum/(nodes_x[mesh_tracker]-2));
}

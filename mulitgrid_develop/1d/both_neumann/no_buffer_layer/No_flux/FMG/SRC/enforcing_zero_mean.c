void enf_zero_mean(double *y,int mesh_tracker)
{
	int i;
	double mean;

	//computing the mean (see "compute_mean.c")
	mean=comp_mean(y,mesh_tracker); 

	for(i=start[mesh_tracker];i<=end[mesh_tracker];i++)
	{
		y[i]-=mean;
	}
}			
	

void enf_zero_mean(double *y)
{
	int i;
	double mean;

	//computing the mean (see "compute_mean.c")
	mean=comp_mean(y); 

	for(i=1;i<nodes_x-1;i++)
	{
		y[i]-=mean;
	}
}			
	

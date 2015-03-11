void enf_zero_mean(double *y)
{
	int i;
	double mean;

	//computing the mean (see "compute_mean.c")
	mean=comp_mean(y); 

	for(i=0;i<nodes_x;i++)
	{
		y[i]-=mean;
	}
}			
	

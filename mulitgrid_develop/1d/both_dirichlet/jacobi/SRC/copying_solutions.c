void copy_soln(double *y,double *y_old)
{
	int i;
	for(i=0;i<nodes_x;i++)
	{
		y_old[i]=y[i];
	}
}

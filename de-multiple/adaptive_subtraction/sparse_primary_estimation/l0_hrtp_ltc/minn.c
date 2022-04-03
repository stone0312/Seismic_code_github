float minn(float *a, int n)
{
	int i;
	float min;
  
	min= a[0];

	for(i=1; i< n; i++)
	{
		if(a[i]< min)min=a[i];
	}

	return min;
}

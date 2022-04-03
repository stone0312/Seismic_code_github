float maxn(float *a, int n)
{
	int i;
	float max;
  
	max= a[0];

	for(i=1; i< n; i++)
	{
		if(a[i]> max)max=a[i];
	}

	return max;
}

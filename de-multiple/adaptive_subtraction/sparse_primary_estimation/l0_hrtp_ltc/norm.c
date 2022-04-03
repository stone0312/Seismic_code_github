float norm(float complex a)
{
	float num;

	 num= crealf(a)*crealf(a)+ cimagf(a)* cimagf(a);

	return num;
}

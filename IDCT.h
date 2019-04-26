#pragma once 
#include<math.h>
static const float cosine[8][8] = {
    {1.00000f, 0.98079f, 0.92388f, 0.83147f, 0.70711f, 0.55557f, 0.38268f, 0.19509f},
    {1.00000f, 0.83147f, 0.38268f, -0.19509f, -0.70711f, -0.98079f, -0.92388f, -0.55557f},
    {1.00000f, 0.55557f, -0.38268f, -0.98079f, -0.70711f, 0.19509f, 0.92388f, 0.83147f},
    {1.00000f, 0.19509f, -0.92388f, -0.55557f, 0.70711f, 0.83147f, -0.38268f, -0.98079f},
    {1.00000f, -0.19509f, -0.92388f, 0.55557f, 0.70711f, -0.83147f, -0.38268f, 0.98079f},
    {1.00000f, -0.55557f, -0.38268f, 0.98079f, -0.70711f, -0.19509f, 0.92388f, -0.83147f},
    {1.00000f, -0.83147f, 0.38268f, 0.19509f, -0.70711f, 0.98079f, -0.92388f, 0.55557f},
    {1.00000f, -0.98079f, 0.92388f, -0.83147f, 0.70711f, -0.55557f, 0.38268f, -0.19509f}
};

void IDCT(Block& res)
{
	static float c[8];
	c[0]=1.0/sqrt(2);
	for(int i=1;i<8;i++)c[i]=1.0;
	static float cc[8][8];
	static float tmp[8][8];
	for(int i=0;i<8;i++)
		for(int j=0;j<8;j++)cc[i][j]=c[i]*c[j];
	for(int i=0;i<8;i++)for(int j=0;j<8;j++)
	{
		tmp[i][j]=0.0;
		for(int n1=0;n1<8;n1++)for(int n2=0;n2<8;n2++)
			tmp[i][j]+=cc[n1][n2]*res.a[n1][n2]*cosine[i][n1]*cosine[j][n2];
		tmp[i][j]*=0.25;
	}
	for(int i=0;i<8;i++)for(int j=0;j<8;j++)
		res.a[i][j]=(int)tmp[i][j];
}

void IDCT2(Block& res)
{
	const int N=8;
	const int idx[]={0,7,1,6,2,5,3,4};
	const float pi=std::acos(-1);
 	static float c[8];
	c[0]=(float)sqrt(1.0/N);
	for(int i=1;i<N;i++)c[i]=(float)sqrt(2.0/N);
	cp w[4*N];
	for(int i=0;i<N;i++)w[i]=cp(cos(2*pi/(4*N)*i),sin(2*pi/(4*N)*i));
	for(int i=0;i<N;i++)w[i]=conj(w[i]);
	FFT solver;
	solver.init(N);
	float tmp[9][9]={0};
	for(int i=0;i<N;i++)for(int j=0;j<N;j++)tmp[i][j]=res.a[i][j];
	for(int i=0;i<N;i++)for(int j=0;j<N;j++)
	{
		tmp[i][j]/=(c[i]*c[j]);
	}
	for(int i=0;i<N;i++)
	{
		cp F[8];
		for(int j=0;j<N;j++)F[j]=cp(tmp[i][j],0)+cp(0,tmp[i][N-j]);
		for(int j=0;j<N;j++)F[j]*=w[j];
		solver.idft(F);
		for(int j=0;j<N;j++)tmp[i][j]=F[idx[j]].real();
	}
	for(int j=0;j<N;j++)
	{
		cp F[8];
		for(int i=0;i<N;i++)F[i]=cp(tmp[i][j],0)+cp(0,tmp[N-i][j]);
		for(int i=0;i<N;i++)F[i]*=w[i];
		solver.idft(F);
		for(int i=0;i<N;i++)tmp[i][j]=F[idx[i]].real();
	}
	for(int i=0;i<N;i++)for(int j=0;j<N;j++)res.a[i][j]=(int)tmp[i][j];
	/*
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)std::cout << std::setfill(' ') << std::setw(4) << res.a[i][j] << ' ';
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/
	}

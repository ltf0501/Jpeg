#pragma once
#include<math.h>
#include<complex>
#define cp std::complex<float>

struct FFT
{
	const double pi=std::acos(-1);
	int n,rev[16];
	cp omega[16],iomega[16];

	void init(int n)
	{
		this->n=n;
		for(int i=0;i<n;i++)
		{
			omega[i]=cp(cos(2*pi/n*i),sin(2*pi/n*i));
			iomega[i]=conj(omega[i]);
		}
		int k=log2(n);
		for(int i=0;i<n;i++)
		{
			int t=0;
			for(int j=0;j<k;j++)if(i&(1<<j))t|=(1<<(k-j-1));
			rev[i]=t;
		}
	}
	void transform(cp* a,cp* omega)
	{
		for(int i=0;i<n;i++)if(i<rev[i])swap(a[i],a[rev[i]]);
		for(int len=2;len<=n;len*=2)
		{
			int mid=len>>1;
			for(cp* p=a;p!=a+n;p+=len)
			{
				for(int i=0;i<mid;i++)
				{
					cp t=omega[n/len*i]*p[mid+i];
					p[mid+i]=p[i]-t;
					p[i]+=t;
				}
			}
		}
	}
	void dft(cp* a)
	{
		transform(a,omega);
	}
	void idft(cp* a)
	{
		transform(a,iomega);
		for(int i=0;i<n;i++)a[i]/=n;
	}
};

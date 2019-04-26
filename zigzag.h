#pragma once

void fillmatrix(int* buffer,Block& tmp)
{
	int now_x=0,now_y=0;
	int flag=-1;
	int* pt=buffer;
	for(int tot=0;tot<15;tot++)
	{
		int t=tot+1;
		if(tot>7)t=15-tot;
		for(int i=(flag==1 ? std::max(0,tot-7) : std::min(7,tot));t;t--,i+=flag)
		{
			int j=tot-i;
			tmp.a[i][j]=*(pt++);
		}
		flag=-flag;
	}
/*	
	for(int i=0;i<8;i++)
	{
		for(int j=0;j<8;j++)
		{
			std::cout << (unsigned int)tmp.a[i][j] << ' ';
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
*/	
}

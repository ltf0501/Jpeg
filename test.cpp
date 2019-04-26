#include<bits/stdc++.h>
using namespace std;
void read(uint8_t* buffer,int num,FILE* file)
{
	fread(buffer,sizeof(char),num,file);
}

int main(int argc,char *argv[])
{	
	FILE* file;
	uint8_t buffer[2];
	file=fopen(argv[1],"rb");	
	uint8_t tmp=0;
	bool b_flag=0;
	bool p_flag=0;
	while(1)
	{
		read(buffer,1,file);
		cout << setfill('0') << setw(2) << hex << (unsigned int)buffer[0] << ' ';
		if(buffer[0]==0xc0 && tmp==0xff)b_flag=1;
		if(buffer[0]==0xc2 && tmp==0xff)p_flag=1;
		if(buffer[0]==0xd9 && tmp==0xff)break ;
		tmp=buffer[0];
	}
	if(b_flag)cout << "~baseline~ " << endl;
	if(p_flag)cout << "~progressive~ " << endl;
	return 0;
}
/// NOW debug on line255

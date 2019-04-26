#pragma once
#include<cstdlib>
#include<cstring>
void bmp_write(std::vector<std::vector<RGB> >& graph,char* filename)
{
	int n=strlen(filename);
	filename[n-1]='p',filename[n-2]='m',filename[n-3]='b';
	FILE* output=fopen(filename,"wb");
	const int header_size=54;
	char header[header_size]={0x42,0x4d,0,0,0,0,0,0,0,0,54,0,0,0,40,0,0,0,0,0,0,0,0,0,0,0,1,0,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int row=graph.size();
	int column=graph[0].size();
	long file_size=(long)row*(long)column*3+header_size;
	header[2]=(unsigned char)(file_size & 0xff);
	header[3]=(file_size>>8) & 0xff;
	header[4]=(file_size>>16) & 0xff;
	header[5]=(file_size>>24) & 0xff;
	
	header[18]=column & 0xff;
	header[19]=(column>>8) & 0xff;
	header[20]=(column>>16) & 0xff;
	header[21]=(column>>24) & 0xff;
	
	header[22]=row & 0xff;
	header[23]=(row>>8) & 0xff;
	header[24]=(row>>16) & 0xff;
	header[25]=(row>>24) & 0xff;
	write(header,header_size,output);
	for(int i=row-1;i>=0;i--)for(int j=0;j<column;j++)
	{
		char buffer[3];
		buffer[0]=graph[i][j].b;
		buffer[1]=graph[i][j].g;
		buffer[2]=graph[i][j].r;
		write(buffer,3,output);	
	}
	fclose(output);
}

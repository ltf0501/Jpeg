#include<cstdio>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>

#include "function.h"
#include "read.h"
#include "zigzag.h"
#include "huffman.h"
#include "FFT.h"
#include "IDCT.h"
#include "BMPoutput.h"
int cal_len(uint8_t* buffer)
{
	int res=buffer[0];
	res<<=8;
	res|=buffer[1];
	return res;
}
void garbage(FILE* file)
{
	uint8_t buffer[2];
	read(buffer,2,file);
	int len=cal_len(buffer);
	uint8_t* gar;
	gar=(uint8_t *)malloc(len-2);
	read(gar,len-2,file);
	check("gar");
}
void readDQT(FILE* file,std::vector<Block>& qtable)
{
	uint8_t buffer[129];
	read(buffer,2,file);
	int len=cal_len(buffer);
	int cnt=2;
	while(cnt<len)
	{
		read(buffer,1,file),cnt++;
		int num=buffer[0]>>4;
		int id=buffer[0] & 0xf;
		num=64*(num+1);
		read(buffer,num,file),cnt+=num;
		Block tmp;
		int buffer1[129];
		for(int i=0;i<num;i++)buffer1[i]=buffer[i];
		fillmatrix(buffer1,tmp);
		qtable[id]=tmp;
	}
	check("DQT");
}
void readSOF0(FILE* file,SOF0data& image)
{
	uint8_t buffer[3];
	read(buffer,2,file);
	int len=cal_len(buffer);
	read(buffer,1,file);
	read(buffer,2,file);
	image.row=cal_len(buffer);
	read(buffer,2,file);
	image.column=cal_len(buffer);
	read(buffer,1,file);
	image.col_num=buffer[0];
	image.Hmax=image.Vmax=0;
	for(int i=0;i<image.col_num;i++)
	{
		read(buffer,3,file);
		image.id.push_back(buffer[0]);
		int HH=buffer[1]>>4;
		int VV=buffer[1] & 0xf;
		image.Hmax=std::max(image.Hmax,HH);
		image.Vmax=std::max(image.Vmax,VV);
		image.horizontal.push_back(HH);
		image.vertical.push_back(VV);
		image.qtable_id.push_back(buffer[2]);
	}
	check("SOF0");
}
void readDHT(FILE* file,std::vector<std::vector<HuffmanTable> >& htable)
{
	uint8_t* buffer;
	buffer=(uint8_t*)malloc(2);
	read(buffer,2,file);
	int len=cal_len(buffer);
	int cur_len=len-2;
	buffer=(uint8_t*)malloc(17);
	// read huffman table
	while(cur_len)
	{
		read(buffer,17,file);
		cur_len-=17;
		bool is_AC=buffer[0]>>4;
		int id=buffer[0] & 0xf;
		int tot=0;
		for(int i=1;i<=16;i++)tot+=buffer[i];
		uint8_t* codes;
		codes=(uint8_t*)malloc(tot);
		read(codes,tot,file);
		cur_len-=tot;
		htable[is_AC][id]=construct(buffer,codes);
	}
	check("DHT");
}
void readDRI(FILE* file)
{
	uint8_t buffer[2];
	read(buffer,2,file);
	read(buffer,2,file);
	check("DRI");
}
bool read0xff(FILE* file,std::string& s)
{
	uint8_t buffer[2];
	while(1)
	{
		read(buffer,1,file);
		if(buffer[0]!=0xff)break ;
	}
	if(buffer[0]==0x00){s+="11111111";return 0;}
	if(buffer[0]==0xd9)return 1;
	if(0xd0<=buffer[0] && buffer[0]<=0xd7)return 0;
	for(int i=7;i>=0;i--)
	{
		if(buffer[0] && (1<<i))s+='1';
		else s+='0';
	}
	return 0;
}
Block make_table(int id,std::vector<std::vector<HuffmanTable> > htable,std::vector<int>& DC_id,std::vector<int>& AC_id,std::string& s)
{
	int tot=(int)s.size();
	static int now=0;
	Block res;
	int buffer[65]={0},p=1;
	int cur=0;
	int len=0;
	HuffmanTable& tree=htable[0][DC_id[id]];
	
	//read DC
	
	while(1)
	{
		cur=(cur<<1)+(s[now++]-'0');len++;
		if(cur-tree.mi[len]<tree.cnt[len])
		{
			int tt=tree.mp[len][cur-tree.mi[len]];
			int T=tt;
			if(!T)break ;
			int flag=s[now]-'0';T--;
			buffer[0]=(buffer[0]<<1)+(s[now++]-'0');
			while(T--)
			{
				buffer[0]=(buffer[0]<<1)+(s[now++]-'0');
			}
			int x=1;
			if(!flag)buffer[0]-=((1<<tt)-1);
			break ;
		}
	}
	tree=htable[1][AC_id[id]];
	cur=0;len=0;

	//read AC
	
	while(1)
	{
		cur=(cur<<1)+(s[now++]-'0');len++;	
		if(cur-tree.mi[len]<tree.cnt[len])
		{
			int tt=tree.mp[len][cur-tree.mi[len]];
			if(tt==0)break ;
			p+=(tt>>4);
			int nxt=tt & 0xf;
			int T=nxt;
			if(T)
			{
				int flag=s[now]-'0';T--;
				buffer[p]=(buffer[p]<<1)+(s[now++]-'0');
				while(T--)
				{
					buffer[p]=(buffer[p]<<1)+(s[now++]-'0');
				}
				if(!flag)buffer[p]-=((1<<nxt)-1);
			}
			p++;
			cur=0,len=0;
			if(p>=64)break ;
		}
	}
	fillmatrix(buffer,res);
	return res;
}
void readSOS(FILE* file,SOF0data image,std::vector<std::vector<HuffmanTable> > htable,std::vector<Block>& qtable,std::vector< std::vector<RGB> >& graph)
{
	uint8_t buffer[2];
	read(buffer,2,file);
	int len=cal_len(buffer);
	read(buffer,1,file);
	int col_num=buffer[0];
	std::vector<int> DC_id(col_num),AC_id(col_num);
	
	//read DC_id AC_id
	
	for(int i=0;i<col_num;i++)
	{
		read(buffer,2,file);
		DC_id[buffer[0]]=buffer[1]>>4;
		AC_id[buffer[0]]=buffer[1] & 0xf;
	}
	read(buffer,3,file);
	
	//read picture 
	std::string s;
	while(1)
	{
		read(buffer,1,file);
		if(buffer[0]!=0xff)
		{
			for(int i=7;i>=0;i--)
			{
				if(buffer[0] & (1<<i))s+='1';
				else s+='0';
			}
		}
		else 
		{
			if(read0xff(file,s))break ;
		}
	}
	// picture basic info
	int row=image.row;
	int column=image.column;
	int b_row=(row+8*image.Vmax-1)/(8*image.Vmax);
	int b_column=(column+8*image.Hmax-1)/(8*image.Hmax);
	row=8*image.Vmax*b_row,column=8*image.Hmax*b_column;
	int rH[3],rV[3];
	for(int i=0;i<3;i++)rH[i]=image.Hmax/image.horizontal[i],rV[i]=image.Vmax/image.vertical[i];
	graph.resize(row,std::vector<RGB>(column));
	
	// read Mcu 
	int pre[4]={0};
	IDCT_init();
	for(int i=0;i<b_row;i++)for(int j=0;j<b_column;j++)
	{
		Mcu mcudata;
		for(int num=0;num<3;num++)
		{
			int id=image.id[num];
			mcudata.v[id].resize(image.vertical[num],std::vector<Block>(image.horizontal[num]));
			int q_id=image.qtable_id[num];		
			for(int x=0;x<image.vertical[num];x++)for(int y=0;y<image.horizontal[num];y++)
			{
				mcudata.v[id][x][y]=make_table(id,htable,DC_id,AC_id,s);
				mcudata.v[id][x][y].a[0][0]+=pre[id];
				pre[id]=mcudata.v[id][x][y].a[0][0];
				for(int aa=0;aa<8;aa++)
				{
					for(int bb=0;bb<8;bb++)mcudata.v[id][x][y].a[aa][bb]*=qtable[q_id].a[aa][bb];
				}	
				IDCT2(mcudata.v[id][x][y]);
			}
		}
		//YCbCr to RGB
		for(int x=0;x<8*image.Vmax;x++)for(int y=0;y<8*image.Hmax;y++)
		{
			int Y=mcudata.v[image.id[0]][x/(8*rV[0])][y/(8*rH[0])].a[(x/rV[0])%8][(y/rH[0])%8];
			int Cb=mcudata.v[image.id[1]][x/(8*rV[1])][y/(8*rH[1])].a[(x/rV[1])%8][(y/rH[1])%8];
			int Cr=mcudata.v[image.id[2]][x/(8*rV[2])][y/(8*rH[2])].a[(x/rV[2])%8][(y/rH[2])%8];

			float R=Y+1.402*Cr+128;
			float G=Y-0.3441416*Cb-0.714136*Cr+128;
			float B=Y+1.772*Cb+128;
			if(R>=255)R=255;
			if(R<=0)R=0;
			if(G>=255)G=255;
			if(G<=0)G=0;
			if(B>=255)B=255;
			if(B<=0)B=0;
			RGB& tmp=graph[8*image.Vmax*i+x][8*image.Hmax*j+y];
			tmp=(RGB){(unsigned int)R,(unsigned int)G,(unsigned int)B};
		}
	}
	check("SOS");
}
int main(int argc,char *argv[])
{	
	FILE* file;
	uint8_t buffer[2];
	if(argc==1)
	{
		puts("No input.");
		return 0;
	}
	file=fopen(argv[1],"rb");
	read(buffer,2,file);
	if(buffer[0]!=start || buffer[1]!=SOI)
	{
		throw std::runtime_error("SOI marker not found");
	}
	std::vector<Block> qtable(4);
	SOF0data image;
	std::vector<std::vector<HuffmanTable> > htable(2,std::vector<HuffmanTable>(4));
	std::vector<std::vector<RGB> > graph;
	bool eoi=0;
	while(!eoi)
	{
		read(buffer,2,file);
		uint8_t marker=buffer[1];
		if(APPO_MIN <= marker && marker <=APPO_MAX){garbage(file);continue;}
		//check the marker
		switch(marker)
		{
			case DQT:
				readDQT(file,qtable);
				break;
			case SOF0:
				readSOF0(file,image);
				break;
			case DHT:
				readDHT(file,htable);
				break; 
			case DRI:
				readDRI(file);
				break;
			case COM:
				break;
			case SOS:
				readSOS(file,image,htable,qtable,graph);
				bmp_write(graph,argv[1]);
				eoi=1;
				break;
		}
	}
	fclose(file);
	return 0;
}

#include<cstdio>
#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<map>

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
//		std::cout << "id = " << (unsigned int)buffer[0] << std::endl;
		int num=buffer[0]>>4;
		int id=buffer[0] & 0xf;
		num=64*(num+1);
		read(buffer,num,file),cnt+=num;
		Block tmp;
		int buffer1[129];
		for(int i=0;i<num;i++)buffer1[i]=buffer[i];
		fillmatrix(buffer1,tmp);
//		std::cout << "iid = " << (unsigned int)buffer[0] << std::endl;
		qtable[id]=tmp;
	}
	check("DQT");
}
void readSOF0(FILE* file,SOF0data& image)
{
	uint8_t buffer[3];
	read(buffer,2,file);
	int len=cal_len(buffer);
	std::cout << "SOF0: len = " << len << std::endl;
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
	
	std::cout << "row = " << (unsigned int)image.row << " column = " << (unsigned int)image.column << " col_num = " << image.col_num << std::endl;
	std::cout << "image.id: ";for(int i : image.id)std::cout << i << ' ';std::cout << std::endl;
	std::cout << "image.H: ";for(int i : image.horizontal)std::cout << i <<' ';std::cout <<std::endl;
	std::cout << "image.V: ";for(int i : image.vertical)std::cout << i <<' ';std::cout <<std::endl;
	std::cout << "image.qtable_id: "; for(int i : image.qtable_id)std::cout << i <<' ';std::cout <<std::endl;
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
		std::cout << "is_AC = " << is_AC << ",id = " << id << std::endl;
		std::cout << "num: ";
		for(int i=1;i<=16;i++)std::cout << (unsigned int)buffer[i] << ' ';std::cout << std::endl; 
		std::cout << "symbol: ";
		for(int i=0;i<tot;i++)std::cout << (unsigned int)codes[i] << ' ';std::cout << std::endl; 
		cur_len-=tot;
		htable[is_AC][id]=construct(buffer,codes);
		//free(buffer);
	}
	check("DHT");
}
void readDRI(FILE* file)
{
	uint8_t buffer[2];
	read(buffer,2,file);
	std::cout << cal_len(buffer) << std::endl;
	read(buffer,2,file);
	std::cout << cal_len(buffer) << std::endl;
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
	if(0xd0<=buffer[0] && buffer[0]<=0xd7)return 0;//RSTn tag, deal later
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
//	std::cout << "is_AC = 0, DC_id[" << id << "] = " << DC_id[id] << std::endl;  
	
	//read DC
	
	while(1)
	{
		cur=(cur<<1)+(s[now++]-'0');len++;
//		std::cout << "XD@@ len = " << len << ", cur = " << cur << std::endl;
		if(tree.mp[len].find(cur)!=tree.mp[len].end())
		{
			int tt=tree.mp[len][cur];
	//		std::cout << "len = " << len << ", cur = " << cur << ", tt = " << tt << std::endl;
			int T=tt;
			if(!T)break ;
			std::string ss;
			ss+=s[now];
			int flag=s[now]-'0';T--;
			buffer[0]=(buffer[0]<<1)+(s[now++]-'0');
			while(T--)
			{
				ss+=s[now];
				buffer[0]=(buffer[0]<<1)+(s[now++]-'0');
			}
			int x=1;
			if(!flag)buffer[0]-=((1<<tt)-1);
//			std::cout << "now = " << buffer[0] << ", ss = " << ss << std::endl;
			break ;
		}
	}
	tree=htable[1][AC_id[id]];
//	std::cout << "is_AC = 1, AC_id[" << id << "] = " << AC_id[id] << std::endl;  
	cur=0;len=0;
	std::string ss;

	//read AC
	
	while(1)
	{
	//	ss+=s[now];
		cur=(cur<<1)+(s[now++]-'0');len++;	
//		std::cout << "XD@@@ len = " << len << ", cur = " << cur << std::endl;
		if(tree.mp[len].find(cur)!=tree.mp[len].end())
		{
			int tt=tree.mp[len][cur];
		//	std::cout << "len = " << len << ", cur = " << cur << ", tt = " << tt << std::endl;
			if(tt==0)break ;
			p+=(tt>>4);
			int nxt=tt & 0xf;
		//	std::cout << "nxt = " << nxt << std::endl;
			int T=nxt;
			if(T)//it might be tt==15*16
			{
			//	ss+=s[now];
				int flag=s[now]-'0';T--;
				buffer[p]=(buffer[p]<<1)+(s[now++]-'0');
				while(T--)
				{
			//		ss+=s[now];
					buffer[p]=(buffer[p]<<1)+(s[now++]-'0');
				}
				if(!flag)buffer[p]-=((1<<nxt)-1);
			}
			p++;
		//	ss.clear();
			cur=0,len=0;
			if(p>64)
			{
				std::runtime_error("Bomb");
				exit(0);
			}
			if(p>=64)break ;
		}
	}
/*	
	std::cout <<"now id is " << id <<std::endl;
	for(int i=0;i<64;i++)
	{
		std::cout << std::setfill(' ') << std::setw(3) << buffer[i] << ((i%8==7) ? '\n' : ' '); 
	}
*/	
//	check("XD");
	
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
		std::cout << "DC_id[" << (unsigned int)buffer[0] << "] = " << DC_id[buffer[0]] << std::endl;
		std::cout << "AC_id[" << (unsigned int)buffer[0] << "] = " << AC_id[buffer[0]] << std::endl;
	}
	read(buffer,3,file);
	std::cout << (unsigned int)buffer[0] << ' ' << (unsigned int)buffer[1] << ' ' << (unsigned int)buffer[2] << std::endl;
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
//	std::cout << s << std::endl;
	int row=image.row;
	int column=image.column;
	//8*image.Vmax-1
	int b_row=(row+8*image.Vmax-1)/(8*image.Vmax);
	int b_column=(column+8*image.Hmax-1)/(8*image.Hmax);
	std::cout << "b_row = " << b_row << " b_column = " << b_column << std::endl;
	row=8*image.Vmax*b_row,column=8*image.Hmax*b_column;
	int rH[3],rV[3];
	for(int i=0;i<3;i++)rH[i]=image.Hmax/image.horizontal[i],rV[i]=image.Vmax/image.vertical[i];
	graph.resize(row,std::vector<RGB>(column));
	
	// read Mcu 
	int pre[4]={0};
	std::vector<std::vector<Mcu> > mcudata(b_row,std::vector<Mcu>(b_column));	
	for(int i=0;i<b_row;i++)for(int j=0;j<b_column;j++)
	{
	//	std::cout << "i = " << i << " j = " << j << std::endl;
		for(int num=0;num<3;num++)
		{
			int id=image.id[num];
			mcudata[i][j].v[id].resize(image.vertical[num],std::vector<Block>(image.horizontal[num]));
			int q_id=image.qtable_id[num];		
	//		std::cout << "id = " << id << std::endl;	
			for(int x=0;x<image.vertical[num];x++)for(int y=0;y<image.horizontal[num];y++)
			{
	//			std::cout << "i = " << i << " j = " << j << " id = " << id << " x = " << x << " y = " << y << std::endl;
				mcudata[i][j].v[id][x][y]=make_table(id,htable,DC_id,AC_id,s);
				mcudata[i][j].v[id][x][y].a[0][0]+=pre[id];
				pre[id]=mcudata[i][j].v[id][x][y].a[0][0];
	/*			
				std::cout << "zigzag:" << std::endl;
				for(int aa=0;aa<8;aa++)
				{
					for(int bb=0;bb<8;bb++)std::cout << std::setfill(' ') << std::setw(4) << mcudata[i][j].v[id][x][y].a[aa][bb] << ' ';
					std::cout << std::endl;
				}
				std::cout << std::endl;
	*/				
				for(int aa=0;aa<8;aa++)
				{
					for(int bb=0;bb<8;bb++)mcudata[i][j].v[id][x][y].a[aa][bb]*=qtable[q_id].a[aa][bb];
					/*			
					for(int bb=0;bb<8;bb++)std::cout << std::setfill(' ') << std::setw(4) << block[id][i+x][j+y].a[aa][bb] << ' ';
					std::cout << std::endl;
					*/
				}	
				IDCT2(mcudata[i][j].v[id][x][y]);
				
		/*
				for(int aa=0;aa<8;aa++)
				{
					for(int bb=0;bb<8;bb++)std::cout << std::setfill(' ') << std::setw(4) << mcudata[i][j].v[id][x][y].a[aa][bb] << ' ';
					std::cout << std::endl;
				}
				std::cout << std::endl;
		*/
			}
		}
		int Y[8*image.Vmax][8*image.Hmax];
		int Cr[8*image.Vmax][8*image.Hmax];
		int Cb[8*image.Vmax][8*image.Hmax];
		for(int x=0;x<8*image.Vmax;x++)for(int y=0;y<8*image.Hmax;y++)
		{
			Y[x][y]=mcudata[i][j].v[image.id[0]][x/(8*rV[0])][y/(8*rH[0])].a[(x/rV[0])%8][(y/rH[0])%8];
			Cb[x][y]=mcudata[i][j].v[image.id[1]][x/(8*rV[1])][y/(8*rH[1])].a[(x/rV[1])%8][(y/rH[1])%8];
			Cr[x][y]=mcudata[i][j].v[image.id[2]][x/(8*rV[2])][y/(8*rH[2])].a[(x/rV[2])%8][(y/rH[2])%8];

			float R=Y[x][y]+1.402*Cr[x][y]+128;
			float G=Y[x][y]-0.3441416*Cb[x][y]-0.714136*Cr[x][y]+128;
			float B=Y[x][y]+1.772*Cb[x][y]+128;
			if(R>=255)R=255;
			if(R<=0)R=0;
			if(G>=255)G=255;
			if(G<=0)G=0;
			if(B>=255)B=255;
			if(B<=0)B=0;
			RGB& tmp=graph[8*image.Vmax*i+x][8*image.Hmax*j+y];
			tmp=(RGB){(unsigned int)R,(unsigned int)G,(unsigned int)B};
		}
	/*	
		std::cout << "i = " << i << " j = " << j << std::endl;
		std :: cout << "Y:" <<std::endl;
		for(int x=0;x<8*image.Vmax;x++)
		{
			for(int y=0;y<8*image.Hmax;y++)std::cout << std::setfill(' ') << std::setw(4) << Y[x][y] << ' ';
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std :: cout << "Cb:" <<std::endl;
		for(int x=0;x<8*image.Vmax;x++)
		{
			for(int y=0;y<8*image.Hmax;y++)std::cout <<std::setfill(' ') << std::setw(4) <<  Cb[x][y] << ' ';
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std :: cout << "Cr:" <<std::endl;
		for(int x=0;x<8*image.Vmax;x++)
		{
			for(int y=0;y<8*image.Hmax;y++)std::cout << std::setfill(' ') << std::setw(4) << Cr[x][y] << ' ';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	*/	
	/*	
		std::cout << "i = " << i << " j = " << j << std::endl;
		std :: cout << "R:" <<std::endl;
		for(int x=0;x<8*image.Vmax;x++)
		{
			for(int y=0;y<8*image.Hmax;y++)std::cout << std::setfill(' ') << std::setw(4) << graph[8*image.Vmax*i+x][8*image.Hmax*j+y].r << ' ';
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std :: cout << "G:" <<std::endl;
		for(int x=0;x<8*image.Vmax;x++)
		{
			for(int y=0;y<8*image.Hmax;y++)std::cout << std::setfill(' ') << std::setw(4) << graph[8*image.Vmax*i+x][8*image.Hmax*j+y].g << ' ';
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std :: cout << "B:" <<std::endl;
		for(int x=0;x<8*image.Vmax;x++)
		{
			for(int y=0;y<8*image.Hmax;y++)std::cout << std::setfill(' ') << std::setw(4) << graph[8*i+x][8*j+y].b << ' ';
			std::cout << std::endl;
		}
		std::cout << std::endl;
	*/	
	}
	check("SOS");
}
int main(int argc,char *argv[])
{	
	FILE* file;
	uint8_t buffer[2];
	file=fopen(argv[1],"rb");
	if(argc==1)
	{
		std::cout << "No input." << std::endl;
		return 0;
	}
	int len=ftell(file);
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
		std::cout << std::hex <<  (unsigned int)buffer[1] << std::endl;
		std::cout << std::dec;
		if(APPO_MIN <= marker && marker <=APPO_MAX){garbage(file);continue;}
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
				if(argc==2)
				{
					std::cout << "No output filename." << std::endl;
					return 0;
				}
				bmp_write(graph,argv[2]);
				eoi=1;
				break;
		}
	}
	fclose(file);
	return 0;
}

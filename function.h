#pragma once
const int maxw=2000;
const int maxh=1500;
const uint8_t start=0xff;
const uint8_t APPO_MIN=0xe0;
const uint8_t APPO_MAX=0xef;
const uint8_t SOI=0xd8;
const uint8_t DQT=0xdb;
const uint8_t SOF0=0xc0;
const uint8_t DHT=0xc4;
const uint8_t DRI=0xdd;
const uint8_t SOS=0xda;
const uint8_t COM=0xfe;
const uint8_t EOI=0xd9;
struct Block
{
	int a[8][8];
};
struct Mcu
{
	std::vector<std::vector<Block> > v[4];
};
struct SOF0data
{
	int row,column;
	int col_num;
	std::vector<int> id,horizontal,vertical,qtable_id;
	int Hmax,Vmax;
};
struct HuffmanTable
{
	int cnt[17],mi[17];
	std::vector<int> mp[17];
};
struct RGB
{
	unsigned int r,g,b;	
};
void check(std::string s)
{
//	std::cout << "I AM A CHECKER " << s << std::endl;
}
void fillmatrix(int* buffer,Block& tmp);
void read(uint8_t* buffer,int num,FILE* file);
void write(char* c,int num,FILE* file);
HuffmanTable construct(uint8_t* count,uint8_t* codes);
void IDCT(Block& res);

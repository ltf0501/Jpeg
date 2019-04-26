#pragma once
void read(uint8_t* buffer,int num,FILE* file)
{
	fread(buffer,sizeof(char),num,file);
}
void write(char* c,int num,FILE* file)
{
	fwrite(c,sizeof(char),num,file);
}


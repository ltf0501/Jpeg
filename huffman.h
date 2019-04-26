#pragma once 
HuffmanTable construct(uint8_t* count,uint8_t* codes)
{
	HuffmanTable res;
	for(int i=1;i<=16;i++)res.cnt[i]=count[i];
	for(int i=1;i<=16;i++)res.mi[i]=-10000;
	uint8_t *pt=codes;
	int pre_val=-1,pre_i=1;
	for(int i=1;i<=16;i++)
	{
		int cnt=count[i];
		while(cnt--)
		{
			pre_val++;
			if((pre_val & ((1<<pre_i)-1))==0)pre_i++;
			pre_val<<=(i-pre_i);
			pre_i=i;
			if(cnt==res.cnt[i]-1)res.mi[i]=pre_val;
//			res.mp[i][pre_val]=*(pt++);
			res.mp[i].push_back(*(pt++));
		}
	}
	
	for(int i=1;i<=16;i++)
	{
		std::cout << "mi[" << i << "] = " << res.mi[i] << std::endl;
		for(int j=0;j<(int)res.mp[i].size();j++)std::cout << "mp[" << i << "][" << j << "] = " << res.mp[i][j] << std::endl;
	}
	std::cout << "XD" << std::endl;
	
	check("construct table");
	return res;
}

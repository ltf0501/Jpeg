all:compile
compile: decoder.cpp 
	g++ -O2 decoder.cpp -o decoder.out

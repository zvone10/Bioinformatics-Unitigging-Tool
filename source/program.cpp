#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "overlapgraph.h"

std::string path;
std::vector<std::string> reads;
std::vector<overlap> overlaps;

std::vector<std::string> split(std::string input, char splitChar)
{
	std::vector<std::string> substrings;
	std::string::iterator it;
	std::string s = "";
	for (it = input.begin(); it<input.end(); it++)
	{
		if (*it == splitChar)
		{
			substrings.push_back(s);
			s = "";
		}
		else
			s += *it;
	}
	substrings.push_back(s);
	return substrings;
}

int main(int argc, char *argv[])
{
	path = argv[1];
	
	std::string readsfile = path + "\\reads.txt";
	std::string overlapsfile = path + "\\overlaps.txt";
	
	std::ifstream readinput(readsfile);
	for (std::string line; std::getline(readinput, line);)
	{
		std::vector<std::string> substrings = split(line, '\t');
		reads.push_back(substrings[1]);
	}
	
	//std::cout << reads.size() << std::endl;
	
	std::ifstream overlapinput(overlapsfile);
	for (std::string line; getline(overlapinput, line);)
	{
		std::vector<std::string> substrings = split(line, '\t');
		overlap o;
		o.read1 = std::stoi(substrings[0]);
		o.read2 = std::stoi(substrings[1]);
		o.ahg = std::stoi(substrings[2]);
		o.bhg = std::stoi(substrings[3]);
		o.orientation = substrings[4][0];
		
		overlaps.push_back(o);
	}
	
}

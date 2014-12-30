#include <stdio.h>
#include <string>
#include <cstdlib>
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

	std::ifstream readinput(readsfile.c_str());
	for (std::string line; std::getline(readinput, line);)
	{
		std::vector<std::string> substrings = split(line, '\t');
		reads.push_back(substrings[1]);
	}

	//std::cout << reads.size() << std::endl;

	std::ifstream overlapinput(overlapsfile.c_str());
	for (std::string line; getline(overlapinput, line);)
	{
		std::vector<std::string> substrings = split(line, '\t');
		overlap o;
		o.read1 = (int)strtol(substrings[0].c_str(),0,10) - 1;
		o.read2 = (int)strtol(substrings[1].c_str(),0,10) - 1;
		o.ahg = (int)strtol(substrings[2].c_str(),0,10);
		o.bhg = (int)strtol(substrings[3].c_str(),0,10);
		o.orientation = substrings[4][0];

		overlaps.push_back(o);
	}

	OverlapGraph o(reads, overlaps);
	o.runUnitigging();

}

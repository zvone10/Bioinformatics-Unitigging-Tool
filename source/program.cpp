#include <stdio.h>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <map>

#include "overlapgraph.h"

using namespace std;

std::string path;
std::vector<std::string> reads;
std::vector<overlap> overlaps;
clock_t start_ts;
clock_t end_ts;


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

/*
Method that loads reads from FASTA file which is produced with ReadSim and 
*/
void loadReads()
{
	ifstream readsMap(path + "/reads.bnk/RED.0.map");
	ifstream readsMap2(path + "/reads.2k.10x.fasta");

	map<string, int> readMap;
	for (string line; getline(readsMap, line);)
	{
		if (line.substr(0, 3) == "RED")
			continue;

		vector<string> tokens = split(line, '\t');
		pair <string, int> p(tokens[2], stoi(tokens[0]));
		readMap.insert(p);
	}

	int no = 0;
	string read;
	for (string line; getline(readsMap2, line);)
	{

		if (line.substr(0, 1) == ">"){
			no = readMap[line.substr(1)];
			continue;
		}
		else
		{
			reads.push_back(line);
		}
	}
}

/*
Method loads overlaps from file overlaps.afg, which is produced by Minimus. Target file is in afg format
and collection of overlaps is populated.
*/
void loadoverlaps()
{
	std::string overlapsfile = path + "/overlaps.afg";

	std::ifstream overlapinput(overlapsfile.c_str());
	overlap *o = NULL;
	string line;
	for (string line; getline(overlapinput, line);)
	{
		if (line[0] == '{')
			o = (overlap *)malloc(sizeof(overlap));
		else if (line.substr(0, 3) == "adj")
			o->orientation = line[4];
		else if (line.substr(0, 3) == "rds")
		{
			string s = line.substr(4);
			vector<string> reads = split(s, ',');
			o->read1 = stoi(reads[0]) - 1;
			o->read2 = stoi(reads[1]) - 1;
		}
		else if (line.substr(0, 3) == "ahg")
		{
			o->ahg = stoi(line.substr(4));
		}
		else if (line.substr(0, 3) == "bhg")
		{
			o->bhg = stoi(line.substr(4));
		}
		else if (line[0] == '}')
		{
			overlaps.push_back(*o);
		}
	}
}

/*
Entry  point of program. Program recieves as argument path to folder where output from Minimus is located.
*/
int main(int argc, char *argv[])
{
	path = argv[1];

	std::cout << "Loading from files started ... " << endl;

	loadReads();
	loadoverlaps();

	std::cout << "Loading from files finished ... " << endl;

	OverlapGraph overlapGraph(reads, overlaps);
	std::cout << "Number of reads: " << reads.size()  << endl;

	std::cout << "Unitigging started ... " << endl;
	start_ts = clock();
	overlapGraph.runUnitigging();
	end_ts = clock();

	cout << (double)(end_ts - start_ts) / CLOCKS_PER_SEC << endl;

	std::cout << "Unitigging finished ... " << endl;

	std::cout << "Dumping to files ... " << endl;
	overlapGraph.printLayouts();
	overlapGraph.unitigsPrinting();
}

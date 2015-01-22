#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;



struct readposition
{
	int clrStart;
	int clrEnd;
	int off;
	int src;
};

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

void loadlayout(string path, vector<readposition> *positions)
{
	std::ifstream overlapinput(path.c_str());
	readposition *rp = NULL;
	string line;

	for (string line; getline(overlapinput, line);)
	{
		if (line == "{TLE")
			rp = (readposition *)malloc(sizeof(readposition));
		else if (line.substr(0, 3) == "clr")
		{
			vector<string> interval = split(line.substr(4), ',');
			rp->clrStart = stoi(interval[0]);
			rp->clrEnd = stoi(interval[1]);
		}
		else if (line.substr(0, 3) == "off")
			rp->off = stoi(line.substr(4));
		else if (line.substr(0, 3) == "src")
			rp->src = stoi(line.substr(4));
		else if (line == "}")
			positions->push_back(*rp);
	}
}

vector<readposition> correctlayout;
vector<readposition> producedlayout;

int main(int argc, char *argv[])
{
	string minimusOutput = argv[1];
	string testingOutput = argv[2];

	loadlayout(minimusOutput, &correctlayout);
	loadlayout(testingOutput, &producedlayout);

	bool isCorrect = true;
	if (correctlayout.size() != producedlayout.size())
	{
		cout << "Output is not correct" << endl;
		return 0;
	}

	for (int i = 0; i < correctlayout.size(); i++)
	{
		if (producedlayout[i].src == correctlayout[i].src)
		{
			readposition r1 = producedlayout[i];
			readposition r2 = correctlayout[i];

			if (r1.clrEnd != r2.clrEnd || r1.clrStart != r2.clrStart || r1.off != r2.off){
				isCorrect = false;
				cout << "error in read " << r1.src << endl;
			}
		}
		else
			isCorrect = false;
	}

	if (isCorrect)
		cout << "Output is correct" << endl;
	else
		cout << "Output is not correct!!!" << endl;
}

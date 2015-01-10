#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "overlapgraph.h"

#define VERTEX_INACTIVE 0
#define VERTEX_ACTIVE 1
#define VERTEX_ELIMINATED -1
#define FUZZ 10

using namespace std;

OverlapGraph::OverlapGraph(std::vector<std::string> readvector, std::vector<overlap> overlapvector)
{
	reads = readvector;
	overlaps = overlapvector;
	initialize();
}

void OverlapGraph::initialize()
{
	graph.resize(reads.size());
	for(int i=0;i<overlaps.size();i++)
	{
		overlap o = overlaps[i];

		if(o.ahg > 0 & o.bhg > 0)
			graph[o.read1].push_back(o);
		else if(o.ahg < 0 & o.bhg < 0)
			graph[o.read2].push_back(o);
	}
}

void OverlapGraph::runUnitigging()
{
	int N = reads.size();
	int *nonContainedReads = (int *)malloc(N*sizeof(int));
	memset(nonContainedReads, 0, N*sizeof(nonContainedReads));

	//removing contained reads
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < graph[i].size(); j++)
		{
			overlap o = graph[i][j];
			nonContainedReads[o.read1] = 1;
			nonContainedReads[o.read2] = 1;
		}
	}

	int *vertexstatus = (int *)malloc(N*sizeof(int));
	memset(vertexstatus, VERTEX_INACTIVE, N*sizeof(int));
	vector< vector<bool> > reduceflags(N);

	for (int i = 0; i < N; i++)
		for (int j = 0; j < graph[i].size(); j++)
			reduceflags[i].push_back(false);

	for (int I = 0; I < N; I++)
	{
		if (graph[I].size() == 0)
			continue;

		for (int j = 0; j < graph[I].size(); j++)
			vertexstatus[graph[I][j].read2] = VERTEX_ACTIVE;

		int longest = 0;
		for (int j = 0; j < graph[I].size(); j++)
		{
			if (graph[I][j].bhg > longest)
				longest = graph[I][j].bhg;
		}

		longest += FUZZ;
		cout << longest << endl;

		for (int j = 0; j < graph[I].size(); j++)
		{
			if (vertexstatus[graph[I][j].read2] == VERTEX_ACTIVE)
			{
				int x = graph[I][j].read2;
				for (int k = 0; k < graph[x].size(); k++)
				{
					if (graph[I][j].bhg + graph[x][k].bhg <= longest)
					{
						if (vertexstatus[graph[x][k].read2] == VERTEX_ACTIVE)
							vertexstatus[graph[x][k].read2] = VERTEX_ELIMINATED;
					}
				}
			}
		}

		for (int j = 0; j < graph[I].size(); j++)
		{
			if (vertexstatus[graph[I][j].read2] == VERTEX_ELIMINATED)
			{
				reduceflags[I][j] = true;
				//cout << "removed: " << graph[I][j].read1 + 1 << "->" << graph[I][j].read2 + 1 << endl;
			}
			vertexstatus[graph[I][j].read2] == VERTEX_INACTIVE;
		}
	}

	for(int i=0;i<graph.size();i++)
    {
        for(int j=0;j<graph[i].size();j++)
        {
			if (!reduceflags[i][j] & nonContainedReads[graph[i][j].read1])
                cout << graph[i][j].read1 << " -> " << graph[i][j].read2 << endl;
        }
    }

	free(vertexstatus);
	free(nonContainedReads);
	return;
}

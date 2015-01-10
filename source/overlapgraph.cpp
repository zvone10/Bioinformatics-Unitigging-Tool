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
}

void OverlapGraph::initialize()
{
	graph.resize(reads.size());
	reducedGraph.resize(reads.size());
	nonContainedReads.resize(reads.size());
	nonContainedReads.assign(reads.size(), 1);

	for(int i=0;i<overlaps.size();i++)
	{
		overlap o = overlaps[i];


		//removing contained reads
		if(o.ahg > 0 & o.bhg > 0)
			graph[o.read1].push_back(o);
		else if (o.ahg < 0 & o.bhg < 0)
			graph[o.read2].push_back(o);
		else if (o.ahg > 0 & o.bhg < 0)
			nonContainedReads[o.read2] = 0;
		else if (o.ahg < 0 & o.bhg>0)
			nonContainedReads[o.read1] = 0;
	}
}

void OverlapGraph::runUnitigging()
{
	int N = reads.size();

	
	initialize();

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
			int x = graph[I][j].read2;
			for (int k = 0; k < graph[x].size(); k++)
			{
				if (graph[x][k].bhg < FUZZ)
				{
					if (vertexstatus[graph[x][k].read2] == VERTEX_ACTIVE)
						vertexstatus[graph[x][k].read2] = VERTEX_ELIMINATED;
				}
			}
		}

		for (int j = 0; j < graph[I].size(); j++)
		{
			if (vertexstatus[graph[I][j].read2] == VERTEX_ELIMINATED)
			{
				reduceflags[I][j] = true;
			}
			vertexstatus[graph[I][j].read2] == VERTEX_INACTIVE;
		}
	}

	for(int i=0;i<graph.size();i++)
    {
        for(int j=0;j<graph[i].size();j++)
        {
			if (!reduceflags[i][j] & nonContainedReads[graph[i][j].read1] & nonContainedReads[graph[i][j].read2]){
				// Add edge to the reduced graph
				reducedGraph[i].push_back(graph[i][j]);
			}
        }
    }
	uniqueJoinCollapsing();

	free(vertexstatus);
	return;
}


void OverlapGraph::uniqueJoinCollapsing(){

	reducedGraphVertexReferenceCounter.assign(reducedGraph.size(), 0);

	// Calculate references to vertexes
	for (int i = 0; i < reducedGraph.size(); i++){
		for (int j = 0; j < reducedGraph[i].size(); j++){
			reducedGraphVertexReferenceCounter[reducedGraph[i][j].read2]++;
		}
	}

	// Find first "juncture" vertex
	int startingVertex = 0;
	bool startingVertexFound = false;
	while (startingVertex < reducedGraph.size()){
		if (reducedGraph[startingVertex].size() > 0){
			startingVertexFound = true;
			break;
		}
		startingVertex++;
	}
	if (startingVertexFound){
		calculateChunks(startingVertex);
	}
}

void OverlapGraph::calculateChunks(int chunkStartingVertexId){
	
	int currentVertexId = chunkStartingVertexId;

	if (vertexAlreadyInChunk(currentVertexId)){
		return;
	}

	// create new chunk
	Chunk chunk;

	while (true){
		if (reducedGraphVertexReferenceCounter[currentVertexId] > 1 && currentVertexId != chunkStartingVertexId){
			// special case: juncture vertex (not part of our chunk)
			calculateChunks(currentVertexId);
			break;
		}

		chunk.members.push_back(currentVertexId);

		if (reducedGraph[currentVertexId].size() < 1){
			// dead end
			break;
		}

		if (reducedGraph[currentVertexId].size() > 1){
			// juncture point - each path vertex represents a new chunk start
			for (int i = 0; i < reducedGraph[currentVertexId].size(); i++){
				calculateChunks(reducedGraph[currentVertexId][i].read2);
			}
			break;
		}

		currentVertexId = reducedGraph[currentVertexId][0].read2;
	}
	collapsedReducedGraph.push_back(chunk);

	return;
}


bool OverlapGraph::vertexAlreadyInChunk(int vertexId){
	for (int i = 0; i < collapsedReducedGraph.size(); ++i){
		for (int j = 0; j < collapsedReducedGraph[i].members.size(); ++j){
			if (vertexId == collapsedReducedGraph[i].members[j]){
				return true;
			}
		}
	}
	return false;
}
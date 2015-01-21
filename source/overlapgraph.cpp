#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include "overlapgraph.h"

#define VERTEX_INACTIVE 0
#define VERTEX_ACTIVE 1
#define VERTEX_ELIMINATED -1
#define FUZZ 10




/*
Constructor for OverlapGraph class.
*/
OverlapGraph::OverlapGraph(std::vector<std::string> readvector, std::vector<overlap> overlapvector)
{
	reads = readvector;
	overlaps = overlapvector;
}

void OverlapGraph::initialize()
{
	int num_of_reads = reads.size();
	graph.resize(num_of_reads);
	completeGraph.resize(num_of_reads);
	reducedGraph.resize(num_of_reads);

	nonContainedReads.resize(num_of_reads);
	nonContainedReads.assign(num_of_reads, 1);

	overlap_flags = (bool *)malloc(num_of_reads*sizeof(bool));
	memset(overlap_flags, 0, num_of_reads*sizeof(bool));

	for(int i=0;i<overlaps.size();i++)
	{
		overlap o = overlaps[i];

		//removing contained reads
		if(o.ahg > 0 & o.bhg > 0)
			graph[o.read1].push_back(o);
		else if (o.ahg < 0 & o.bhg < 0)
			graph[o.read2].push_back(o);
		else if (o.ahg >= 0 & o.bhg <= 0)
			nonContainedReads[o.read2] = 0;
		else if (o.ahg <= 0 & o.bhg>=0)
			nonContainedReads[o.read1] = 0;

		completeGraph[o.read1].push_back(o);
		overlap_flags[o.read1] = true;
		overlap_flags[o.read2] = true;
	}

	int *indegree = (int *)malloc(num_of_reads*sizeof(int));
	memset(indegree, 0, num_of_reads*sizeof(int));
	for (int i = 0; i < graph.size(); i++)
	{
		for (int j = 0; j < graph[i].size(); j++)
		{
			indegree[graph[i][j].read2]++;
		}
	}

	for (int i = 0; i < reads.size(); i++)
	{
		if (indegree[i] == 0 & graph[i].size() == 0)
			nonContainedReads[i] = 0;
	}

	free(indegree);
}

//////////////////////////////////////////////////////////////////////////////

void OverlapGraph::runUnitigging()
{
	using namespace std;
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
		//cout << longest << endl;

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

	cout << "Vertices and edges in unitig:" << endl;
	for(int i=0;i<graph.size();i++)
    {
        for(int j=0;j<graph[i].size();j++)
        {
			if (!reduceflags[i][j] & nonContainedReads[graph[i][j].read1] & nonContainedReads[graph[i][j].read2]){
				// Add edge to the reduced graph
				reducedGraph[i].push_back(graph[i][j]);
				cout << graph[i][j].read1 + 1<< " -> " << graph[i][j].read2 + 1<< endl;
			}
        }
    }

	uniqueJoinCollapsing();

	free(vertexstatus);
	return;
}



void OverlapGraph::unitigsPrinting()
{
	using namespace std;

	ofstream output("unitigs.afg");
	
	int counter = 1;

	output << "{unitigs" << endl;

	for (Chunk c : collapsedReducedGraph)
	{
		output << "{" << "unitig" << endl;
		output << "id: " << counter << endl;
		output << "members: ";
		for (int i = 0; i < c.members.size();i++)
		{
			if (i < c.members.size() - 1)
				output << c.members[i] + 1 << ",";
			else
				output << c.members[i] + 1 << endl;
		}
		output << "}" << endl;
		counter++;
	}

	output << "}" << endl;
	output.close();
}

void OverlapGraph::printLayouts()
{
	using namespace std;

	ofstream output("layouts.afg");

	int *offsets = (int *)malloc(reads.size()*sizeof(int));
	memset(offsets, 0, reads.size()*sizeof(int));

	calculateOffsets(offsets);

	output << "{LAY" << endl;
	for (int i = 0; i < reads.size(); i++)
	{
		if (overlap_flags[i]){
			output << "{TLE" << endl;
			output << "clr:0," << reads[i].length() << endl;
			output << "off:" << offsets[i] << endl;
			output << "src:" << i + 1 << endl;
			output << "}" << endl;
		}
	}

	output << "}" << endl;
	output.close();

	free(offsets);
}


void OverlapGraph::calculateOffsets(int *offsets)
{
	//array that stores id of non-containted read to which 
	//contained read is mapped. mapping[19] = 14 means that 
	//contained read 19 is mapped to non-contained read 14.
	int *mapping = (int *)malloc(reads.size()*sizeof(int));

	//smallest right overhang for contained read
	int *bhglen = (int *)malloc(reads.size()*sizeof(int));
	memset(mapping, -1, reads.size() * sizeof(int));

	for (int i = 0; i < reads.size(); i++)
	{
		for (int j = 0; j < completeGraph[i].size(); j++)
		{
			overlap o = completeGraph[i][j];
			if (nonContainedReads[o.read1] & nonContainedReads[o.read2]) {
				//calculating offset for non-contained reads that are part of final contig.
				offsets[o.read2] = offsets[o.read1] + o.ahg;
			}
			else if (nonContainedReads[o.read1] & !nonContainedReads[o.read2])
			{
				//mapping contained read to smallest right overhang
				if ((mapping[o.read2] == -1 || abs(o.bhg) < bhglen[o.read2]) && o.bhg <= 0)
				{
					mapping[o.read2] = o.read1;
					bhglen[o.read2] = abs(o.bhg);
				}
			}
		}
	}

	for (int i = 0; i < reads.size(); i++)
	{
		if (!nonContainedReads[i] && overlap_flags[i])
		{
			int m = mapping[i];
			overlap *o = find_overlap(mapping[i], i);
			offsets[i] = offsets[m] + o->ahg;
		}
	}

	free(mapping);
	free(bhglen);
}

overlap* OverlapGraph::find_overlap(int read1, int read2)
{
	for (int j = 0; j < completeGraph[read1].size(); j++)
		if (completeGraph[read1][j].read2 == read2)
			return &completeGraph[read1][j];
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
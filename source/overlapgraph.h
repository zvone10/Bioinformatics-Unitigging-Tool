#include<vector>
#include<string>

struct overlap
{
	int read1;
	int read2;
	char orientation;
	int ahg;
	int bhg;
};

class Chunk{
	public: 
		std::vector< int > members;
};

class OverlapGraph
{
	private:
		std::vector<std::string> reads;
		std::vector<overlap> overlaps;
		std::vector<int> nonContainedReads;
		std::vector< std::vector<overlap> > completeGraph;
		std::vector< std::vector<overlap> > graph;
		std::vector< std::vector<overlap> > reducedGraph;
		std::vector< Chunk > collapsedReducedGraph;
		std::vector<int> reducedGraphVertexReferenceCounter;
		bool *overlapFlags;
		void initialize();
		void uniqueJoinCollapsing();
		void calculateChunks(int chunkStartingVertexId);
		bool vertexAlreadyInChunk(int vertexId);
		void calculateOffsets(int *offsets);
		overlap* find_overlap(int read1, int read2);
	public:
		OverlapGraph(std::vector<std::string> reads, std::vector<overlap> overlaps);
		void runUnitigging();
		void printLayouts();
		void unitigsPrinting();
};

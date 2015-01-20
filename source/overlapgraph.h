#include<vector>
#include<string>

/**
Struct that represents overlap of two reads.

read1 -> id of first read
read2 -> id of second read
orientation -> orientation of reads. Two posibilities: 'N' for normal, and 'I' for innie orientation)
ahg -> length ofleft overhang
bhg -> length of right overhang

Example:

read1:==============================>   <- bhg ->
       <-ahg->       read2:===========================>

*/
struct overlap
{
	int read1;
	int read2;
	char orientation;
	int ahg;
	int bhg;
};

/**
Class stores all unitigs that are result of overlap graph reduction algorithm.
*/
class Chunk{
	public: 
		std::vector< int > members;
};


/**
Class provides functionality of storing overlap graph and running unitigging algorithm.
Class also prints layout to layouts.afg file and unitigs that are result of unitigging.

Example of usage:

OverlapGraph graph(std::vector<std::string> reads, std::vector<overlap> overlaps);
graph.runUnitigging(); 
graph.printLayouts(); //this line prints layout to layouts.afg.
graph.unitigsPrinting();

*/
class OverlapGraph
{
	public:
		/**
		Constructor for class.
		*/
		OverlapGraph(std::vector<std::string> reads, std::vector<overlap> overlaps);

		/**
		Method in which unitigging algorithm is implemented. This method will print on standard output all edges
		that are result of algorithm.
		*/
		void runUnitigging();

		/**
		Method prints layout of contained and non-contained reads to layouts.afg file.
		*/
		void printLayouts();

		/**
		Method creates afg file with unitig(s).
		*/
		void unitigsPrinting();
	private:

		/**
		Method which initializes and prepares data structures for unitigging algorithm. This method is called from constructor.
		In this method is performed removal of contained reads. 
		*/
		void initialize();

		/**
		Method creates union of vertices and edges into one big contig after reduction algorithm.
		*/
		void uniqueJoinCollapsing();

		/**
		Helper method for contig calculation.
		*/
		void calculateChunks(int chunkStartingVertexId);

		/**
		Method that returns true if vertex is already in chunk, false otherwise.
		*/
		bool vertexAlreadyInChunk(int vertexId);

		/**
		Method that calculates offset of vertices in final layout. This method creates same output as if Minimus assembler
		was used on same dataset.
		*/
		void calculateOffsets(int *offsets);

		/**
		Method finds overlap with read1 and read2 in overlap graph.
		*/
		overlap* find_overlap(int read1, int read2);

		/**
		Collection of reads.
		*/
		std::vector<std::string> reads;
		std::vector<overlap> overlaps;
		std::vector<int> nonContainedReads;
		std::vector< std::vector<overlap> > completeGraph;
		std::vector< std::vector<overlap> > graph;
		std::vector< std::vector<overlap> > reducedGraph;
		std::vector< Chunk > collapsedReducedGraph;
		std::vector<int> reducedGraphVertexReferenceCounter;

		/**
		Array that is dynamically allocated and stores overlap flag for each read. This flag
		determines whether read is overlap with some other read or not. Flag is true if read is overlaped with another read, 
		false otherwise.
		*/
		bool *overlap_flags;
};

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

class OverlapGraph
{
	private:
		std::vector<std::string> reads;
		std::vector<overlap> overlaps;
		std::vector< std::vector<overlap> > graph;
		void initialize();
		void containedReadRemoval();
		void transitiveEdgeRemoval();
		void uniqueJoinCollapsing();
	public:
		OverlapGraph(std::vector<std::string> reads, std::vector<overlap> overlaps);
		void runUnitigging();
};

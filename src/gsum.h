#pragma once

#include "parallel_hashmap/btree.h"
#include "parallel_hashmap/phmap.h"

#include "graph.h"

/// all algorithms below return a vector sid
/// sid[v] := the id of the super-node that contains v

/// our proposed algorithms
VI mags(Graph &g, int hash_num = 40, int round = 50);
VI mags_dm(Graph &g, int hash_num = 40, int round = 50);
/// competitor algorithms
VI randonmized(Graph &g);
VI greedy(Graph &g);
/// end algorithms

struct SuperNode {
	/// @brief the vertices in this super-node
	VI vertices;
	/// @brief the super-node neighbors of this super-node
	VI snbr;
};

struct Correction {
	/// @brief the positive correction
	VI Cp;
	/// @brief the negative correction
	VI Cm;
};

class SuperGraph {
public:
/// @brief build the optimal encoding based on sid
	SuperGraph(const VI &sid, Graph &g);
	SuperGraph() {}
	~SuperGraph();

	/// @return the degree of v
	int deg(int v);
	/// @return the neighbor set of v
	VI nbr(int v);
	/// @return the encoding cost of the super-graph
	LL get_cost();

protected:
	void __build(Graph &g);

	VI sid_;
	std::vector<SuperNode*> P;
	std::vector<Correction*> C;
};
#pragma once

#include "parallel_hashmap/btree.h"
#include "parallel_hashmap/phmap.h"
#include <omp.h>

#include "graph.h"
#include "gsum.h"


class ParaSuperGraph : public SuperGraph {
public:
	ParaSuperGraph(const VI &sid, Graph &g);

private:
	void __build(Graph &g);
};

VI parallel_mags_dm(Graph &g, int hash_num = 40, int round = 50);
VI parallel_mags(Graph &g, int hash_num = 40, int round = 50);

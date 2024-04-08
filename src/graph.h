#pragma once

#include "global.h"

struct Vertex {
	int degree;
	VI nbr;
};

class Graph {
public:
	Graph();
	~Graph();

	/// @brief sort all adjacency lists by id, and set is_nbr_ordered_ to true
	void order_nbrs();

	/// @brief adjacency lists
	std::vector<Vertex*> vertices_;
	
	/// @brief graph information
	int vertices_num_;
	LL edges_num_;
	int max_deg_;
	bool is_nbr_ordered_;
};

/// @brief create graph from a set of edges
Graph* from_edges(std::vector<PII> &edges);

/// @brief create graph from a text file, each line is an edge separated by space
Graph* from_text(std::string input_text_path);

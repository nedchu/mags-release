#include "graph.h"

#include <fstream>
#include <cstdio>
#include <algorithm>
#include <cstring>

Graph::Graph()
{
	is_nbr_ordered_ = false;
}

Graph::~Graph()
{
	if (!vertices_.empty()) {
		for (auto u : vertices_)
			delete u;
	}
}

void Graph::order_nbrs()
{
	for (int i = 0; i < vertices_.size(); i++) {
		VI &nbr = vertices_[i]->nbr;
		std::sort(nbr.begin(), nbr.end());
	}
	is_nbr_ordered_ = true;
}

Graph * from_edges(std::vector<PII>& edges)
{
	Graph *g = new Graph;
	// build n
	int n = 0;
	for (PII &e : edges) {
		n = std::max(n, std::max(e.first, e.second));
	}
	n++;

	// init graph, vn, en, vert
	g->vertices_num_ = n;
	g->edges_num_ = edges.size();
	g->vertices_ = std::vector<Vertex*>(n, nullptr);
	for (int i = 0; i < n; i++) g->vertices_[i] = new Vertex;
	
	// add edges
	for (PII &e : edges) {
		g->vertices_[e.first]->nbr.push_back(e.second);
		g->vertices_[e.second]->nbr.push_back(e.first);
	}

	// init deg and max_deg
	int max_deg = 0;
	for (int i = 0; i < n; i++) {
		g->vertices_[i]->degree = g->vertices_[i]->nbr.size();
		max_deg = std::max(max_deg, g->vertices_[i]->degree);
	}
	g->max_deg_ = max_deg;
	
	return g;
}

Graph *from_text(std::string input_text_path)
{
    std::ifstream fin(input_text_path, std::ios::in);
    if (!fin.is_open()) {
		printf("Cannot open file: %s\n", input_text_path.c_str());
		exit(1);
		return nullptr;
    }

	std::vector<PII> edges;
	int u, v;
	while (fin >> u >> v) {
		edges.push_back({u, v});
	}
	fin.close();
    return from_edges(edges);
}
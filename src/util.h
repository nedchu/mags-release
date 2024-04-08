#pragma once
#include "global.h"

/// @param x, y: sorted vectors
/// @return the number of common neighbors of x and y 
int num_common_nbr(const VI &x, const VI &y);

/// @param x, y: sorted vectors
/// @return the Jaccard similarity of x and y
double Jaccard(const VI &x, const VI &y);

/// @return the time duration between s and t (in seconds)
double get_duration(hclock::time_point s, hclock::time_point t);

struct DisjointSetUnion {
public:
	DisjointSetUnion(int n);

	int find(int x);
	void init(int n);
	void joint(int u, int v);
	bool same(int u, int v);

	VI fa;
};

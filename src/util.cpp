#include "util.h"

int num_common_nbr(const VI &x, const VI &y)
{
	int ix = 0, iy = 0;
	int szx = x.size(), szy = y.size();
	int ans = 0;
	while (ix < szx && iy < szy) {
		if (x[ix] < y[iy]) {
			ix++;
		}
		else if (x[ix] == y[iy]) {
			ans++, ix++, iy++;
		}
		else {
			iy++;
		}
	}
	return ans;
}

double Jaccard(const VI &x, const VI &y)
{
	int inter = num_common_nbr(x, y);
	return double(inter) / (x.size() + y.size() - inter);
}

double get_duration(hclock::time_point s, hclock::time_point t)
{
	return std::chrono::duration_cast<std::chrono::duration<double>>(t - s).count();
}

DisjointSetUnion::DisjointSetUnion(int n) : fa(n + 1, 0)
{
	for (int i = 0; i <= n; i++) fa[i] = i;
}

int DisjointSetUnion::find(int x)
{
	return (x == fa[x]) ? x : (fa[x] = find(fa[x]));
}

void DisjointSetUnion::init(int n)
{
	fa = VI(n + 1, 0);
	for (int i = 0; i <= n; i++) fa[i] = i;
}

void DisjointSetUnion::joint(int u, int v)
{
	u = find(u), v = find(v);
	if (u != v) fa[v] = u;
}

bool DisjointSetUnion::same(int u, int v)
{
	return find(u) == find(v);
}

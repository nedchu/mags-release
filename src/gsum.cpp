#include <random>
#include <set>
#include <map>

#include "gsum.h"
#include "util.h"


template <typename T>
int __get_cost(T const& u_A, int u, VI &sz)
{
	int cu = 0;
	for (auto pr : u_A) {
		if (pr.first == u) {
			LL pi = (LL)sz[u] * (sz[u] - 1) / 2;
			cu += std::min(pr.second / 2LL, 1 + pi - pr.second / 2);
		}
		else {
			LL pi = (LL)sz[u] * sz[pr.first];
			cu += std::min((LL)pr.second, 1 + pi - pr.second);
		}
	}
	return cu;
}

template <typename T>
int __get_merge_cost(T const& u_A, int u, int v, VI &sz)
{
	int cw = 0;
	LL szw = sz[u] + sz[v];
	for (auto pr : u_A) {
		if (pr.first == u) {
			LL pi = szw * (szw - 1) / 2;
			cw += std::min(pr.second / 2LL, 1 + pi - pr.second / 2);
		}
		else {
			LL pi = szw * sz[pr.first];
			cw += std::min((LL)pr.second, 1 + pi - pr.second);
		}
	}
	return cw;
}

VI randonmized(Graph & g)
{
	int n = g.vertices_num_;
	// random number
	std::mt19937 rng(2333);
	std::uniform_int_distribution<std::mt19937::result_type> dist;
	
	DisjointSetUnion dsu(n);
	VI cand(n, 0);
	VI nxt(n, -1);
	VI end(n, 0);
	VI sz(n, 1);
	for (int i = 0; i < n; i++) {
		cand[i] = end[i] = i;
	}

	while (!cand.empty()) {
		int pos = dist(rng) % cand.size();
		int u = cand[pos];
		u = dsu.find(u);

		int best_cw = -1e9, best_v = -1;
		double best_suv = -1e9;
		// find <= two hop neighbors
		phmap::flat_hash_set<int> onehop;
		for (int i = u; i != -1; i = nxt[i]) {
			for (auto nbr : g.vertices_[i]->nbr) {
				onehop.insert(nbr);
			}
		}
		for (int i = u; i != -1; i = nxt[i]) onehop.erase(i);

		phmap::flat_hash_set<int> sn_set;
		for (int i : onehop) {
			sn_set.insert(dsu.find(i));
		}
		for (auto cur : onehop) {
			for (auto nbr : g.vertices_[cur]->nbr) {
				sn_set.insert(dsu.find(nbr));
			}
		}
		sn_set.erase(u);

		// try merge <= two hop neighbors
		phmap::flat_hash_map<int, int> u_A;
		for (int i = u; i != -1; i = nxt[i]) {
			for (auto nbr : g.vertices_[i]->nbr) {
				u_A[dsu.find(nbr)]++;
			}
		}
		int cu = __get_cost(u_A, u, sz);

		for (auto v : sn_set) {
			phmap::flat_hash_map<int, int> v_A;
			for (int i = v; i != -1; i = nxt[i]) {
				for (auto nbr : g.vertices_[i]->nbr) {
					v_A[dsu.find(nbr)]++;
				}
			}
			int cv = __get_cost(v_A, v, sz);

			// build cw
			for (auto pr : u_A) {
				v_A[pr.first] += pr.second;
			}
			if (v_A.find(u) != v_A.end()) {
				v_A[v] += v_A[u];
				v_A.erase(u);
			}
			int cw = __get_merge_cost(v_A, v, u, sz);
			double new_suv = double(cu + cv - cw) / (cu + cv);
			if (new_suv > best_suv) {
				best_v = v;
				best_cw = cw;
				best_suv = new_suv;
			}
		}

		// merge u v if suv > 0, otherwise discard u
		if (best_suv > 0) {
			dsu.joint(u, best_v);
			//c[u] = best_cw;
			nxt[end[u]] = best_v;
			end[u] = end[best_v];
			sz[u] += sz[best_v];
		}
		std::swap(cand[pos], cand[cand.size() - 1]);
		cand.pop_back();
	}

	// build sid
	VI index;
	for (int i = 0; i < n; i++) {
		if (dsu.fa[i] == i) index.push_back(i);
	}
	VI ans(n, 0);
	for (int i = 0; i < n; i++) {
		int u = dsu.find(i);
		ans[i] = std::lower_bound(index.begin(), index.end(), u) - index.begin();;
	}
	return ans;
}

inline 
PII min_pair(int u, int v)
{
	return { std::min(u, v), std::max(u, v) };
}

VI greedy(Graph & g)
{
	int n = g.vertices_num_;

	DisjointSetUnion dsu(n);
	VI nxt(n, -1);
	VI end(n, 0);
	VI sz(n, 1);
	for (int i = 0; i < n; i++) {
		end[i] = i;
	}
	phmap::btree_set<PD_PII, std::greater<>> pq;
	phmap::btree_map<PII, double> W;

	for (int u = 0; u < n; u++) {
		VI &nbr = g.vertices_[u]->nbr;
		phmap::flat_hash_set<int> onehop(nbr.begin(), nbr.end());
		phmap::flat_hash_set<int> twohop;
		for (auto v : nbr) {
			twohop.insert(v);
			for (auto w : g.vertices_[v]->nbr) {
				twohop.insert(w);
			}
		}
		twohop.erase(u);

		int cu = nbr.size();
		for (auto v : twohop) {
			if (u > v) continue;
			int cv = g.vertices_[v]->degree;
			int common_neighbor = 0;
			for (auto w : g.vertices_[v]->nbr) {
				if (onehop.contains(w)) common_neighbor++;
			}
			int cw = onehop.contains(v) ? cu + cv - common_neighbor - 1 : cu + cv - common_neighbor;
			double suv = (cu + cv - cw) / double(cu + cv);
			if (suv > 0) {
				pq.insert({ suv, {u, v} });
				W[{u, v}] = suv;
				W[{v, u}] = suv;
			}
		}
	}

	while (!pq.empty()) {
		double score = pq.begin()->first;
		int merge_u, merge_v;
		std::tie(merge_u, merge_v) = pq.begin()->second;

		// find candidates to change
		phmap::flat_hash_set<int> cand;
		cand.insert(merge_u);
		cand.insert(merge_v);
		for (int i = merge_u; i != -1; i = nxt[i]) {
			for (int w : g.vertices_[i]->nbr) {
				cand.insert(dsu.find(w));
			}
		}
		for (int i = merge_v; i != -1; i = nxt[i]) {
			for (int w : g.vertices_[i]->nbr) {
				cand.insert(dsu.find(w));
			}
		}
		
		// remove old suv
		std::vector<PII> to_remove;
		for (int w : cand) {
			for (auto it = W.lower_bound({ w, 0 }); it != W.end() && it->first.first == w; it++) {
				to_remove.push_back(it->first);
			}
		}
		for (auto e : to_remove) {
			auto it = W.find(e);
			if (it != W.end()) {
				pq.erase(pq.find({ it->second, min_pair(e.first, e.second) }));
				W.erase(it);
				W.erase(W.find({ e.second, e.first }));
			}
		}

		// merge u v
		dsu.joint(merge_u, merge_v);
		nxt[end[merge_u]] = merge_v;
		end[merge_u] = end[merge_v];
		sz[merge_u] += sz[merge_v];

		cand.erase(merge_v);
		for (int u : cand) {
			// find <= two hop neighbors
			phmap::flat_hash_set<int> onehop;
			for (int i = u; i != -1; i = nxt[i]) {
				for (auto nbr : g.vertices_[i]->nbr) {
					onehop.insert(nbr);
				}
			}
			for (int i = u; i != -1; i = nxt[i]) onehop.erase(i);

			phmap::flat_hash_set<int> sn_set;
			for (int i : onehop) {
				sn_set.insert(dsu.find(i));
			}
			for (auto cur : onehop) {
				for (auto nbr : g.vertices_[cur]->nbr) {
					sn_set.insert(dsu.find(nbr));
				}
			}
			sn_set.erase(u);

			// try merge <= two hop neighbors
			phmap::flat_hash_map<int, int> u_A;
			for (int i = u; i != -1; i = nxt[i]) {
				for (auto nbr : g.vertices_[i]->nbr) {
					u_A[dsu.find(nbr)]++;
				}
			}
			int cu = __get_cost(u_A, u, sz);

			for (auto sn_v : sn_set) {
				phmap::flat_hash_map<int, int> v_A;
				for (int i = sn_v; i != -1; i = nxt[i]) {
					for (auto nbr : g.vertices_[i]->nbr) {
						v_A[dsu.find(nbr)]++;
					}
				}
				int cv = __get_cost(v_A, sn_v, sz);

				// build cw
				for (auto pr : u_A) {
					v_A[pr.first] += pr.second;
				}
				if (v_A.find(u) != v_A.end()) {
					v_A[sn_v] += v_A[u];
					v_A.erase(u);
				}
				int cw = __get_merge_cost(v_A, sn_v, u, sz);

				double new_suv = double(cu + cv - cw) / (cu + cv);
				if (new_suv > 0) {
					W[{u, sn_v}] = new_suv;
					W[{sn_v, u}] = new_suv;
					pq.insert({ new_suv, min_pair(u, sn_v) });
				}
			}
		}
	}

	// build sid
	VI index;
	for (int i = 0; i < n; i++) {
		if (dsu.fa[i] == i) index.push_back(i);
	}
	VI ans(n, 0);
	for (int i = 0; i < n; i++) {
		int u = dsu.find(i);
		ans[i] = std::lower_bound(index.begin(), index.end(), u) - index.begin();;
	}
	return ans;
}

inline
PD_PII trans(int u, const PID &pid)
{
	return { pid.second, {u, pid.first} };
}

VI mags(Graph & g, int hash_num, int round)
{
	hclock::time_point s1 = hclock::now();
	g.order_nbrs();
	int n = g.vertices_num_;
	std::mt19937 rng(2333);
	int hit_limit = 30;
	int sample_num = 5;

	DisjointSetUnion dsu(n);
	VI sz(n, 1);
	std::vector<phmap::flat_hash_map<int, int>> cnt;
	phmap::btree_set<PD_PII, std::greater<>> pq_all;
	std::vector<phmap::flat_hash_map<int, double>> W;

	cnt.resize(n);
	for (int i = 0; i < n; i++) {
		for (int v : g.vertices_[i]->nbr) {
			cnt[i][v] = 1;
		}
	}

	VVI min_hash(n, VI(hash_num, -1));
	for (int hh = 0; hh < hash_num; hh++) {
		VI arr(n, 0);
		for (int j = 0; j < n; j++) arr[j] = j;
		std::shuffle(arr.begin(), arr.end(), rng);
		for (int i = 0; i < n; i++) {
			int u = arr[i];
			for (int v : g.vertices_[u]->nbr) {
				if (min_hash[v][hh] == -1) {
					min_hash[v][hh] = i;
				}
			}
		}
	}

	W.resize(n);
	for (int u = 0; u < n; u++) {
		VI &nbr = g.vertices_[u]->nbr;
		phmap::flat_hash_set<int> twohop;
		for (auto v : nbr) {
			if (v > u) twohop.insert(v);
		}

		int ccnt = 0;
		for (auto v : nbr) {
			for (auto w : g.vertices_[v]->nbr) {
				if (w > u) twohop.insert(w);
			}
			if (++ccnt == sample_num) break;
		}
		twohop.erase(u);

		phmap::btree_set<PII> pq0;
		for (auto v : twohop) {
			int hit = 0;
			for (int hh = 0; hh < hash_num; hh++) hit += min_hash[v][hh] == min_hash[u][hh];
			if (pq0.size() < hit_limit) {
				pq0.insert({ hit, v });
			}
			else if (pq0.begin()->first < hit) {
				pq0.erase(pq0.begin());
				pq0.insert({ hit, v });
			}
		}

		int cu = nbr.size();
		for (auto pii : pq0) {
			int v = pii.second;
			int cv = g.vertices_[v]->degree;
			int common_neighbor = num_common_nbr(nbr, g.vertices_[v]->nbr);
			int cw = cnt[u].contains(v) ? cu + cv - common_neighbor - 1 : cu + cv - common_neighbor;
			double suv = (cu + cv - cw) / double(cu + cv);
			pq_all.insert({ suv, min_pair(u, v) });
			W[u].insert({ v, suv });
			W[v].insert({ u, suv });
		}
	}

	hclock::time_point s2 = hclock::now();
	double ratio = std::pow(0.005, 1. / (round - 1));
	double t = 0.5 / ratio;
	for (int cur_r = 0; cur_r < round; cur_r++) {
		t *= ratio;
		if (cur_r == round - 1) t = 1e-3;
		phmap::flat_hash_set<int> to_update;

		VI to_remove;
		for (auto &pd_pii : pq_all) {
			int u, v;
			std::tie(u, v) = pd_pii.second;
			u = dsu.find(u), v = dsu.find(v);
			if (u == v) continue;
			if (pd_pii.first < t) break;

			int cu = __get_cost(cnt[u], u, sz);
			int cv = __get_cost(cnt[v], v, sz);

			phmap::flat_hash_map<int, int> u_A(cnt[u]);
			for (auto pr : cnt[v]) {
				u_A[pr.first] += pr.second;
			}
			if (u_A.find(v) != u_A.end()) {
				u_A[u] += u_A[v];
				u_A.erase(v);
			}
			int cw = __get_merge_cost(u_A, u, v, sz);
			double score = double(cu + cv - cw) / (cu + cv);

			if (score < t) continue;
			// merge u v
			dsu.joint(u, v);
			sz[u] += sz[v];

			for (auto pr : cnt[v]) {
				cnt[pr.first].erase(v);
			}
			cnt[v].clear();
			cnt[u] = u_A;
			for (auto pr : u_A) {
				if (pr.first != u) {
					cnt[pr.first][u] = pr.second;
				}
			}
			to_remove.push_back(v);
			to_update.insert(u);
		}
		for (int v : to_remove) {
			int u = dsu.find(v);
			// remove old suv
			for (auto pid : W[v]) {
				int w = pid.first;
				pq_all.erase(pq_all.find({ pid.second, min_pair(v, w) }));
				W[w].erase(v);
				if (w != u && W[u].find(w) == W[u].end()) {
					W[u][w] = W[w][u] = -1;
					pq_all.insert({ -1, min_pair(u, w) });
				}
			}
			W[v].clear();
		}

		phmap::flat_hash_set<int> to_update2(std::move(to_update));
		for (int uuu : to_update2) {
			if (dsu.find(uuu) != uuu) continue;
			to_update.insert(uuu);
			for (auto pid : cnt[uuu]) {
				to_update.insert(pid.first);
			}
		}

		phmap::flat_hash_map<int, int> cus;
		for (int uuu : to_update) {
			if (cus.find(uuu) == cus.end()) {
				cus[uuu] = __get_cost(cnt[uuu], uuu, sz);
			}
			int cu = cus[uuu];

			VI to_remove;
			std::vector<PID> to_update_suv;
			for (auto &pid : W[uuu]) {
				int vvv = pid.first;
				double old_suv = pid.second;
				if (to_update.contains(vvv) && vvv < uuu) continue;
				if (cus.find(vvv) == cus.end()) {
					cus[vvv] = __get_cost(cnt[vvv], vvv, sz);
				}
				int cv = cus[vvv];

				phmap::flat_hash_map<int, int> u_A(cnt[uuu]);
				for (auto pr : cnt[vvv]) {
					u_A[pr.first] += pr.second;
				}
				if (u_A.find(vvv) != u_A.end()) {
					u_A[uuu] += u_A[vvv];
					u_A.erase(vvv);
				}
				int cw = __get_merge_cost(u_A, uuu, vvv, sz);

				double new_suv = double(cu + cv - cw) / (cu + cv);
				if (old_suv == new_suv) continue;

				if (new_suv <= -0.03) {
					to_remove.push_back(vvv);
				}
				else {
					to_update_suv.push_back({ vvv, old_suv });
					pid.second = new_suv;
				}
			}
			for (auto vvv : to_remove) {
				auto it = pq_all.find({ W[uuu][vvv], min_pair(uuu, vvv) });
				pq_all.erase(it);
				W[uuu].erase(vvv);
				W[vvv].erase(uuu);
			}
			for (auto pid : to_update_suv) {
				int vvv = pid.first;
				W[vvv][uuu] = W[uuu][vvv];
				auto it = pq_all.find({ pid.second, min_pair(uuu, vvv) });
				pq_all.erase(it);
				pq_all.insert({ W[uuu][vvv], min_pair(uuu, vvv) });
			}
		}
	}

	// build sid
	VI index;
	for (int i = 0; i < n; i++) {
		if (dsu.fa[i] == i) index.push_back(i);
	}
	VI ans(n, 0);
	for (int i = 0; i < n; i++) {
		int u = dsu.find(i);
		ans[i] = std::lower_bound(index.begin(), index.end(), u) - index.begin();;
	}
	return ans;
}


const int max_group = 500;
const int max_cnt_down = 10;

void dfs_sep(std::vector<PII> &to_sort, VVI &min_hash, VI &h_choice, VI &ans,
	int d, int l, int r)
{
	if (d == max_cnt_down || r - l <= max_group) {
		ans.push_back(l);
		return;
	}
	int hh = h_choice[d];
	for (int i = l; i < r; i++) {
		int u = to_sort[i].second;
		to_sort[i].first = min_hash[u][hh];
	}
	std::sort(to_sort.begin() + l, to_sort.begin() + r);

	int prev = l;
	for (int i = l + 1; i < r; i++) {
		if (to_sort[i].first != to_sort[i - 1].first) {
			dfs_sep(to_sort, min_hash, h_choice, ans, d + 1, prev, i);
			prev = i;
		}
	}
	dfs_sep(to_sort, min_hash, h_choice, ans, d + 1, prev, r);
}

VI get_sep(std::vector<PII> &to_sort, VVI &min_hash, VI &h_choice)
{
	VI sep;
	dfs_sep(to_sort, min_hash, h_choice, sep, 0, 0, to_sort.size());
	sep.push_back(to_sort.size());
	return sep;
}

VI mags_dm(Graph & g, int hash_num, int round)
{
	hclock::time_point s1 = hclock::now();
	int n = g.vertices_num_;
	std::mt19937 rng(2333);
	int hit_limit = 5;

	DisjointSetUnion dsu(n);
	VI sz(n, 1);
	std::vector<phmap::flat_hash_map<int, int>> cnt;
	std::vector<phmap::flat_hash_set<int>> exclude;
	cnt.resize(n);
	exclude.resize(n);
	for (int i = 0; i < n; i++) {
		for (int v : g.vertices_[i]->nbr) {
			cnt[i][v] = 1;
		}
	}

	VVI min_hash(n, VI(hash_num, -1));
	for (int hh = 0; hh < hash_num; hh++) {
		VI arr(n, 0);
		for (int j = 0; j < n; j++) arr[j] = j;
		std::shuffle(arr.begin(), arr.end(), rng);
		for (int i = 0; i < n; i++) {
			int u = arr[i];
			for (int v : g.vertices_[u]->nbr) {
				if (min_hash[v][hh] == -1) {
					min_hash[v][hh] = i;
				}
			}
		}
	}
	hclock::time_point s2 = hclock::now();

	VI h_choice;
	for (int i = 0; i < hash_num; i++) h_choice.push_back(i);
	double ratio = std::pow(0.005, 1. / (round - 1));
	double t = 0.5 / ratio;
	for (int cur_r = 0; cur_r < round; cur_r++) {
		t *= ratio;
		int hh = cur_r % hash_num;
		std::vector<PII> to_sort;
		for (int i = 0; i < n; i++) {
			if (dsu.fa[i] == i) {
				to_sort.push_back({ 0, i });
			}
		}

		std::shuffle(h_choice.begin(), h_choice.end(), rng);
		VI sep = get_sep(to_sort, min_hash, h_choice);
		int sum = 0;
		for (int i = 0; i + 1 < sep.size(); i++) {
			int l = sep[i], r = sep[i + 1];
			sum += (r - l) * (r - l);
			VI group;
			for (int j = l; j < r; j++) {
				group.push_back(to_sort[j].second);
			}
			std::shuffle(group.begin(), group.end(), rng);
			while (group.size() > 1) {
				int u = group.back();
				group.pop_back();
				std::set<PII> pq;
				for (int j = 0; j < group.size(); j++) {
					int v = group[j];
					int hit = 0;
					for (int k = 0; k < hash_num; k++) hit += min_hash[u][k] == min_hash[v][k];
					if (pq.size() < hit_limit) {
						pq.insert({ hit, j });
					}
					else if (hit > pq.begin()->first) {
						pq.erase(pq.begin());
						pq.insert({ hit, j });
					}
				}

				int cu = __get_cost(cnt[u], u, sz);
				int v_j = -1;
				double max_score = -1e9;

				for (auto &pr : pq) {
					int v = group[pr.second];
					if (exclude[u].find(v) != exclude[u].end()) continue;
					int cv = __get_cost(cnt[v], v, sz);
					// build cw
					phmap::flat_hash_map<int, int> u_A(cnt[u]);
					for (auto pr : cnt[v]) {
						u_A[pr.first] += pr.second;
					}
					if (u_A.find(v) != u_A.end()) {
						u_A[u] += u_A[v];
						u_A.erase(v);
					}
					int cw = __get_merge_cost(u_A, u, v, sz);
					double new_suv = double(cu + cv - cw) / (cu + cv);
					if (new_suv > max_score) {
						v_j = pr.second;
						max_score = new_suv;
					}
				}

				if (v_j == -1) continue;
				int v = group[v_j];
				// try merge <= two hop neighbors
				if (max_score >= t - 1e-5) {
					dsu.joint(u, v);
					sz[u] += sz[v];
					group[v_j] = u;
					for (int k = 0; k < hash_num; k++) {
						min_hash[u][k] = std::min(min_hash[u][k], min_hash[v][k]);
					}

					for (auto pr : cnt[v]) {
						cnt[u][pr.first] += pr.second;
					}
					if (cnt[u].find(v) != cnt[u].end()) {
						cnt[u][u] += cnt[u][v];
					}
					for (auto pr : cnt[v]) {
						if (pr.first != v) cnt[pr.first].erase(v);
					}
					cnt[v].clear();
					for (auto pr : cnt[u]) {
						if (pr.first != u) {
							cnt[pr.first][u] = pr.second;
						}
					}

					exclude[u].insert(exclude[v].begin(), exclude[v].end());
					exclude[v].clear();
				} else if (max_score <= -0.03) {
					exclude[u].insert(v);
					exclude[v].insert(u);
				}
			}
		}
	}

	// build sid
	VI index;
	for (int i = 0; i < n; i++) {
		if (dsu.fa[i] == i) index.push_back(i);
	}
	VI ans(n, 0);
	for (int i = 0; i < n; i++) {
		int u = dsu.find(i);
		ans[i] = std::lower_bound(index.begin(), index.end(), u) - index.begin();;
	}
	
	return ans;
}

SuperGraph::SuperGraph(const VI & sid, Graph & g) : sid_(sid)
{
	__build(g);
}

SuperGraph::~SuperGraph()
{
	for (auto p : P) delete p;
	for (auto p : C) delete p;
}

int SuperGraph::deg(int v)
{
	int sid = sid_[v];
	int ans = C[v]->Cp.size() - C[v]->Cm.size();
	for (auto sn : P[sid]->snbr) {
		ans += P[sn]->vertices.size();
		if (sn == sid) ans--;
	}
	return ans;
}

VI SuperGraph::nbr(int v)
{
	phmap::flat_hash_set<int> se;
	se.insert(v);
	for (int v : C[v]->Cm) {
		se.insert(v);
	}
	VI ans(C[v]->Cp);
	for (int sn : P[sid_[v]]->snbr) {
		for (int u : P[sn]->vertices) {
			if (!se.contains(u)) ans.push_back(u);
		}
	}
	return ans;
}

LL SuperGraph::get_cost()
{
	LL cost = 0;
	for (auto p : P) {
		cost += p->snbr.size();
	}
	for (auto p : C) {
		cost += p->Cm.size() + p->Cp.size();
	}
	return cost;
}

void SuperGraph::__build(Graph & g)
{
	assert(g.vertices_num_ == sid_.size());
	int n = g.vertices_num_;
	if (!g.is_nbr_ordered_) g.order_nbrs();

	// initialize P and C
	for (int i = 0; i < n; i++) {
		int s = sid_[i];
		while (s >= P.size()) {
			P.push_back(new SuperNode);
		}
		P[s]->vertices.push_back(i);
	}
	C = std::vector<Correction*>(n, nullptr);
	for (int i = 0; i < n; i++) {
		C[i] = new Correction;
	}

	// build A_uv
	phmap::flat_hash_map<PII, int> ma;
	for (int i = 0; i < n; i++) {
		int su = sid_[i];
		for (int nbr : g.vertices_[i]->nbr) {
			if (nbr > i) break;
			int sv = sid_[nbr];
			if (sv < su) {
				ma[{sv, su}]++;
			}
			else {
				ma[{su, sv}]++;
			}
		}
	}

	// build edges in P
	for (auto it = ma.begin(); it != ma.end(); it++) {
		int su, sv;
		std::tie(su, sv) = it->first;
		int cnt = it->second;
		int susz = P[su]->vertices.size();
		if (su == sv) {
			if (2 * cnt > susz * (susz - 1) / 2) {
				P[su]->snbr.push_back(su);
			}
		}
		else {
			if (2 * cnt > susz * P[sv]->vertices.size() + 1) {
				P[su]->snbr.push_back(sv);
				P[sv]->snbr.push_back(su);
			}
		}
	}

	// build C, neighbor is ordered
	for (int s = 0; s < P.size(); s++) {
		VI temp;
		for (int snbr : P[s]->snbr) {
			VI &snbrV = P[snbr]->vertices;
			temp.insert(temp.end(), snbrV.begin(), snbrV.end());
		}
		std::sort(temp.begin(), temp.end());
		for (int v : P[s]->vertices) {
			Correction &cur = *C[v];
			VI &nbr = g.vertices_[v]->nbr;
			int ix = 0, iy = 0;
			int szx = temp.size(), szy = nbr.size();
			while (ix < szx && iy < szy) {
				if (temp[ix] < nbr[iy]) {
					if (temp[ix] != v) cur.Cm.push_back(temp[ix]);
					ix++;
				}
				else if (temp[ix] == nbr[iy]) {
					ix++, iy++;
				}
				else {
					cur.Cp.push_back(nbr[iy++]);
				}
			}
			while (ix < szx) {
				if (temp[ix] != v) cur.Cm.push_back(temp[ix]);
				ix++;
			}
			while (iy < szy) cur.Cp.push_back(nbr[iy++]);
		}
	}
}
#include <omp.h>
#include <map>
#include <random>
#include <set>

#if defined(__linux__)
#include <parallel/algorithm>
#endif

#include "pgsum.h"
#include "util.h"

#pragma omp declare reduction (merge_vpii : std::vector<PII> : omp_out.insert(omp_out.end(),	std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (merge_vi : std::vector<int> : omp_out.insert(omp_out.end(),	std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (merge_set_i : std::set<int> : omp_out.insert(std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))
#pragma omp declare reduction (merge_vpd_pii : std::vector<PD_PII> : omp_out.insert(omp_out.end(),	std::make_move_iterator(omp_in.begin()), std::make_move_iterator(omp_in.end())))

ParaSuperGraph::ParaSuperGraph(const VI & sid, Graph & g)
{
	sid_ = sid;
	__build(g);
}

void ParaSuperGraph::__build(Graph & g)
{
	assert(g.vertices_num_ == sid_.size());
	int n = g.vertices_num_;
	if (!g.is_nbr_ordered_) {
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			auto &nbr = g.vertices_[i]->nbr;
			std::sort(nbr.begin(), nbr.end());
		}
		g.is_nbr_ordered_ = true;
	}

	// initialize P and C
	int mx = 0;
#pragma omp parallel for reduction(max : mx)
	for (int i = 0; i < n; i++) {
		mx = std::max(mx, sid_[i]);
	}
	P.resize(mx + 1);
	C.resize(n);

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i <= mx; i++) {
			P[i] = new SuperNode;
		}
#pragma omp for
		for (int i = 0; i < n; i++) {
			C[i] = new Correction;
		}
#pragma omp single
		{
			for (int i = 0; i < n; i++) {
				int s = sid_[i];
				P[s]->vertices.push_back(i);
			}
		}
#pragma omp for
		for (int su = 0; su < P.size(); su++) {
			std::map<int, int> ma;
			for (auto u : P[su]->vertices) {
				for (int nbr : g.vertices_[u]->nbr) {
					int sv = sid_[nbr];
					ma[sv]++;
				}
			}
			for (auto &pr : ma) {
				int sv = pr.first;
				int cnt = pr.second;
				LL susz = P[su]->vertices.size();
				if (su == sv) {
					if (cnt > susz * (susz - 1) / 2) {
						P[su]->snbr.push_back(su);
					}
				}
				else {
					if (2 * cnt > susz * P[sv]->vertices.size() + 1) {
						P[su]->snbr.push_back(sv);
					}
				}
			}
		}
#pragma omp barier
#pragma omp for
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
}

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

const int max_group = 500;
const int max_cnt_down = 10;

void dfs_sep_para(std::vector<PII> &to_sort, VVI &min_hash, VI &h_choice, VI &ans,
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
#if defined(__linux__)
	__gnu_parallel::sort(to_sort.begin() + l, to_sort.begin() + r);
#else
	std::sort(to_sort.begin() + l, to_sort.begin() + r);
#endif

	int prev = l;
	for (int i = l + 1; i < r; i++) {
		if (to_sort[i].first != to_sort[i - 1].first) {
			dfs_sep_para(to_sort, min_hash, h_choice, ans, d + 1, prev, i);
			prev = i;
		}
	}
	dfs_sep_para(to_sort, min_hash, h_choice, ans, d + 1, prev, r);
}

VI get_sep_para(std::vector<PII> &to_sort, VVI &min_hash, VI &h_choice)
{
	VI sep;
	dfs_sep_para(to_sort, min_hash, h_choice, sep, 0, 0, to_sort.size());
	sep.push_back(to_sort.size());
	return sep;
}

VI parallel_mags_dm(Graph & g, int hash_num, int round)
{
	int n = g.vertices_num_;
	std::vector<std::mt19937*> rngs;
	int hit_limit = 5;
	int thread_num;
	double ratio = std::pow(0.005, 1. / (round - 1));

	DisjointSetUnion dsu(n);
	VI sz(n, 1);
	VVI min_hash(n, VI(hash_num, -1));
	std::vector<phmap::flat_hash_map<int, int>> cnt;
	std::vector<std::mutex*> locks;
	cnt.resize(n);
	locks.resize(n);

#pragma omp parallel
	{
#pragma omp single
		{
			thread_num = omp_get_num_threads();
			rngs.resize(thread_num);
			std::seed_seq seq{ time(0) };
			std::vector<std::uint32_t> seeds(thread_num);
			seq.generate(seeds.begin(), seeds.end());
			for (int i = 0; i < thread_num; i++) {
				rngs[i] = new std::mt19937(seeds[i]);
			}
		}
#pragma omp for 
		for (int i = 0; i < n; i++) {
			for (int v : g.vertices_[i]->nbr) {
				cnt[i][v] = 1;
			}
			locks[i] = new std::mutex;
		}
#pragma omp for
		for (int hh = 0; hh < hash_num; hh++) {
			VI arr(n, 0);
			for (int j = 0; j < n; j++) arr[j] = j;
			int tn = omp_get_thread_num();
			std::shuffle(arr.begin(), arr.end(), *rngs[tn]);
			for (int i = 0; i < n; i++) {
				int u = arr[i];
				for (int v : g.vertices_[u]->nbr) {
					if (min_hash[v][hh] == -1) {
						min_hash[v][hh] = i;
					}
				}
			}
		}
	}

	VI h_choice;
	for (int i = 0; i < hash_num; i++) h_choice.push_back(i);
	double t = 0.5 / ratio;
	for (int cur_r = 0; cur_r < round; cur_r++) {
		t *= ratio;
		int hh = cur_r % hash_num;
		std::vector<PII> to_sort;
		VI to_remove;

		// build to_sort
#pragma omp parallel for reduction(merge_vpii: to_sort)
		for (int i = 0; i < n; i++) {
			if (dsu.fa[i] == i) {
				to_sort.push_back({ 0, i });
			}
		}
		std::shuffle(h_choice.begin(), h_choice.end(), *rngs[0]);
		VI sep = get_sep_para(to_sort, min_hash, h_choice);
		int sep_sz = (int)sep.size() - 1;
#pragma omp parallel for schedule(dynamic) reduction(merge_vi: to_remove)
		for (int i = 0; i < sep_sz; i++) {
			int l = sep[i], r = sep[i + 1];
			VI group;
			for (int j = l; j < r; j++) {
				group.push_back(to_sort[j].second);
			}
			int tn = omp_get_thread_num();
			std::shuffle(group.begin(), group.end(), *rngs[tn]);

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

				locks[u]->lock();
				auto tmp_cntu = cnt[u];
				locks[u]->unlock();
				int cu = __get_cost(tmp_cntu, u, sz);
				int v_j = -1;
				double max_score = -1e9;

				for (auto &pr : pq) {
					int v = group[pr.second];
					// build cw
					phmap::flat_hash_map<int, int> u_A(tmp_cntu);

					locks[v]->lock();
					int cv = __get_cost(cnt[v], v, sz);
					for (auto pr : cnt[v]) {
						u_A[pr.first] += pr.second;
					}
					locks[v]->unlock();

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
					group[v_j] = u;
					to_remove.push_back(v);

#pragma omp critical(merge)
					{
						dsu.joint(u, v);
						sz[u] += sz[v];
						sz[v] = 0;
					}

					for (int k = 0; k < hash_num; k++) {
						min_hash[u][k] = std::min(min_hash[u][k], min_hash[v][k]);
					}

					locks[v]->lock();
					auto tmpcntv = cnt[v];
					for (auto &pr : cnt[v]) {
						pr.second = 0;
					}
					locks[v]->unlock();

					for (auto &pr : tmpcntv) {
						if (pr.first == v || pr.first == u) continue;
						locks[pr.first]->lock();
						cnt[pr.first][v] = 0;
						cnt[pr.first][u] += pr.second;
						locks[pr.first]->unlock();
					}

					locks[u]->lock();
					for (auto pr : tmpcntv) {
						cnt[u][pr.first] += pr.second;
					}
					if (cnt[u].find(v) != cnt[u].end()) {
						cnt[u][u] += cnt[u][v];
						cnt[u].erase(v);
					}
					locks[u]->unlock();
				}
			}
		}
		for (auto v : to_remove) {
			for (auto pr : cnt[v]) {
				cnt[pr.first].erase(v);
			}
			cnt[v].clear();
		}
	}

	for (int i = 0; i < n; i++) {
		delete locks[i];
	}
	for (int i = 0; i < thread_num; i++) {
		delete rngs[i];
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


inline
PII min_pair(int u, int v)
{
	return { std::min(u, v), std::max(u, v) };
}

VI parallel_mags(Graph & g, int hash_num, int round)
{
	hclock::time_point s1 = hclock::now();
	int n = g.vertices_num_;
	std::vector<std::mt19937*> rngs;
	int avg_deg = std::round((double)g.edges_num_ / g.vertices_num_);
	int hit_limit = std::min(30, avg_deg * 10);
	double ratio = std::pow(0.005, 1. / (round - 1));
	int thread_num;

	DisjointSetUnion dsu(n);
	VI sz(n, 1);
	std::vector<phmap::flat_hash_map<int, int>> cnt;
	std::set<PD_PII, std::greater<>> pq_all;
	std::vector<phmap::flat_hash_map<int, double>> W;

	VVI min_hash(n, VI(hash_num, -1));
	cnt.resize(n);
	W.resize(n);

	std::vector<std::mutex*> locks;
	locks.resize(n);

#pragma omp parallel
	{
#pragma omp single
		{
			thread_num = omp_get_num_threads();
			rngs.resize(thread_num);
			std::seed_seq seq{ time(0) };
			std::vector<std::uint32_t> seeds(thread_num);
			seq.generate(seeds.begin(), seeds.end());
			for (int i = 0; i < thread_num; i++) {
				rngs[i] = new std::mt19937(seeds[i]);
			}
		}
#pragma omp for 
		for (int i = 0; i < n; i++) {
			auto &nbr = g.vertices_[i]->nbr;
			std::sort(nbr.begin(), nbr.end());
			for (int v : nbr) {
				cnt[i][v] = 1;
			}
			locks[i] = new std::mutex;
		}
		// build min hash
#pragma omp for
		for (int hh = 0; hh < hash_num; hh++) {
			VI arr(n, 0);
			for (int j = 0; j < n; j++) arr[j] = j;
			int tn = omp_get_thread_num();
			std::shuffle(arr.begin(), arr.end(), *rngs[tn]);
			for (int i = 0; i < n; i++) {
				int u = arr[i];
				for (int v : g.vertices_[u]->nbr) {
					if (min_hash[v][hh] == -1) {
						min_hash[v][hh] = i;
					}
				}
			}
		}
		std::vector<PD_PII> pq_all_pri;
#pragma omp for schedule(dynamic) nowait
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
				if (++ccnt == 5) break;
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
			std::vector<PID> to_insert;
			for (auto pii : pq0) {
				int v = pii.second;
				int cv = g.vertices_[v]->degree;
				int common_neighbor = num_common_nbr(nbr, g.vertices_[v]->nbr);
				int cw = cnt[u].contains(v) ? cu + cv - common_neighbor - 1 : cu + cv - common_neighbor;
				double suv = (cu + cv - cw) / double(cu + cv);
				to_insert.push_back({ v, suv });
				pq_all_pri.push_back({ suv, min_pair(u, v) });
			}
			locks[u]->lock();
			for (auto pid : to_insert) {
				W[u][pid.first] = pid.second;
			}
			locks[u]->unlock();

			for (auto pid : to_insert) {
				int v = pid.first;
				locks[v]->lock();
				W[v][u] = pid.second;
				locks[v]->unlock();
			}
		}
#pragma omp critical
		{
			for (auto &pd_pii : pq_all_pri) {
				pq_all.insert(pd_pii);
			}
		}
#pragma omp barrier
	}
	min_hash.clear();

	double t = 0.5 / ratio;
	DisjointSetUnion dsu_merge(n);
	for (int cur_r = 0; cur_r < round; cur_r++) {
		t *= ratio;
		if (cur_r == round - 1) t = 1e-3;

		int cc = 0;
		for (auto &pd_pii : pq_all) {
			if (pd_pii.first < t) break;
			PII uv = pd_pii.second;
			dsu_merge.joint(uv.first, uv.second);
			cc++;
		}
		phmap::flat_hash_map<int, int> disc;
		VI delta;
		for (auto &pd_pii : pq_all) {
			if (pd_pii.first < t) break;
			PII uv = pd_pii.second;
			int num = dsu_merge.find(uv.first);
			if (!disc.contains(num)) {
				int _ = disc.size();
				disc.insert({ num, _ });
				delta.push_back(1);
			}
			else {
				delta[disc[num]]++;
			}
		}
		delta.push_back(0);
		for (int i = 1; i < delta.size(); i++) {
			delta[i] += delta[i - 1];
		}
		
		std::vector<PII> to_merge;
		to_merge.resize(cc);
		for (auto &pd_pii : pq_all) {
			if (pd_pii.first < t) break;
			PII uv = pd_pii.second;
			int num = disc[dsu_merge.find(uv.first)];
			to_merge[--delta[num]] = uv;
		}
		dsu_merge.init(n);
		
		std::set<int> to_update;
		VI to_remove;
#pragma omp parallel for reduction(merge_set_i: to_update) reduction(merge_vi: to_remove)
		for (int st = 1; st < delta.size(); st++) {
			for (int ss = delta[st] - 1; ss >= delta[st - 1]; ss--) {
				int u, v;
				std::tie(u, v) = to_merge[ss];
#pragma omp critical(dsu)
				u = dsu.find(u), v = dsu.find(v);
				if (u == v) continue;

				locks[u]->lock();
				phmap::flat_hash_map<int, int> u_A(cnt[u]);
				locks[u]->unlock();
				int cu = __get_cost(u_A, u, sz);

				locks[v]->lock();
				int cv = __get_cost(cnt[v], v, sz);
				for (auto pr : cnt[v]) {
					u_A[pr.first] += pr.second;
				}
				locks[v]->unlock();
				if (u_A.find(v) != u_A.end()) {
					u_A[u] += u_A[v];
					u_A.erase(v);
				}
				int cw = __get_merge_cost(u_A, u, v, sz);
				double score = double(cu + cv - cw) / (cu + cv);

				if (score < t || dsu.fa[u] != u || dsu.fa[v] != v) continue;
				// merge u v
#pragma omp critical(merge)
				{
					dsu.joint(u, v);
					sz[u] += sz[v];
					sz[v] = 0;
				}

				locks[v]->lock();
				auto tmpcntv = cnt[v];
				for (auto &pr : cnt[v]) {
					pr.second = 0;
				}
				locks[v]->unlock();

				for (auto &pr : tmpcntv) {
					if (pr.first == v || pr.first == u) continue;
					locks[pr.first]->lock();
					cnt[pr.first][v] = 0;
					cnt[pr.first][u] += pr.second;
					locks[pr.first]->unlock();
				}

				locks[u]->lock();
				for (auto pr : tmpcntv) {
					cnt[u][pr.first] += pr.second;
				}
				if (cnt[u].find(v) != cnt[u].end()) {
					cnt[u][u] += cnt[u][v];
					cnt[u].erase(v);
				}
				locks[u]->unlock();
				to_remove.push_back(v);

				to_update.insert(u);
			}
		}
		for (auto v : to_remove) {
			int u = dsu.find(v);
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

			for (auto pr : cnt[v]) {
				cnt[pr.first].erase(v);
			}
			cnt[v].clear();
		}

		phmap::flat_hash_set<int> to_update2(to_update.begin(), to_update.end());
		to_update.clear();
		for (int uuu : to_update2) {
			if (dsu.find(uuu) != uuu) continue;
			to_update.insert(uuu);
			for (auto pid : cnt[uuu]) {
				to_update.insert(pid.first);
			}
		}

		VI to_update_v(to_update.begin(), to_update.end());
		std::vector<PD_PII> to_update_suv;
#pragma omp parallel for reduction(merge_vpd_pii: to_update_suv)
		for (int uuu : to_update_v) {
			int cu = __get_cost(cnt[uuu], uuu, sz);

			VI to_remove_pri;
			for (auto & pid : W[uuu]) {
				double old_suv = pid.second;
				int vvv = pid.first;
				if (to_update.count(vvv) && vvv < uuu) continue;

				int cv = __get_cost(cnt[vvv], vvv, sz);

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

				to_update_suv.push_back({ old_suv, {uuu, vvv} });
				if (new_suv <= -0.03) {
					to_remove_pri.push_back(vvv);
				}
				W[uuu][vvv] = new_suv;
			}
			for (auto vvv : to_remove_pri) {
				W[uuu].erase(vvv);
			}
		}
		for (auto &pd_pii : to_update_suv) {
			int uuu, vvv;
			std::tie(uuu, vvv) = pd_pii.second;
			auto it = pq_all.find({ pd_pii.first, min_pair(uuu, vvv) });
			if (it == pq_all.end()) continue;
			pq_all.erase(it);
			W[vvv].erase(uuu);
			if (!W[uuu].contains(vvv)) continue;
			W[vvv][uuu] = W[uuu][vvv];
			pq_all.insert({ W[uuu][vvv], min_pair(uuu, vvv) });
		}
	}

	for (int i = 0; i < n; i++) {
		delete locks[i];
	}
	for (int i = 0; i < thread_num; i++) {
		delete rngs[i];
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

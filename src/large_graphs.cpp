#include "large_graphs.h"

#include <iostream>

#include "bab.h"
#include "upper_bound.h"
#include "lower_bound.h"

using namespace std;

/********** Close merge sparse ***************/
pair<int, int> random_neighbors(int x, const Graph &g, RNG &rng)
{
	auto nbs = g.neighbors(x) | g.red_neighbors(x);
	vector<int> res;
	sample(nbs.begin(), nbs.end(), 
			back_inserter(res), 2 /*two elements*/, rng);
	return make_pair(res[0], res[1]);
}

pair<int, contr_seq> close_merge_sparse(const Graph &g_init, RNG &rng, int outer_it, int inner_it)
{
	contr_seq best_sol;
	int best_cost = INFTY;

	for (int oit = 0; oit < outer_it; ++oit)
	{
		auto g = g_init;
		contr_seq cur_sol;

		vector<int> deg_gt_2;
		while (g.actual_n() > 1)
		{
			deg_gt_2.clear();
			for (int u: g.vertices())
				if (g.total_deg(u) > 1)
					deg_gt_2.push_back(u);

			if (deg_gt_2.empty())
			{
				// Max degree is at most 1
				contr_seq matching;
				for (int u: g.vertices())
				{
					for (int v: g.neighbors(u))
						if (u < v)
							matching.emplace_back(u, v);

					for (int v: g.red_neighbors(u))
						if (u < v)
							matching.emplace_back(u, v);
				}

				for (auto [u, v]: matching)
				{
					g.merge_nohint(u, v);
					cur_sol.emplace_back(u, v);
				}

				// Max degree is now 0:
				// just merge all vertices
				vector<int> vxs(g.vertices().begin(), g.vertices().end());
				for (size_t i = 1; i < vxs.size(); ++i)
				{
					g.merge_nohint(vxs[0], vxs[i]);
					cur_sol.emplace_back(vxs[0], vxs[i]);
				}

				continue;
			}

			pair<int,int> best_uv(-1, -1);
			Si best_hint;
			for (int i = 0; i < inner_it; ++i)
			{
				int x = random_from(deg_gt_2, rng);
				auto [u, v] = random_neighbors(x, g, rng);

				auto tmp_hint = g.merge_cost(u, v);
				if (best_uv.first == -1 || tmp_hint.size() < best_hint.size())
				{
					best_uv = make_pair(u, v);
					best_hint = move(tmp_hint);
				}
			}
			
			g.merge(best_uv.first, best_uv.second, move(best_hint));
			cur_sol.push_back(best_uv);

			// auto tmp = g.kernelize();
			// cur_sol.insert(cur_sol.end(), tmp.begin(), tmp.end());
		}
	
		if (g.full_width() < best_cost)
		{
			best_cost = g.full_width();
			best_sol = move(cur_sol);
		}
	}

	return make_pair(best_cost, move(best_sol));
}

// Large graph heuristics
pair<int, contr_seq> best_heur_sparse(const Graph &g)
{
	static random_device rd;
	static RNG rng(rd());

	auto &&res = close_merge_sparse(g, rng, DEFAULT_OUTER_SP, DEFAULT_INNER_SP);
	// auto &&res2 = greedy_mincost_local(g);
	return res;
	// return min(res, res2);
}


contr_seq solve_large(const Graph &g)
{
	pair<int, contr_seq> ub = best_heur_sparse(g);
	int lb_size = 25;
	int lb = subgraph_lb(g, lb_size);
	while (ub.first > lb)
	{
		cerr 
			<< "ub: " << ub.first
			<< ", lb: " << lb
			<< ", lb_size: " << lb_size
			<< endl;
		ub = min(ub, best_heur_sparse(g));
		lb = max(lb, subgraph_lb(g, lb_size, lb));
		lb_size = min(lb_size + 1, BitGraph::VxContainer::MAX_SIZE);
	}

	return ub.second;
}
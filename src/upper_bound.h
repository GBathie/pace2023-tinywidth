#pragma once

#include "common.h"

#include <queue>

#include "params.h"

using std::vector;
// Tree heuristic: process the graph as if it were a tree.
// Works well for sparse, tree-like graphs
template <class G>
int aux_merge(int u, int v, G &g, contr_seq &sol)
{
	if (u == -1)
		return v;
	else if (v == -1)
		return u;
	else
	{
		g.merge_nohint(u, v);
		sol.emplace_back(u, v);
		return u;
	}
}

template <class G>
int tree_merge_aux(int u, vector<bool> &seen, G &g, contr_seq &sol)
{
	if (seen[u])
		return -1;
	seen[u] = true;


	int r = -1;
	// Todo: add random order here.
	for (int v: g.neighbors(u))
	{
		v = tree_merge_aux(v, seen, g, sol);
		r = aux_merge(r, v, g, sol);
	}

	for (int v: g.red_neighbors(u))
	{
		v = tree_merge_aux(v, seen, g, sol);
		r = aux_merge(r, v, g, sol);
	}

	aux_merge(u, r, g, sol);

	return u;
}

template <class G>
std::pair<int, contr_seq> tree_merge_no_copy(G &g, RNG &rng)
{
	vector<bool> seen(g.n, false);
	contr_seq sol;

	vector<int> order;
	for (int u: g.vertices())
		order.push_back(u);
	std::shuffle(order.begin(), order.end(), rng);

	for (int u: order)
		tree_merge_aux(u, seen, g, sol);


	// Merge the remaining isolated vertices
	order.clear();
	for (int u: g.vertices())
		order.push_back(u);

	for (size_t i = 1; i < order.size(); ++i)
	{
		g.merge_nohint(order[0], order[i]);
		sol.emplace_back(order[0], order[i]);
	}

	return make_pair(g.full_width(), sol);
}

// TODO: add iterations
template <class G>
std::pair<int, contr_seq> tree_merge(const G &g_init, RNG &rng, int it)
{
	auto g = g_init;
	auto &&[score, sol] = tree_merge_no_copy(g, rng);

	for (int i = 1; i < it; ++i)
	{
		g = g_init;
		auto &&[score2, sol2] = tree_merge_no_copy(g, rng);
		if (score2 < score)
		{
			score = score2;
			sol = sol2;
		}
	}
	return make_pair(score, sol);
}


// Close-merge heuristic: pick a vertex of degree at least 2,
// and merge two of its neigbhors.
// Works well because two vertices that share many neighbors have a good chance of being selected that way.

inline int random_from(const vector<int> &v, RNG &rng)
{
	std::uniform_int_distribution<int> unif(0, v.size() - 1);
	return v[unif(rng)];
}

template <class G>
std::pair<int, int> random_neighbors(int x, const G &g, RNG &rng)
{
	std::uniform_real_distribution unif(0.0, 1.0);
	auto nbs = g.neighbors(x);
	if (nbs.size() < 2 || g.red_deg(x) * unif(rng) >= nbs.size())
		nbs = g.red_neighbors(x);

	// TODO: fix iterators on long_bitsets
	vector<int> tmp;
	for (int u: nbs)
		tmp.push_back(u);

	vector<int> res;
	std::sample(tmp.begin(), tmp.end(), 
			std::back_inserter(res), 2 /*two elements*/, rng);
	return std::make_pair(res[0], res[1]);
}

template <class G>
std::pair<int, contr_seq> close_merge(const G &g_init, RNG &rng, int outer_it, int inner_it)
{	
	contr_seq best_sol;
	int best_cost = INFTY;

	for (int oit = 0; oit < outer_it; ++oit)
	{
		auto g = g_init;
		contr_seq cur_sol = g.kernelize();

		vector<int> deg_gt_2;
		while (g.actual_n() > 1)
		{
			deg_gt_2.clear();
			for (int u: g.vertices())
				if (g.deg(u) > 1 || g.red_deg(u) > 1)
					deg_gt_2.push_back(u);

			if (deg_gt_2.empty())
			{
				// Degree max is 2, use tree method
				auto &&[score_tree, sol_tree] = tree_merge_no_copy(g, rng);
				cur_sol.insert(cur_sol.end(), sol_tree.begin(), sol_tree.end());
				break;
			}

			std::pair<int,int> best_uv(-1, -1);
			typename G::VxContainer best_hint;
			for (int i = 0; i < inner_it; ++i)
			{
				int x = random_from(deg_gt_2, rng);
				auto [u, v] = random_neighbors(x, g, rng);

				auto tmp_hint = g.merge_cost(u, v);
				if (best_uv.first == -1 || tmp_hint.size() < best_hint.size())
				{
					best_uv = std::make_pair(u, v);
					best_hint = std::move(tmp_hint);
				}
			}
			
			g.merge(best_uv.first, best_uv.second, std::move(best_hint));
			cur_sol.push_back(best_uv);

			auto tmp = g.kernelize();
			cur_sol.insert(cur_sol.end(), tmp.begin(), tmp.end());
		}
	
		if (g.full_width() < best_cost)
		{
			best_cost = g.full_width();
			best_sol = std::move(cur_sol);
		}
	}

	return make_pair(best_cost, std::move(best_sol));
}

/***** Greedy mincost : merge the pair of vertices with smallest approximate fusion cost ****/
template <class G>
std::pair<int, contr_seq> greedy_mincost(const G &g_init)
{	
	auto g = g_init;
	contr_seq cur_sol = g.kernelize();

	while (g.actual_n() > 1)
	{
		std::pair<int,int> best_uv(-1, -1);
		typename G::VxContainer best_hint;
		for (int u: g.vertices())
		{
			for (int v: g.vertices())
			{
				if (u >= v)
					continue;
				auto tmp_hint = g.merge_cost(u, v);
				if (best_uv.first == -1 || tmp_hint.size() < best_hint.size())
				{
					best_uv = std::make_pair(u, v);
					best_hint = std::move(tmp_hint);
				}
			}
		}
		
		g.merge(best_uv.first, best_uv.second, std::move(best_hint));
		cur_sol.push_back(best_uv);

		auto tmp = g.kernelize();
		cur_sol.insert(cur_sol.end(), tmp.begin(), tmp.end());
	}

	return make_pair(g.full_width(), std::move(cur_sol));
}


// Grid heuristic
template <class G>
vector<vector<int>> floyd_warshall(const G &g)
{
	vector<vector<int>> res(g.n, vector<int>(g.n, INFTY));
	for (int u: g.vertices())
	{
		res[u][u] = 0;
		for (int v: g.neighbors(u))
			res[u][v] = 1;
		for (int v: g.red_neighbors(u))
			res[u][v] = 1;
	}
	for (int k = 0; k < g.actual_n(); ++k)
		for (int u: g.vertices())
			for (int v: g.vertices())
				res[u][v] = std::min(res[u][v], res[u][k] + res[k][v]);

	return res;
}


template <class G>
void update_floyd_warshall(vector<vector<int>> &res, int u, int v, const G &g)
{
    for (int w: g.vertices())
        res[u][w] = res[w][u] = std::min(res[u][w], res[v][w]);
    for (int w1: g.vertices())
		for (int w2: g.vertices())
			res[w1][w2] = std::min(res[w1][w2], res[w1][u] + res[u][w2]);
}

template <class G>
std::pair<int, contr_seq> greedy_mincost_local(const G &g_init)
{
	auto g = g_init;
	contr_seq cur_sol;

	auto ap_dist = floyd_warshall(g);

	auto cost = [&](int u, int v) -> std::pair<int,int> {
        auto merge_hint = g.merge_cost(u, v);
		auto new_red = merge_hint - g.red_neighbors(u) - g.red_neighbors(v);
        int resuv = 0;
        int resvu = 0;
        for (int w: new_red) {
            resuv += ap_dist[u][w];
            resvu += ap_dist[v][w];
        }
        int res = std::min(resuv, resvu);
        return std::make_pair(res + merge_hint.size(), merge_hint.size());
	};

	int prev = -1;
	while (g.actual_n() > 1)
	{
		std::pair<int, int> best_uv(-1, -1);
		std::pair<int, int> best_cost;
		if (prev < 0)
		{
			for (int u: g.vertices())
			{
				for (int v: g.vertices())
				{
					if (u >= v)
						continue;
					auto tmp_cost = cost(u, v);
					if (best_uv.first == -1 || tmp_cost < best_cost)
					{
						best_uv = std::make_pair(u, v);
						best_cost = tmp_cost;
					}
				}
			}
		}
		else
		{
            for (int u : g.vertices()) {
                for (int v : g.vertices()) {
                    if (u >= v) continue;
                    if (ap_dist[prev][u] > 1 || ap_dist[prev][v] > 1) {
                        if (prev == u || prev == v) { // When we are at the limit of a line, we continue
                            if (ap_dist[prev][u] + ap_dist[prev][v] > 2) continue;
                        } else { // Otherwise we pick vertices at distance 2
                            if (ap_dist[u][v] != 2) continue;
                        }
                    }
                    auto tmp_cost = cost(u, v);
                    if (best_uv.first == -1 || tmp_cost < best_cost)
                    {
                        best_uv = std::make_pair(u, v);
                        best_cost = tmp_cost;
                    }
                }
            }
		}
		
		g.merge_nohint(best_uv.first, best_uv.second);
		cur_sol.push_back(best_uv);

        prev = best_uv.first;
        update_floyd_warshall(ap_dist, best_uv.first, best_uv.second, g);

	}

	return make_pair(g.full_width(), std::move(cur_sol));
}


// Aggregate heuristics
template <class G>
std::pair<int, contr_seq> apply_heur(const G &g,
									int tree_it = DEFAULT_TREE, 
									int outer_it = DEFAULT_OUTER,
									int inner_it = DEFAULT_INNER)
{
	static std::random_device rd;
	static RNG rng(rd());
	auto &&res = greedy_mincost(g);
	auto &&p1 = tree_merge(g, rng, tree_it);
	if (p1.first < res.first)
		res = p1;

	auto &&p2 = close_merge(g, rng, outer_it, inner_it);
	if (p2.first < res.first)
		res = p2;

	auto &&p3 = greedy_mincost_local(g);
	if (p3.first < res.first)
		res = p3;

	return std::move(res);
}


template <class G>
std::pair<int, contr_seq> best_heur(const G &g)
{
	auto x1 = apply_heur(g);

	auto g2 = g;
	auto kernel_op = g2.kernelize_heur();
	auto x2 = apply_heur(g2);
	kernel_op.insert(kernel_op.end(), x2.second.begin(), x2.second.end());
	std::swap(kernel_op, x2.second);
	
	return std::min(x1, x2);
}
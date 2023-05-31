#pragma once

#include <vector>
#include <queue>
#include <iostream>
#include <cassert>

#include "common.h"
#include "params.h"

using std::cerr;
using std::endl;

template<class T, class Mem>
RetValue mem_bab_aux_lb_init(T& g, int lb, int &min_score, Mem &mem)
{
	int full_width = g.full_width();
	if (full_width >= min_score)
		return make_pair(INFTY, contr_seq());

	int w = g.cur_width();
	auto kernel_moves = g.kernelize();
	kernel_moves = contr_seq(kernel_moves.rbegin(), kernel_moves.rend());
	if (g.actual_n() == 1)
	{
		if (full_width < min_score)
			min_score = full_width;
		return make_pair(w, kernel_moves);
	}

	auto [it, success] = mem.try_emplace(g.get_key(), INFTY, contr_seq());
	if (success) // key was not already present
	{
		auto moves = g.options();
		assert(!moves.empty());
		for (auto &[u, v] : moves)
		{
			T &gp = g.contract(u, v);
			auto &&[score, sol] = mem_bab_aux_lb_init(gp, lb, min_score, mem);
			score = std::max(w, score);
			if (score < it->second.first)
			{
				it->second.first = score;	
				it->second.second = std::move(sol);	
				it->second.second.emplace_back(u, v);
				it->second.second.insert(it->second.second.end(), kernel_moves.begin(), kernel_moves.end());
			}
			if (full_width >= min_score)
				break;

			if (min_score <= lb)
				break;
		}
	}

	return it->second;
}



template<class T>
RetValue mem_bab_heur_with_ub_lb(T &g, int ub, int lb)
{
	int min_score = ub;
	MemType mem;
	auto &&[score, sol] = mem_bab_aux_lb_init(g, lb, min_score, mem);

	contr_seq res = contr_seq(sol.rbegin(), sol.rend());

	return make_pair(min_score, res);
}

template<class G>
contr_seq cc_bab_with_lb(const G &g)
{
	contr_seq res = g.sol();

	std::vector<bool> seen(g.n, false);
	std::vector<int> q;

	int lb = g.full_width();
	// We need to keep one representative for each CC to merge them at the end
	std::vector<int> repr;
	std::unordered_set<int> cc;
	// m[i] is the index in g of the i-th vertex of the CC,
	// m_inv[m[i]] is i
	std::vector<int> m(g.n), m_inv(g.n);
	auto kernelize_if_lb_gt2 = [&](BitGraph &h) {
		// Add some kernelization if available
		if (lb >= 2)
			for (auto [u, v]: h.kernelize_tww_gt2())
				res.emplace_back(m[u], m[v]);
	};
	for (int u: g.vertices())
	{
		cc.clear();

		q.push_back(u);
		while (!q.empty())
		{
			int v = q.back(); q.pop_back();
			if (seen[v])
				continue;
			seen[v] = true;

			m[cc.size()] = v;
			m_inv[v] = cc.size();
			cc.insert(v);

			for (int w: g.neighbors(v))
				q.push_back(w);

			for (int w: g.red_neighbors(v))
				q.push_back(w);
		}

		if (cc.size() > 0)
		{
			auto h = g.dense_subgraph(cc, m_inv);
			kernelize_if_lb_gt2(h);

			// Compute upper bound
			auto&& [ub, sol] = best_heur(h);
			
			// Compute lb with time proportional to the size
			int cc_lb_time = (h.actual_n() * LB_TIME_S) / g.actual_n();
			int cc_lb = timed_iter_subgraph_lb_early_exit(h, ub, lb, LB_K, cc_lb_time);
			lb = std::max(lb, cc_lb);

			contr_seq sol2;
			if (lb >= 2) 
				sol2 = h.kernelize_tww_gt2();

			std::cerr << "n: " << h.actual_n()
					  << ", Ub: " << ub
					  << ", cc_lb: " << cc_lb
					  << ", lb: " << lb << std::endl;

			auto&& [h_score, h_res] = mem_bab_heur_with_ub_lb(h, ub, lb);
			sol2.insert(sol2.end(), h_res.begin(), h_res.end());

			std::cerr << "Ub: " << ub
					  << ", cc_lb: " << cc_lb
					  << ", lb: " << lb
					  << ", bab: " << h_score << std::endl;
			if (h_score < ub)
				sol = sol2;
			
			// Propagate lb to other ccs:
			// if some CC has an optimal value of W,
			// other CCs do not need to do better.
			lb = std::max(lb, std::min(h_score, ub));

			for (auto [u, v]: sol)
				res.emplace_back(m[u], m[v]);

			if (cc.size() == 1)
				repr.push_back(*cc.begin());
			else
				repr.push_back(res.back().first);
		}
	}

	// Merge all ccs
	for (size_t i = 1; i < repr.size(); ++i)
		res.emplace_back(repr[i], repr[i - 1]);

	return res;
}

#pragma once

#include <random>
#include <queue>
#include <vector>
#include <algorithm>
#include <cassert>
#include <chrono>

#include "common.h"

template <class ItType, class RNG>
int reservoir_sampling(ItType begin, ItType end, RNG &rng)
{
	std::uniform_real_distribution unif(0.0, 1.0);
	int c = -1;
	double w = 0;
	for (ItType it = begin; it != end; ++it)
	{
        w += 1;
        if (w * unif(rng) < 1)
            c = *it;

	}
	assert(c != -1);
	return c;
}

/*
 * Property: for any subgraph H of G, tww(H) <= tww(G).
 * We use this as a lower bound: take a random subgraph of size k
 * and compute its twin-width.
 */
template <class T>
int subgraph_lb(const T &g, int k, int prev_lb = 0)
{
	static std::random_device rd;
	static RNG rng(rd());
	std::uniform_real_distribution unif(0.0, 1.0);

	std::vector<bool> seen(g.n, false);
	typename T::VxContainer h;
	std::priority_queue<std::pair<float,int>> q;

	// start with a random vertex
	const auto &tmp = g.vertices();
	int s = reservoir_sampling(tmp.begin(), tmp.end(), rng);
	q.emplace(unif(rng), s);

	while ((int)h.size() < k && !q.empty())
	{
		auto [_, v] = q.top(); q.pop();
		if (seen[v])
			continue;
		h.insert(v);
		seen[v] = true;

		for (int u: g.neighbors(v))
			q.emplace(unif(rng), u); // add u with random prioriy
		
		for (int u: g.red_neighbors(v))
			q.emplace(unif(rng), u); // add u with random prioriy
	}

	// Solve on subgraph containing vertices of h
	auto H = g.subgraph(h);
	int min_score = INFTY;
	MemType mem;
	// Here we use the following trick: prev_lb can be seen as an lower bound for
	// this graph: if we find a better smaller solution for H
	// its optimal solution will not be a better lower bound.
	auto &&[res, _] = mem_bab_aux_lb_init(H, prev_lb, min_score, mem);

	return (res == INFTY) ? prev_lb : res;
}

/*
 * Return the best of `it` runs of subgraph_lb of size k.
 */
template <class T>
int iter_subgraph_lb(const T &g, int k, int it)
{
	int best = 0;
	for (int i = 0; i < it; ++i)
		best = std::max(best, subgraph_lb(g, k));
	return best;
}


/*
 * Return the best subgraph_lb of size k found within `sec_max` seconds.
 */

using namespace std::chrono;

template <class T>
int timed_iter_subgraph_lb(const T &g, int k, int sec_max)
{
	int best = 0;
	auto start = high_resolution_clock::now();
	while (duration_cast<seconds>(high_resolution_clock::now() - start).count() < sec_max)
	{
		best = std::max(best, subgraph_lb(g, k, best));
	}
	return best;
}

template <class T>
int timed_iter_subgraph_lb_early_exit(const T &g, int best_score, int lb0, int k, int sec_max)
{
	int best_lb = lb0;
	auto start = high_resolution_clock::now();
	while (duration_cast<seconds>(high_resolution_clock::now() - start).count() < sec_max)
	{
		best_lb = std::max(best_lb, subgraph_lb(g, k));
		if (best_lb >= best_score)
			break;
	}
	return best_lb;
}

/*
 * Apply subgraph_lb starting from a small size, for at most `sec_max` seconds.
 * When the lb is equal to the best so far, we increase the size.
 * Otherwise, we just do another run.
 */
template <class T>
int timed_growing_subgraph_lb(const T &g, int s0, int sec_max)
{
	int best = 0;
	auto start = high_resolution_clock::now();
	while (duration_cast<seconds>(high_resolution_clock::now() - start).count() < sec_max)
	{
		int lb = subgraph_lb(g, s0, best);
		if (lb <= best)
		{
			s0 += 1;
			std::cerr << "New lb size: " << s0 << std::endl;
		}
		else
		{
			best = lb;
			std::cerr << "New lb value: " << lb << std::endl;
		}

		if (s0 >= g.actual_n())
			break;
	}
	return best;
}


template<class T>
int greedy_lb(const T &g)
{
	int lb = INFTY;
	for (int u: g.vertices())
		for (int v: g.vertices())
			if (u != v)
			{
				auto tmp = (g.neighbors(u) ^ g.neighbors(v));
				tmp.erase(u);
				tmp.erase(v);
				if (tmp.size() < lb)
					lb = tmp.size();
			}

	return lb;
}
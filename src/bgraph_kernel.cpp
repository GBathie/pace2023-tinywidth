#include "bgraph.h"

#include "union_find.hpp"

using namespace std;

contr_seq BitGraph::kernelize()
{
	contr_seq res;
	bool cont = true;
	while (cont)
	{
		cont = find_dominating(res);
	}
	return res;
}

// Tests whether u dominates v,
// i.e. whether merging v into u creates new red edges or not.
// This happens when redN(v) \cup (N(u) \Delta N(v)) is a subset of redN(u)
// and N(u) is a subset of N(v)
// This generalizes twinness to vertices with red edges
bool BitGraph::dominates(int u, int v) const
{
	VxContainer nu = adj[u];
	nu.erase(v);
	VxContainer rnu = red_adj[u];
	rnu.erase(v);

	VxContainer nv = adj[v];
	nv.erase(u);
	VxContainer rnv = red_adj[v];
	rnv.erase(u);
	
	return (nu <= nv) && (((nv - nu) | rnv) <= rnu);
}

/*
 * Faster version of find_one_dominating:
 * we take a vertex u, and take the intersection of the neighborhoods
 * of its neighbors. Any vertex v != u in the resulting set contains N[u]:
 * it is a candidate for domination.
 */
bool BitGraph::find_one_dominating(contr_seq &seq)
{
	for (int u: vertices())
	{
		VxContainer candidates = vertices() - VxContainer::singleton(u);
		for (int w: neighbors(u))
			candidates &= neighbors(w);

		bool erased = false;
		for (int v: candidates)
			if (dominates(u, v))
			{
				merge_dominating(u, v, seq);
				erased = true;
			}
		if (erased)
			return true;
	}
	return false;
}

bool BitGraph::find_dominating(contr_seq &seq)
{

	bool cont = true, res = false;
	while (cont)
	{
		cont = find_one_dominating(seq);
		res = res || cont;
	}

	return res;
}


void BitGraph::merge_dominating(int u, int v, contr_seq &seq)
{
	erase(v);
	update_key(u, v);
	// Update sol
	seq.emplace_back(u, v);
}



/********** Approximate kernelization ***********/
contr_seq BitGraph::kernelize_tww_gt2()
{
	contr_seq res;
	
	reduce_trees(res);

	return res;
}


void BitGraph::reduce_trees(contr_seq &seq)
{
	vector<int> parent(n, -1);
	vector<vector<int>> children(n);

	// We first build the tree:
	// loop over candidate leaves
	// add ther only non-leaf neighbor as their parent,
	// and add it to the stack.
	vector<int> q1, q2;
	for (int u: vertices())
		if (total_deg(u) == 1)
			q1.push_back(u);

	while (!q1.empty())
	{
		for (int u: q1)
		{
			if (parent[u] != -1)
				continue;
			if (total_deg(u) != ((int)children[u].size() + 1))
				continue;
			for (int v: neighbors(u))
				if (parent[v] == -1)
				{
					children[v].push_back(u);
					parent[u] = v;
					q2.push_back(v);
					break;
				}
			for (int v: red_neighbors(u))
				if (parent[v] == -1)
				{
					children[v].push_back(u);
					parent[u] = v;
					q2.push_back(v);
					break;
				}
		}

		swap(q1, q2);
		q2.clear();
	}

	auto aux_merge = [&](int u, int v) {
		if (u == -1)
			return v;
		else if (v == -1)
			return u;
		merge_nohint(u, v);
		seq.emplace_back(u, v);
		return u;
	};

	// Reduces a node to have a single child,
	// and return the id of that child.
	function<int(int)> merge_tree = [&](int u) -> int {
		int c = -1, gc = -1;
		// for each children v, reduce it to depth 1,
		// merge its grandchild w with gc,
		// and merge v with c
		for (int v: children[u])
		{
			int w = merge_tree(v);
			gc = aux_merge(gc, w);
			c = aux_merge(c, v);
		}

		// if reducing to depth 1,
		// merge grandchildren with children
		if (parent[u] != -1)
		{
			c = aux_merge(c, gc);
			return c;
		}
		else
			return -1;
	};

	for (int u = 0; u < n; ++u)
		if (!is_deleted(u) && parent[u] == -1)
			merge_tree(u);
}


/******* Path reduction *******/
contr_seq BitGraph::kernelize_heur()
{
	contr_seq res;
	reduce_paths(res);
	reduce_trees(res);
	reduce_paths(res);
	return res;
}

void BitGraph::reduce_paths(contr_seq &seq)
{
	vector<bool> seen(n, false);
	vector<vector<int>> to_merge(n);
	for (int u: vertices())
		reduce_paths_aux(u, seen, to_merge);

	OrderedUnionFind uf(n);
	for (int u = 0; u < n; ++u)
	{
		for (int v: to_merge[u])
		{
			merge_nohint(uf.find(u), uf.find(v));
			seq.emplace_back(uf.find(u), uf.find(v));
			uf.merge(u, v);
		}
	}
}

// Returns true iff u is a mid-path vertex.
// A vertex is *mid-path* if and only if
// it has degree at most 2 and so do 
// all of its (at most 2) neighbors.
bool BitGraph::reduce_paths_aux(int u, vector<bool> &seen,
		vector<vector<int>> &to_merge)
{
	if (total_deg(u) > 2)
		return false;
	if (seen[u])
		return false;
	seen[u] = true;

	bool midpoint = true;
	for (int v: neighbors(u))
	{
		if (reduce_paths_aux(v, seen, to_merge))
			to_merge[u].push_back(v);
		if (total_deg(v) > 2)
			midpoint = false;
	}

	for (int v: red_neighbors(u))
	{
		if (reduce_paths_aux(v, seen, to_merge))
			to_merge[u].push_back(v);
		if (total_deg(v) > 2)
			midpoint = false;
	}
	if (!midpoint)
		to_merge[u].clear();

	return midpoint;
}
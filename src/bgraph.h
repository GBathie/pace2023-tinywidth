#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include <set>
#include <random>

#include "common.h"

#include "long_bitset.hpp"

class BitGraph
{
public:
	int n;
    static std::vector<BitGraph> instances;
	explicit BitGraph(int n, int past_tww = 0);

	using VxContainer = LongBitset<1>;

	inline const std::string &get_key() const { return key; }
	// BaB functions
	inline int full_width() const { return full_tww; }
	inline int cur_width() const { return cur_tww; }
	BitGraph &contract(int u, int v) const;

	contr_seq options() const;

	VxContainer merge_cost(int u, int v) const;
	void merge(int u, int v, VxContainer &&hint);
	inline void merge_nohint(int u, int v) { merge(u, v, merge_cost(u, v)); };

	BitGraph subgraph(VxContainer h) const;

	contr_seq kernelize();
	contr_seq kernelize_tww_gt2();
	contr_seq kernelize_heur();

	inline VxContainer vertices() const { return vertex_mask; }
	inline int actual_n() const { return vertex_mask.size(); }
	inline int deg(int u) const { return adj[u].size(); }
	inline int red_deg(int u) const { return red_adj[u].size(); }
	inline int total_deg(int u) const { return deg(u) + red_deg(u); }
	inline bool adjacent(int u, int v) const { return adj[u].contains(v) || red_adj[u].contains(v); }
	inline VxContainer neighbors(int u)     const { return adj[u]; }
	inline VxContainer red_neighbors(int u) const { return red_adj[u]; }
	inline bool is_deleted(int u) const { return !vertex_mask.contains(u); }
	inline VxContainer non_neighbors(int u) const { return ((~adj[u]) & vertex_mask) - VxContainer::singleton(u); }

    inline static void init(int n) { instances = std::vector<BitGraph>(n, BitGraph(n)); }

	friend std::ostream &operator<<(std::ostream &os, const BitGraph &g);
	friend class Graph;

	template<class Func>
	inline void iter_nodes(Func &&f) const
	{
		for (int u: vertex_mask)
			f(u);
	}

	template<class Func>
	inline void iter_pair_diff_lt_nodes(Func &&f) const
	{
		for (int u: vertex_mask)
			for (auto it = vertex_mask.begin(u+1); it != vertex_mask.end(); ++it)
				f(u, *it);
	}
private:
	int full_tww;
    int cur_tww;
	std::vector<VxContainer> adj;
	std::vector<VxContainer> red_adj;
	VxContainer vertex_mask;
	std::string key;

	void add_edge(int u, int v, bool red = false);
	void erase_edge(int u, int v, bool red = false);
	void erase(int u);
	void update_key(int u, int v);
	void compute_width();
    void copy(const BitGraph &from);
	
	bool dominates(int u, int v) const;
	bool find_one_dominating(contr_seq &seq);
	bool find_dominating(contr_seq &seq);
	void merge_dominating(int u, int v, contr_seq &seq);
	bool kernelize_stars();

	void reduce_trees(contr_seq &seq);
	void reduce_paths(contr_seq &seq);
	bool reduce_paths_aux(int u, std::vector<bool> &seen, std::vector<std::vector<int>> &to_merge);
};

inline std::ostream &operator<<(std::ostream &os, const BitGraph &g)
{
	os << "----" << std::endl;
	for (int i = 0; i < g.n; ++i)
	{
		if (g.is_deleted(i))
			continue;
		os << i << ", " << g.adj[i] << std::endl; 
		os << i << ", " << g.red_adj[i] << std::endl; 
		os << "----" << std::endl;
	}
	return os;
}

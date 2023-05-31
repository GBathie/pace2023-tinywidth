#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <unordered_set>

#include "common.h"
#include "bgraph.h"

using Vb = std::vector<bool>;
using Si = std::unordered_set<int>;

class Graph
{
public:
	int n;
	Graph(int n);

	using VxContainer = Si;
	int full_width() const { return full_tww; }
	const contr_seq &sol() const { return seq; }

	inline const Si &vertices() const { return vertex_mask; }
	inline int actual_n() const { return vertex_mask.size(); }
	inline int deg(int u) const { return adj[u].size(); }
	inline int red_deg(int u) const { return red_adj[u].size(); }
	inline int total_deg(int u) const { return adj[u].size() + red_adj[u].size(); }
	inline bool adjacent(int u, int v) const { return adj[u].contains(v) || red_adj[u].contains(v); }
	inline const Si &neighbors(int u)     const { return adj[u]; }
	inline const Si &red_neighbors(int u) const { return red_adj[u]; }
	inline bool is_deleted(int u) const { return !vertex_mask.contains(u); }

	Si merge_cost(int u, int v) const;
	void merge(int u, int v, Si &&hint);
	void merge_nohint(int u, int v);

	// Kernelization
	contr_seq kernelize();
	void kernelize_safe();
	contr_seq kernelize_heur();

	static Graph from_istream(std::istream &is);
	inline static Graph from_file(const std::string &fname) { std::ifstream ifs(fname); return from_istream(ifs); }
	inline static Graph from_cin() { return from_istream(std::cin); }


	std::vector<std::pair<std::vector<int>, BitGraph>> connected_components() const;
	BitGraph subgraph(const Si& vx) const;
	BitGraph dense_subgraph(const Si& vx, const std::vector<int> &m_inv) const;
    int largest_cc_size() const;

	friend std::ostream &operator<<(std::ostream &os, const Graph &g);
private:
	int full_tww;
	// int cur_tww;
	std::vector<Si> adj;
	std::vector<Si> red_adj;
	Si vertex_mask;
	contr_seq seq;

	void add_edge(int u, int v);
	void erase_edge(int u, int v, bool red = false);
	void erase(int u);

	// Kernelization aux
	void merge_twins();
	bool reduce_paths_aux(int u, std::vector<bool> &seen,
		std::vector<std::vector<int>> &to_merge);
	void reduce_paths();
	void reduce_trees();
	void reduce_bulls();
};


inline std::ostream &operator<<(std::ostream &os, const Graph &g)
{
	os << "----" << std::endl;
	for (int i : g.vertices())
	{
		os << i << ", " << g.adj[i] << std::endl; 
		os << i << ", " << g.red_adj[i] << std::endl; 
		os << "----" << std::endl;
	}
	return os;
}
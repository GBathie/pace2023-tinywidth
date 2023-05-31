#pragma once

#include <vector>

class OrderedUnionFind
{
public:
	OrderedUnionFind(int n):
		parent(n),
		size(n, 1) {
		for (int i = 0; i < n; ++i)
			parent[i] = i;
	}

	int find(int u)	{
		if (parent[u] != u)
			parent[u] = find(parent[u]);
		return parent[u];
	}

	bool merge(int u, int v) {
		u = find(u);
		v = find(v);
		if (u == v)
			return false;

		parent[v] = u;
		size[u] += size[v];
	
		return true;
	}

	std::vector<int> parent;	
	std::vector<int> size;
};
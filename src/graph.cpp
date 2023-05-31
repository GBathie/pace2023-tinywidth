#include "graph.h"

#include <cassert>
#include <unordered_map>
#include <stack>

using namespace std;

Graph::Graph(int n):
	n(n),
	full_tww(0), 
	adj(n), red_adj(n),
	vertex_mask(), seq()
{
	for (int i = 0; i < n; ++i)
		vertex_mask.insert(i);
}

void Graph::merge_nohint(int u, int v)
{
	assert(u != v);
	// Remove edge u-v if it exists, to avoid self-loops
	erase_edge(u, v);
	erase_edge(u, v, true);

	red_adj[u] = red_adj[u] | red_adj[v] | (adj[u] ^ adj[v]);
	adj[u] = adj[u] & adj[v];

	// Update adjacencies on the other end of edges
	for (int w: red_adj[u])
	{
		red_adj[w].insert(u);
		adj[w].erase(u);
	}

	erase(v);

	// Update width
	full_tww = max(full_tww, red_deg(u));
	for (int w: red_adj[u])
		full_tww = max(full_tww, red_deg(w));
	// compute_width();
	// for (int w: vertex_mask)
		// assert(full_tww >= red_deg(w));

	// Update sol
	seq.emplace_back(u, v);
}


Si Graph::merge_cost(int u, int v) const
{
	auto tmp = red_adj[u] | red_adj[v] | (adj[u] ^ adj[v]);
	tmp.erase(u);
	tmp.erase(v);
	return tmp;
}

void Graph::merge(int u, int v, Si &&hint)
{
	// Remove edge u-v if it exists, to avoid self-loops
	erase_edge(u, v);
	erase_edge(u, v, true);

	red_adj[u] = std::move(hint);
	adj[u] = adj[u] & adj[v];

	// Update adjacencies on the other end of edges
	for (int w: red_adj[u])
	{
		red_adj[w].insert(u);
		adj[w].erase(u);
	}

	erase(v);

	// Update width
	full_tww = max(full_tww, red_deg(u));
	for (int w: red_adj[u])
		full_tww = max(full_tww, red_deg(w));
}


/*************** Private Methods ***************/
void Graph::add_edge(int u, int v)
{
	adj[u].insert(v);
	adj[v].insert(u);
}

void Graph::erase_edge(int u, int v, bool red)
{
	(red ? red_adj[u] : adj[u]).erase(v);
	(red ? red_adj[v] : adj[v]).erase(u);
}

void Graph::erase(int u)
{
	for (int v: adj[u])
		adj[v].erase(u);

	for (int v: red_adj[u])
		red_adj[v].erase(u);
	
	vertex_mask.erase(u);
}

/*************** Parsing Graphs ***************/
Graph Graph::from_istream(istream &is)
{
	string a;
	int n, m;
	is >> a >> a >> n >> m;
	Graph res(n);
	
	int u,v;
	for (int i = 0; i < m; ++i)
	{
		is >> u >> v;
		res.add_edge(u-1, v-1);
	}

	return res;
}

BitGraph Graph::subgraph(const Si& vx) const
{
	vector<int> m_inv(n);
	int i = 0;
	for (auto u: vx)
		m_inv[u] = i++;

	return dense_subgraph(vx, m_inv);
}

int Graph::largest_cc_size() const
{
    vector<bool> visited(n, false);
    int res = 0;
    int current_res = 0;
    stack<int> to_visit;
    for (int u : vertices())
    {
        if (visited[u]) continue;
        current_res = 0;
        to_visit.push(u);
        visited[u] = true;
        while (!to_visit.empty())
        {
            current_res++;
            auto v = to_visit.top();
            to_visit.pop();
            for (auto w : neighbors(v) | red_neighbors(v))
            {
                if (visited[w]) continue;
                visited[w] = true;
                to_visit.push(w);
            }
        }
        res = max(res, current_res);
    }
    return res;
}


BitGraph Graph::dense_subgraph(const Si& vx, const vector<int> &m_inv) const
{
	BitGraph res(vx.size(), full_tww);
	for (int u: vx)
	{
		for (int v: neighbors(u))
			if (vx.count(v) > 0)
				res.add_edge(m_inv[u], m_inv[v]);

		for (int v: red_neighbors(u))
			if (vx.count(v) > 0)
				res.add_edge(m_inv[u], m_inv[v], true);
	}
	res.compute_width();

	return res;
}

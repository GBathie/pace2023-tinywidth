#include "bgraph.h"

#include <cassert>

using namespace std;

static std::random_device rd;
static RNG rng(rd());

vector<BitGraph> BitGraph::instances = vector<BitGraph>();

BitGraph::BitGraph(int n, int past_tww):
	n(n),
	full_tww(past_tww), cur_tww(0),
	adj(n), red_adj(n),
	vertex_mask(VxContainer::full(n)),
	key(n, 0)
{
	for (int i = 0; i < n; ++i)
		key[i] = i;
	assert(n <= VxContainer::MAX_SIZE);
}

void BitGraph::copy(const BitGraph &from) {
    full_tww = from.full_tww;
    cur_tww = from.cur_tww;
    std::copy(from.adj.begin(), from.adj.end(), adj.begin());
    std::copy(from.red_adj.begin(), from.red_adj.end(), red_adj.begin());
    vertex_mask = from.vertex_mask;
    std::copy(from.key.begin(), from.key.end(), key.begin());
}

BitGraph &BitGraph::contract(int u, int v) const
{
	BitGraph &res = BitGraph::instances[actual_n()-1];
    res.copy(*this);
	res.cur_tww = 0;
	res.merge_nohint(u, v);
	return res;
}

contr_seq BitGraph::options() const
{
	vector<tuple<int,int,int>> tmp;

	iter_pair_diff_lt_nodes([&](int u, int v) { tmp.emplace_back(-red_deg(u) - red_deg(v), u, v); });
	std::sort(tmp.begin(), tmp.end());

	contr_seq res;
	for (auto &&[_, u, v]: tmp)
		res.emplace_back(u, v);

	return res;
}


BitGraph::VxContainer BitGraph::merge_cost(int u, int v) const
{
	auto tmp = red_adj[u] | red_adj[v] | (adj[u] ^ adj[v]);
	tmp.erase(u);
	tmp.erase(v);
	return tmp;
}



BitGraph BitGraph::subgraph(BitGraph::VxContainer h) const
{
	BitGraph res = *this;
	res.vertex_mask = h;
	for (auto &a: res.adj)
		a &= h;
	for (auto &a: res.red_adj)
		a &= h;

	res.cur_tww = 0;
	res.compute_width();
	
	return res;
}

void BitGraph::compute_width()
{
	for (int u: vertex_mask)
		cur_tww = max(cur_tww, red_deg(u));
	full_tww = max(full_tww, cur_tww);
}
/*************** Private Methods ***************/
void BitGraph::add_edge(int u, int v, bool red)
{
	assert(u != v);
	((red) ? red_adj[u] : adj[u]).insert(v);
	((red) ? red_adj[v] : adj[v]).insert(u);
}


void BitGraph::erase_edge(int u, int v, bool red)
{
	assert(u != v);
	((red) ? red_adj[u] : adj[u]).erase(v);
	((red) ? red_adj[v] : adj[v]).erase(u);
}

void BitGraph::erase(int u)
{
	for (int v: adj[u])
		adj[v].erase(u);

	for (int v: red_adj[u])
		red_adj[v].erase(u);
	
	vertex_mask.erase(u);
}

void BitGraph::merge(int u, int v, VxContainer &&hint)
{
	// Remove edge u-v if it exists, to avoid self-loops
	erase_edge(u, v);
	erase_edge(u, v, true);

	red_adj[u] = std::move(hint);
	adj[u] &= adj[v];

	// Update adjacencies on the other end of edges
	for (int w: red_adj[u])
	{
		red_adj[w].insert(u);
		adj[w].erase(u);
	}

	erase(v);

	// Update width
	cur_tww = max(cur_tww, red_deg(u));
	for (int w: red_adj[u])
		cur_tww = max(cur_tww, red_deg(w));
	full_tww = max(full_tww, cur_tww);

	update_key(u, v);
}

void BitGraph::update_key(int u, int v)
{
	char cu = key[u];
	char cv = key[v];
	if (cu > cv)
		swap(cu, cv);
	assert(cu != cv);
	// cu is now smaller
	for (char &c: key)
		if (c == cv)
			c = cu;
}
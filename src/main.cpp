#include "common.h"
#include "long_bitset.hpp"
#include "bgraph.h"
#include "graph.h"
#include "bab.h"
#include "upper_bound.h"
#include "lower_bound.h"
#include "large_graphs.h"

using namespace std;

void print_sol(const contr_seq &sol)
{
	for (auto &[u, v]: sol)
		cout << u+1 << " " << v+1 << "\n";

	cout << flush;
}

void solve_cin()
{
	auto g = Graph::from_cin();
	g.kernelize_safe();
    BitGraph::init(min(g.n, BitGraph::VxContainer::MAX_SIZE));

	int max_cc_size = g.largest_cc_size();

	contr_seq sol;
	if (max_cc_size > BitGraph::VxContainer::MAX_SIZE)
	{
		cerr << "Starting large graphs branch" << endl;
		sol = solve_large(g);
	}
	else
	{
		cerr << "Starting dense graphs branch" << endl;
		sol = cc_bab_with_lb(g);
	}
	print_sol(sol);
}

int main()
{
	solve_cin();

	return 0;
}
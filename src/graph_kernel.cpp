#include "graph.h"

#include <cassert>
#include <functional>
#include <unordered_map>

#include "union_find.hpp"

using namespace std;

/******* Kernelization ********/
void Graph::kernelize_safe()
{
	merge_twins();
}

/******* Fast twins finding *******/
struct TrieNode
{
	unordered_map<int, TrieNode*> children;
	vector<int> twins;

	TrieNode(): children(), twins() { }
	
	~TrieNode()
	{
		for (auto [k, v]: children)
			delete v;
	}

	TrieNode *child(int u)
	{
		auto &tmp = children[u];
		if (tmp == nullptr)
			tmp = new TrieNode();
		return tmp;
	}
	
	// Iterative version of insertion
	static void insert(TrieNode* root, int u, const vector<int> &adj)
	{
		for (int x: adj)
			root = root->child(x);

		root->twins.push_back(u);
	}

	void merge_twins(Graph &g)
	{
		for (auto [k, v]: children)
			v->merge_twins(g);
		for (size_t	i = 1; i < twins.size(); ++i)
			g.merge_nohint(twins[0], twins[i]);
	}
};

void Graph::merge_twins()
{
	int starting_n;
	do 
	{
		starting_n = actual_n();
		// True twins
		TrieNode root;
		for (int u: vertices())
		{
			vector<int> nb(adj[u].begin(), adj[u].end());
			nb.push_back(u); // To handle true twins
			sort(nb.begin(), nb.end());
			TrieNode::insert(&root, u, nb);
		}
		root.merge_twins(*this);

		// False twins
		TrieNode root2;
		for (int u: vertices())
		{
			vector<int> nb(adj[u].begin(), adj[u].end());
			sort(nb.begin(), nb.end());
			TrieNode::insert(&root2, u, nb);
		}
		root2.merge_twins(*this);

	} while (actual_n() < starting_n);
}
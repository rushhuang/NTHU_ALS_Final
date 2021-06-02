#include <iostream>
#include <climits> // INT_MAX
using namespace std;

#include <LEDA/graph/graph.h>
#include <LEDA/graph/min_cut.h>
#include <LEDA/graph/min_cost_flow.h>
using namespace leda;

int main()
{
	// Graph construction
	graph G;
	node n0 = G.new_node();
	node n1 = G.new_node();
	node n2 = G.new_node();
	node n3 = G.new_node();
	edge e0 = G.new_edge(n0, n1);
	edge e1 = G.new_edge(n1, n3);
	edge e2 = G.new_edge(n0, n2);
	edge e3 = G.new_edge(n2, n3);

	// Assign edge weights
	edge_array<int> weight(G);
	weight[e0] = 1;
	weight[e1] = 3;
	weight[e2] = 2;
	weight[e3] = 2;

	// Min cut algorithm
	list<node> cut;
	int cut_value = MIN_CUT(G, weight, cut);

	cout << "The minimum cut has value: " << cut_value << endl;
	cout << "cut:";
	node v;
	forall(v, cut)
		G.print_node(v);
	cout << endl;

	return 0;
}

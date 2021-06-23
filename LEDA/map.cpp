// LEDA Manual
// http://www.algorithmic-solutions.info/leda_manual/MANUAL.html

#include <iostream>
#include <string>
#include <sstream> // convert sting to int and vice versa
#include <fstream>
#include <vector>
#include <map>
#include <queue>
#include <stack>
#include <climits> // INT_MAX
#include <algorithm> // std::remove
using namespace std;

#include <LEDA/core/list.h>
#include <LEDA/graph/graph.h>
#include <LEDA/graph/min_cut.h>
#include <LEDA/graph/min_cost_flow.h>
using namespace leda;

// ASSERT macro for assert with error message
// From: https://stackoverflow.com/questions/3767869/adding-message-to-assert
#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

struct NODE{
	int NODE_level;

	std::string NODE_name;

	// Define NODE type
	// 0: PI, 1: AND, 2: OR, 3: INVERTER, 4: CONSTANT_0, 5: CONSTANT_1, 6: BUFFER, 7: UNKNOWN
	std::string NODE_type;

	// Input nodes
	std::vector<std::string> InNODEs;

	// Output nodes
	std::vector<std::string> OutNODEs;
};

// Reference:
// https://www.geeksforgeeks.org/cpp-program-for-topological-sorting/
// A recursive function used by topologicalSort
void topologicalSortUtil(std::map<std::string, NODE>& NODES_map, std::vector<NODE>& Total_NODES,
						 std::string v, std::map<std::string, bool>& visited, std::stack<std::string>& Stack){
    // Mark the current node as visited.
    visited[v] = true;
  
    // Recur for all the vertices adjacent to this vertex
    for (std::vector<std::string>::iterator i = NODES_map[v].OutNODEs.begin(); i != NODES_map[v].OutNODEs.end(); ++i)
        if (!visited[*i])
            topologicalSortUtil(NODES_map, Total_NODES, *i, visited, Stack);
  
    // Push current vertex to stack which stores result
    Stack.push(v);
}

// The function to do Topological Sort. It uses recursive topologicalSortUtil()
void topologicalSort(std::stack<std::string>& Stack, std::map<std::string, NODE>& NODES_map, std::vector<NODE>& Total_NODES){  
    // Mark all the vertices as not visited
    // std::map<std::string, bool>* visited = new std::map<std::string, bool>;
    std::map<std::string, bool> visited;
    for (int i = 0; i < Total_NODES.size(); i++)
        visited[Total_NODES[i].NODE_name] = false;
  
    // Call the recursive helper function to store Topological
    // Sort starting from all vertices one by one
    for (int i = 0; i < Total_NODES.size(); i++)
        if (visited[Total_NODES[i].NODE_name] == false)
            topologicalSortUtil(NODES_map, Total_NODES, Total_NODES[i].NODE_name, visited, Stack);
  	
    // Print contents of stack
    // while (Stack.empty() == false) {
    //     cout << Stack.top() << " ";
    //     Stack.pop();
    // }
    // cout << '\n';
}

int main(int argc, char **argv){
	// cout << "Total: " << argc << " arguments fed!\n";
	// for(int i = 0; i < argc; i++){
	// 	cout << argv[i] << '\n';
	// }

	// Read in the arguments
	int k;
	istringstream iss(argv[2]);
	iss >> k;
	ASSERT(!iss.fail(), "Input k should be given!\nPlease follow the format below:\n./map -k 4 map01.blif output.blif");
	char* input_file_name = argv[3];
	char* output_file_name = argv[4];
	cout << "k: " << k << ", input file name: " << input_file_name << ", output file name: " << output_file_name << '\n';

	// Read in input file
	ifstream infile(input_file_name);

	// Handle model name, PI, PO
	// Use std::string to avoid error for declaration of string datatype
	std::string tmp_raw;
	std::string model_name_raw;
	std::string primary_input_raw;
	std::string primary_output_raw;
	char delimiter = ' ';
	getline(infile, model_name_raw);
	getline(infile, tmp_raw);
	if(tmp_raw[tmp_raw.size()-1] == '\\'){ // Multiple line
		while(tmp_raw[tmp_raw.size()-1] == '\\'){
			primary_input_raw.append(tmp_raw.begin(), tmp_raw.end()-1);
			getline(infile, tmp_raw);
		}
	}
	primary_input_raw.append(tmp_raw);
	getline(infile, tmp_raw);
	
	if(tmp_raw[tmp_raw.size()-1] == '\\'){
		while(tmp_raw[tmp_raw.size()-1] == '\\'){ // Multiple line
			primary_output_raw.append(tmp_raw.begin(), tmp_raw.end()-1);
			getline(infile, tmp_raw);
		}
	}
	primary_output_raw.append(tmp_raw);
	
	// getline(infile, primary_input_raw);
	// getline(infile, primary_output_raw);
	cout << "Model name: " << model_name_raw << ", PI: " << primary_input_raw << ", PO: " << primary_output_raw << "\n";

	// +1 for skipping white space of the delimiter
	std::string model_name = model_name_raw.substr(model_name_raw.find(delimiter,0)+1);
	std::vector<std::string> PIs;
	std::vector<std::string> POs;
	
	// Parse primary_input_raw for PIs with delimiter
	size_t last = 0;
	size_t next = 0;
	while((next = primary_input_raw.find(delimiter, last)) != std::string::npos){
		if(last != 0){ // Skip the first '.inputs'
			PIs.push_back(primary_input_raw.substr(last, next-last));
		}
		last = next + 1;
	}
	PIs.push_back(primary_input_raw.substr(last)); // For the element after the last delimiter

	// Parse primary_output_raw for POs with delimiter
	last = 0;
	next = 0;
	while((next = primary_output_raw.find(delimiter, last)) != std::string::npos){
		if(last != 0){ // Skip the first '.outputs'
			POs.push_back(primary_output_raw.substr(last, next-last));
		}
		last = next + 1;
	}
	POs.push_back(primary_output_raw.substr(last)); // For the element after the last delimiter

	// List of total NODEs
	std::vector<NODE> Total_NODES;
	// Mapping between node names and struct node
	std::map<std::string, NODE> NODES_map;

	// Initial nodes to do BFS (PI + CONSTANT)
	std::vector<std::string> BFS_seed;

	// cout << "Model name:\n" << model_name << '\n';

	// cout << "PIs: \n";
	for(std::vector<std::string>::const_iterator it = PIs.begin(); it != PIs.end(); ++it){
		// cout << *it << '\n';

		// NODE initialization
		NODE tmp_node;
		tmp_node.NODE_level = 0; // Level for PI is initialized as 0
		tmp_node.NODE_name = *it;
		tmp_node.NODE_type = "0"; // PI

		// PUT new NODE into total NODE list
		Total_NODES.push_back(tmp_node);
		
		// Put new NODE into mapp
		NODES_map[*it] = tmp_node;

		// Put PI into BFS_seed
		BFS_seed.push_back(*it);
	}

	// cout << "POs: \n";
	for(std::vector<std::string>::const_iterator it = POs.begin(); it != POs.end(); ++it){
		// cout << *it << '\n';

		// NODE initialization
		NODE tmp_node;
		tmp_node.NODE_level = -1; // Level for PO -1s initialized as -1 to tell if the NODE is traversed
		tmp_node.NODE_name = *it;
		tmp_node.NODE_type = "0"; // PO

		// PUT new NODE into total NODE list
		Total_NODES.push_back(tmp_node);
		
		// Put new NODE into NODE map
		NODES_map[*it] = tmp_node;
	}


	// for(std::map<std::string, NODE>::iterator it = NODES_map.begin(); it != NODES_map.end(); ++it){
	// 	cout << it->first << " -> (" << it->second.NODE_name << ", " << it->second.NODE_type << ")\n";
	// }
	// cout << "-----------------------------------------------\n";
	
	int count = 0; // Relation counts for every .names
	std::string node_type = "0"; // Record node type for each node
	std::string names_raw;
	std::string relations_raw;
	std::vector<std::vector<std::string> > tree_nodes;
	std::vector<std::string> names;
	getline(infile, names_raw);
	while(names_raw.find(".end") == std::string::npos){ // Not .end

		// cout << names_raw << '\n';

		// Parse names_raw for names of nodes with delimiter
		last = 0;
		next = 0;
		while((next = names_raw.find(delimiter, last)) != std::string::npos){
			if(last != 0){ // Skip the first element '.names'
				names.push_back(names_raw.substr(last, next-last)); // These are input nodes
			}
			last = next + 1;
		}
		names.push_back(names_raw.substr(last)); // For the element after the last delimiter, which is the output node

		// Parse gate type
		count = 0;
		node_type = "0"; // 1: AND, 2: OR, 3: INVERTER, 4: CONSTANT_0, 5: CONSTANT_1, 6: BUFFER, 7: UNKNOWN
		getline(infile, relations_raw);
		while(relations_raw.find(".names") == std::string::npos &&\
			  relations_raw.find(".end") == std::string::npos){ // Not .names and .end
			count++; // Record how many relation lines are read to parse constant_0, which has no relation line.
			// cout << count << ": " << relations_raw << '\n';
			for(int i = 0; i < names.size()-1; ++i){ // Go through every input
				if(relations_raw[i] == '-'){ // OR gate
					node_type = "2";
				}
				else if(relations_raw[i] == '0'){ // INVERTER gate
					node_type = "3";
					if(names.size() > 2) node_type = "7"; // INVERTER with more than 1 input
				}
				else{ // AND gate
					if(node_type != "0") continue; // Avoid changing determined node type
					node_type = "1";
					if(names.size() == 2) node_type = "6"; // BUFFER gate
				}
			}
			if(names.size() < 2) node_type = "5"; // CONSTANT_1 gate

			getline(infile, relations_raw);
		}
		// CONSTANT_0 gate
		if(count == 0) node_type = "4";

		// Insert node for circuit G
		// for(std::vector<std::string>::const_iterator itr = names.begin(); itr != names.end(); ++itr){
		// 	cout << *itr << ' ';
		// }
		// cout << '\n';

		// NODE creation
		int this_gate_idx = names.size()-1;
		// Create input NODES if needed
		for(int i = 0; i < this_gate_idx; ++i){
			if(NODES_map.count(names[i]) > 0){ // names[i] already exists.
				// Put this gate into its input nodes' output NODE list
				NODES_map[names[i]].OutNODEs.push_back(names[this_gate_idx]);
			}
			else{ // names[i] not exists.
				NODE tmp_node;
				tmp_node.NODE_level = -1;  // Initialize NODE level to -1 to tell if the NODE is traversed
				tmp_node.NODE_name = names[i];
				tmp_node.NODE_type = "0"; // Initialize its node type to PI/PO
				
				// PUT new NODE into total NODE list
				Total_NODES.push_back(tmp_node);
		
				// Put new NODE into NODE map
				NODES_map[names[i]] = tmp_node;

				// Add gate NODE to its input NODEs' output NODE list
				NODES_map[names[i]].OutNODEs.push_back(names[this_gate_idx]);
			}
		}

		// Create gate NODE if needed
		if(NODES_map.count(names[this_gate_idx]) > 0){ // names[-1] already exists.
			NODES_map[names[this_gate_idx]].NODE_type = node_type; // Update node type
			for(int i = 0; i < this_gate_idx; ++i){ // Update input NODE list
				NODES_map[names[this_gate_idx]].InNODEs.push_back(names[i]);
			}
		}
		else{
			NODE tmp_node;
			tmp_node.NODE_level = -1;  // Initialize NODE level to -1 to tell if the NODE is traversed
			tmp_node.NODE_name = names[this_gate_idx];
			if(node_type == "4" || node_type == "5"){
				tmp_node.NODE_level = 0; // Initialize NODE level to 0 if it is CONSTANT gate
				// Put CONSTANT into BFS_seed
				BFS_seed.push_back(names[this_gate_idx]);
			}
			tmp_node.NODE_type = node_type;
				
			// PUT new NODE into total NODE list
			Total_NODES.push_back(tmp_node);
		
			// Put new NODE into NODE map
			NODES_map[names[this_gate_idx]] = tmp_node;

			for(int i = 0; i < this_gate_idx; ++i){ // Update input NODE list
				NODES_map[names[this_gate_idx]].InNODEs.push_back(names[i]);
			}
		}

		names.push_back(node_type);
		tree_nodes.push_back(names);
		names.clear();

		names_raw = relations_raw;
	}

	// DFS stack for topological order
	std::stack<std::string> Labeling_Stack;
	topologicalSort(Labeling_Stack, NODES_map, Total_NODES);

	// Traverse the graph in topological order and assign FlowMap label (level) to each NODE.
	while (!Labeling_Stack.empty()){
		std::string NODE_front = Labeling_Stack.top();
        Labeling_Stack.pop();

        // Update NODE level for level undecided NODE
		if(NODES_map[NODE_front].NODE_level == -1){

			// Breaking down the nodes with more than two inputs
			while(NODES_map[NODE_front].InNODEs.size() > 2){
				int smallest = INT_MAX;
				int small_second = INT_MAX;
				int smallest_idx = -1;
				int small_second_idx = -1;
				NODE tmp_node;
				NODE smallest_NODE;
				NODE small_second_NODE;

				// For new NODE name
				stringstream i2s;
				i2s << NODES_map.size();
				std::string new_name = i2s.str();

				// Find the smallest two
				for(int k = 0; k < NODES_map[NODE_front].InNODEs.size(); ++k){
					if(NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level < smallest){
						if(smallest != small_second){
							small_second = smallest;
							small_second_idx = smallest_idx;
						}
						smallest = NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level;
						smallest_idx = k;
					}
					else{
						if(NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level != smallest || smallest != small_second){
							if(NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level < small_second){
								small_second = NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level;
								small_second_idx = k;
							}
						}
					}
				}
				
				// Create a new NODE
				new_name = "[" + new_name + "]";
				tmp_node.NODE_name = new_name;
				// Copy the original gate type to new NODE
				tmp_node.NODE_type = NODES_map[NODE_front].NODE_type;
				// Add the smallest two NODES into new NODE's input NODE list
				tmp_node.InNODEs.push_back(NODES_map[NODE_front].InNODEs[smallest_idx]);
				tmp_node.InNODEs.push_back(NODES_map[NODE_front].InNODEs[small_second_idx]);
				// Add the gate NODE into new NODE's output NODE list
				tmp_node.OutNODEs.push_back(NODE_front);

				smallest_NODE = NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]];
				small_second_NODE = NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]];

				// Use the max level of the two input NODES +1 as the level of new NODE
				tmp_node.NODE_level = ((smallest_NODE.NODE_level > small_second_NODE.NODE_level) ?
										smallest_NODE.NODE_level : small_second_NODE.NODE_level) + 1;

				// Delete gate NODE from the smallest two NODEs' output NODE list
				// CANNOT use new NODE variable (smallest_NODE, small_second_NODE) to adjust original output NODE lists
				NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.erase(
					std::remove(NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.begin(),
								NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.end(),
								NODE_front),
					NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.end());
				NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.erase(
					std::remove(NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.begin(),
								NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.end(),
								NODE_front),
					NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.end());

				// Add the new NODE to the smallest two NODEs' output NODE list
				NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.push_back(new_name);
				NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.push_back(new_name);

				// Delete the smallest two NODEs from gate NODE's input NODE list
				NODES_map[NODE_front].InNODEs.erase(NODES_map[NODE_front].InNODEs.begin() + smallest_idx);
				if(smallest_idx < small_second_idx){ // To prevent from mis-indexing caused by deleting from previous line
					NODES_map[NODE_front].InNODEs.erase(NODES_map[NODE_front].InNODEs.begin() + small_second_idx-1);
				}
				else{
					NODES_map[NODE_front].InNODEs.erase(NODES_map[NODE_front].InNODEs.begin() + small_second_idx);
				}

				// Add the new NODE to gate NODE input NODE list
				NODES_map[NODE_front].InNODEs.push_back(new_name);

				// PUT new NODE into total NODE list
				Total_NODES.push_back(tmp_node);
		
				// Put new NODE into NODE map
				NODES_map[new_name] = tmp_node;
			}

			// Update the front NODE level to the max of the two input NODEs
			int max_level = -1;
			for(int j = 0; j < NODES_map[NODE_front].InNODEs.size(); ++j){
				if(NODES_map[NODES_map[NODE_front].InNODEs[j]].NODE_level == -1){ // Input NODE's level info is not complete
					max_level = -2; // To make the front NODE's level be -1 after this round
					break;
				}
				if(NODES_map[NODES_map[NODE_front].InNODEs[j]].NODE_level > max_level){
					max_level = NODES_map[NODES_map[NODE_front].InNODEs[j]].NODE_level;
				}
			}
			NODES_map[NODE_front].NODE_level = max_level+1;
		}
    }

    // Get the topological order for newly decomposed newwork
    std::stack<std::string> FlowMap_cut_Stack;
	topologicalSort(FlowMap_cut_Stack, NODES_map, Total_NODES);

	// Create a map to record l(v) for a given node name
	std::map<std::string, int> FlowMap_NODE_label;
	// Initialize l(v)=0 for v in PI list
	for(std::vector<std::string>::const_iterator it = PIs.begin(); it != PIs.end(); ++it){
		FlowMap_NODE_label[*it] = 0;
	}

	// Save all the collapsed network N_t' as a map
	std::map<std::string, std::map<std::string, node> > N_t_prime_map;

    // FlowMap cut
    cout << "FlowMap cut sequence:\n";
    while (!FlowMap_cut_Stack.empty()){
		std::string NODE_front = FlowMap_cut_Stack.top();
        FlowMap_cut_Stack.pop();

        // Skip the PI nodes, CONSTANT as well
        if(NODES_map[NODE_front].InNODEs.size()==0) continue;

        cout << "Subgraph ending at " << NODE_front << ": \n";

        // Construct subnetwork N_t backtracking from NODE t (NODE_front)
		graph N_t;
		node tmp_node;
		edge tmp_edge;
		std::queue<std::string> N_t_traversal; // For BFS traversal
		std::map<std::string, node> N_t_node_mapping; // Recording N_t nodes

		// Create node t and map it with its name
		node t = N_t.new_node();
		N_t_node_mapping[NODE_front] = t;

		// Set the front NODE as BFS traversal root
		N_t_traversal.push(NODE_front);

		// Subnetwork N_t construction
		while(!N_t_traversal.empty()){
			std::string predecessor = N_t_traversal.front();
			N_t_traversal.pop();

			// Build edges to connect this node with all its descendant
			for(int i = 0; i < NODES_map[predecessor].InNODEs.size(); ++i){
				// NODES_map[predecessor].InNODEs[i] has already been visited
				if(N_t_node_mapping.count(NODES_map[predecessor].InNODEs[i]) > 0) continue;

				// NODES_map[predecessor].InNODEs[i] has not been visited
				N_t_traversal.push(NODES_map[predecessor].InNODEs[i]);

				// Create a new node for predecessor's input node that has not been visited
				tmp_node = N_t.new_node();

				// Map the new node with its name
				N_t_node_mapping[NODES_map[predecessor].InNODEs[i]] = tmp_node;

				// Create a edge to connect the new input node and predecessor node
				tmp_edge = N_t.new_edge(N_t_node_mapping[NODES_map[predecessor].InNODEs[i]], N_t_node_mapping[predecessor]);
				cout << NODES_map[predecessor].InNODEs[i] << " -> " << predecessor << '\n';
			}
		}

		// Let p = max{l(u): u in input(t)}
		int p = -1;
		for(int i = 0; i < NODES_map[NODE_front].InNODEs.size(); ++i){
			if(FlowMap_NODE_label[NODES_map[NODE_front].InNODEs[i]] > p){
				p = FlowMap_NODE_label[NODES_map[NODE_front].InNODEs[i]];
			}
		}
		FlowMap_NODE_label[NODE_front] = p;

		cout << "\np(" << NODE_front << "): " << p << '\n';

		// Create a new graph N_t_copied
		graph N_t_copied;
		N_t_copied = N_t; // Copy N_t graph and memory addresses of its nodes have also changed.

		// Copplapse all nodes in N_t with label p into t as new node t_prime
		std::vector<std::string> Collapsed_nodes;
		std::string t_prime_name = NODE_front+"_prime";
		node t_prime = N_t.new_node();
		for(std::map<std::string, node>::iterator it = N_t_node_mapping.begin(); it != N_t_node_mapping.end(); ++it){
			if(FlowMap_NODE_label[it->first] == p){
				Collapsed_nodes.push_back(it->first);
			}
		}
		for(std::vector<std::string>::iterator it = Collapsed_nodes.begin(); it != Collapsed_nodes.end(); ++it){
			// N_t node with label p
			if(FlowMap_NODE_label[*it] == p){
				cout << *it << "(";
				N_t.print_node(N_t_node_mapping[*it]);
				cout << ") has label: " << FlowMap_NODE_label[*it] << ", and p is: " << p << '\n';

				for(int i = 0; i < NODES_map[*it].InNODEs.size(); ++i){
					// If the input node of collapsed node is also need to be collapsed together with t
					if(FlowMap_NODE_label[NODES_map[*it].InNODEs[i]] == p) continue;

					// Connect input nodes of collapsed node to t_prime
					N_t.new_edge(N_t_node_mapping[NODES_map[*it].InNODEs[i]], t_prime);
				}

				// Delete the collapsed node
				N_t.hide_node(N_t_node_mapping[*it]);
				N_t_node_mapping.erase(*it);
			}
		}
		N_t_node_mapping[t_prime_name] = t_prime;

		node print_node;
		edge print_edge;
		cout << "Collapsed network:\n";
		cout << "\tNodes: ";
		forall_nodes(print_node, N_t)
			N_t.print_node(print_node);
		cout << "\n\tEdges:\n";
		forall_edges(print_edge, N_t){
			N_t.print_edge(print_edge);
			cout << '\n';
		}

		cout << "Copied original network:\n";
		cout << "\tNodes: ";
		forall_nodes(print_node, N_t_copied)
			N_t_copied.print_node(print_node);
		cout << "\n\tEdges:\n";
		forall_edges(print_edge, N_t_copied){
			N_t_copied.print_edge(print_edge);
			cout << '\n';
		}
		cout << "-----------------------------------------------\n";

		// edge e;
		// forall_edges(e, N_t)
		// 	N_t.print_edge(e);

		// node n0 = G.new_node();
		// node n1 = G.new_node();
		// node n2 = G.new_node();
		// node n3 = G.new_node();
		// edge e0 = G.new_edge(n0, n1);
		// edge e1 = G.new_edge(n1, n3);
		// edge e2 = G.new_edge(n0, n2);
		// edge e3 = G.new_edge(n2, n3);
		// edge e4 = G.new_edge(n0, n3);

		// node tmp_node;
		// std::vector<node> node_list;
		// for(int i = 0; i < 10; ++i){
		// 	tmp_node = G.new_node();
		// 	node_list.push_back(tmp_node);
		// }

		// edge tmp_edge;
		// for(int i = 0; i < 9; ++i){
		// 	tmp_edge = G.new_edge(node_list[i], node_list[node_list.size()-1]);
		// }

		// Assign edge weights
		// edge_array<int> weight(G);
		// weight[e0] = 1;
		// weight[e1] = 3;
		// weight[e2] = 2;
		// weight[e3] = 2;
		// weight[e4] = 5;


		// for(int i = 0; i < 10; ++i){
		// 	cout << "Node [" << i << "]:\n";
		// 	cout << "Address: " << node_list[i] << ", Value: ";
		// 	G.print_node(node_list[i]);
		// 	cout << ", Degree: " << G.indeg(node_list[i]);
		// 	cout << '\n';
		// }

		// cout << "Degree [n0]: " << G.indeg(n0);
		// cout << '\n';
		// cout << "Degree [n1]: " << G.indeg(n1);
		// cout << '\n';
		// cout << "Degree [n2]: " << G.indeg(n2);
		// cout << '\n';
		// cout << "Degree [n3]: " << G.indeg(n3);
		// cout << '\n';

		// forall_in_edges(tmp_edge, node_list[node_list.size()-1])
		// 	weight[tmp_edge] = 2;

		// int total_weight = 0;
		// forall_in_edges(tmp_edge, node_list[node_list.size()-1])
		// 	total_weight += weight[tmp_edge];
		// cout << total_weight << '\n';

		// const list<node> All_nodes = G.all_nodes();
		// node n;
		// forall(n, All_nodes)
		// 	G.print_node(n);
		// cout << endl;


		// Min cut algorithm
		// list<node> cut;
		// int cut_value = MIN_CUT(G, weight, cut);

		// cout << "The minimum cut has value: " << cut_value << endl;
		// cout << "cut:";
		// node v;
		// forall(v, cut)
		// 	G.print_node(v);
		// cout << endl;

		// G.del_all_nodes();
    }

	// // BFS Traversal to go through the graph in topological order
	// int BFS_count = 0; // To check how many unique nodes that are traversed
	// std::queue<std::string> NODE_queue;
	// // Starting BFS from BFS_seed list
	// for(std::vector<std::string>::const_iterator it = BFS_seed.begin(); it != BFS_seed.end(); ++it){
	// 	NODE_queue.push(*it);
	// 	BFS_count++;
	// }

	// while(!NODE_queue.empty()){
	// 	std::string NODE_front = NODE_queue.front();
	// 	NODE_queue.pop();

	// 	// Add output NODEs into the queue
	// 	for(int i = 0; i < NODES_map[NODE_front].OutNODEs.size(); ++i){
	// 		bool add_output_list = true;
	// 		for(int j = 0; j < NODES_map[NODES_map[NODE_front].OutNODEs[i]].InNODEs.size(); ++j){
	// 			// Check if output NODE's input NODEs have complete level info: for topological order with BFS
	// 			if(NODES_map[NODES_map[NODES_map[NODE_front].OutNODEs[i]].InNODEs[j]].NODE_level == -1){
	// 				add_output_list = false;
	// 				break;
	// 			}
	// 		}
	// 		if(add_output_list){ // Add the output NODE only if its input NODEs have complete level info.
	// 			NODE_queue.push(NODES_map[NODE_front].OutNODEs[i]);
	// 		}
	// 		else{ // Add the front NODE back to the queue for later ouput NODEs traversal
	// 			NODE_queue.push(NODE_front);
	// 		}
	// 	}

	// 	// Update NODE level for level undecided NODE
	// 	if(NODES_map[NODE_front].NODE_level == -1){
	// 		BFS_count++;

	// 		// Breaking down the nodes with more than two inputs
	// 		while(NODES_map[NODE_front].InNODEs.size() > 2){
	// 			int smallest = INT_MAX;
	// 			int small_second = INT_MAX;
	// 			int smallest_idx = -1;
	// 			int small_second_idx = -1;
	// 			NODE tmp_node;
	// 			NODE smallest_NODE;
	// 			NODE small_second_NODE;

	// 			// For new NODE name
	// 			stringstream i2s;
	// 			i2s << NODES_map.size();
	// 			std::string new_name = i2s.str();

	// 			// Find the smallest two
	// 			for(int k = 0; k < NODES_map[NODE_front].InNODEs.size(); ++k){
	// 				if(NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level < smallest){
	// 					if(smallest != small_second){
	// 						small_second = smallest;
	// 						small_second_idx = smallest_idx;
	// 					}
	// 					smallest = NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level;
	// 					smallest_idx = k;
	// 				}
	// 				else{
	// 					if(NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level != smallest || smallest != small_second){
	// 						if(NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level < small_second){
	// 							small_second = NODES_map[NODES_map[NODE_front].InNODEs[k]].NODE_level;
	// 							small_second_idx = k;
	// 						}
	// 					}
	// 				}
	// 			}
				
	// 			// Create a new NODE
	// 			new_name = "[" + new_name + "]";
	// 			tmp_node.NODE_name = new_name;
	// 			// Copy the original gate type to new NODE
	// 			tmp_node.NODE_type = NODES_map[NODE_front].NODE_type;
	// 			// Add the smallest two NODES into new NODE's input NODE list
	// 			tmp_node.InNODEs.push_back(NODES_map[NODE_front].InNODEs[smallest_idx]);
	// 			tmp_node.InNODEs.push_back(NODES_map[NODE_front].InNODEs[small_second_idx]);
	// 			// Add the gate NODE into new NODE's output NODE list
	// 			tmp_node.OutNODEs.push_back(NODE_front);

	// 			smallest_NODE = NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]];
	// 			small_second_NODE = NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]];

	// 			// Use the max level of the two input NODES +1 as the level of new NODE
	// 			tmp_node.NODE_level = ((smallest_NODE.NODE_level > small_second_NODE.NODE_level) ?
	// 									smallest_NODE.NODE_level : small_second_NODE.NODE_level) + 1;

	// 			// Delete gate NODE from the smallest two NODEs' output NODE list
	// 			// CANNOT use new NODE variable (smallest_NODE, small_second_NODE) to adjust original output NODE lists
	// 			NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.erase(
	// 				std::remove(NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.begin(),
	// 							NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.end(),
	// 							NODE_front),
	// 				NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.end());
	// 			NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.erase(
	// 				std::remove(NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.begin(),
	// 							NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.end(),
	// 							NODE_front),
	// 				NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.end());

	// 			// Add the new NODE to the smallest two NODEs' output NODE list
	// 			NODES_map[NODES_map[NODE_front].InNODEs[smallest_idx]].OutNODEs.push_back(new_name);
	// 			NODES_map[NODES_map[NODE_front].InNODEs[small_second_idx]].OutNODEs.push_back(new_name);

	// 			// Delete the smallest two NODEs from gate NODE's input NODE list
	// 			NODES_map[NODE_front].InNODEs.erase(NODES_map[NODE_front].InNODEs.begin() + smallest_idx);
	// 			if(smallest_idx < small_second_idx){ // To prevent from mis-indexing caused by deleting from previous line
	// 				NODES_map[NODE_front].InNODEs.erase(NODES_map[NODE_front].InNODEs.begin() + small_second_idx-1);
	// 			}
	// 			else{
	// 				NODES_map[NODE_front].InNODEs.erase(NODES_map[NODE_front].InNODEs.begin() + small_second_idx);
	// 			}

	// 			// Add the new NODE to gate NODE input NODE list
	// 			NODES_map[NODE_front].InNODEs.push_back(new_name);

	// 			// PUT new NODE into total NODE list
	// 			Total_NODES.push_back(tmp_node);
		
	// 			// Put new NODE into NODE map
	// 			NODES_map[new_name] = tmp_node;
	// 		}

	// 		// Get the level for each node in BFS traversal without transform to two inputs graph
	// 		int max_level = -1;
	// 		for(int j = 0; j < NODES_map[NODE_front].InNODEs.size(); ++j){
	// 			if(NODES_map[NODES_map[NODE_front].InNODEs[j]].NODE_level == -1){ // Input NODE's level info is not complete
	// 				NODE_queue.push(NODE_front); // Re-add front NODE to BFS queue
	// 				max_level = -2; // To make the front NODE's level be -1 after this round
	// 				break;
	// 			}
	// 			if(NODES_map[NODES_map[NODE_front].InNODEs[j]].NODE_level > max_level){
	// 				max_level = NODES_map[NODES_map[NODE_front].InNODEs[j]].NODE_level;
	// 			}
	// 		}
	// 		NODES_map[NODE_front].NODE_level = max_level+1;
	// 	}
	// }

	// To check how many unique nodes that are traversed
	// cout << "BFS count: " << BFS_count << ", Total nodes count: " << NODES_map.size() << '\n';

	// for(std::map<std::string, NODE>::iterator it = NODES_map.begin(); it != NODES_map.end(); ++it){
	// 	cout << "Index: " << it->first << " -> NODE: (" << it->second.NODE_name << ", ";
	// 	if(it->second.NODE_level == -1) cout << "Level: UNDECIDED, ";
	// 	else cout << "Level: " << it->second.NODE_level << ", ";
	// 	if(it->second.NODE_type == "1") cout << "AND";
	// 	else if(it->second.NODE_type == "2") cout << "OR";
	// 	else if(it->second.NODE_type == "3") cout << "INVERTER";
	// 	else if(it->second.NODE_type == "4") cout << "CONSTANT_0";
	// 	else if(it->second.NODE_type == "5") cout << "CONSTANT_1";
	// 	else if(it->second.NODE_type == "6") cout << "BUFFER";
	// 	else if(it->second.NODE_type == "7") cout << "UNKNOWN";
	// 	else cout << "PI/PO";
	// 	cout << ")\n";
	// 	cout << "Input NODEs (" << it->second.InNODEs.size() << "): \n";
	// 	for(int i = 0; i < it->second.InNODEs.size(); ++i){
	// 		cout << it->second.InNODEs[i] << ' ';
	// 	}
	// 	cout << '\n';

	// 	cout << "Output NODEs (" << it->second.OutNODEs.size() << "): \n";
	// 	for(int i = 0; i < it->second.OutNODEs.size(); ++i){
	// 		cout << it->second.OutNODEs[i] << ' ';
	// 	}
	// 	cout << '\n';
	// 	cout << "-----------------------------------------------\n";
	// }

	return 0;
}

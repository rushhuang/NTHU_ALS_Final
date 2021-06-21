// LEDA Manual
// http://www.algorithmic-solutions.info/leda_manual/MANUAL.html

#include <iostream>
#include <string>
#include <sstream> // convert sting to int and vice versa
#include <fstream>
#include <vector>
#include <map>
#include <climits> // INT_MAX
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

	std::string node_name;

	// node_type
	// 0: PI, 1: AND, 2: OR, 3: INVERTER, 4: CONSTANT_0, 5: CONSTANT_1, 6: BUFFER, 7: UNKNOWN
	std::string node_type;

	// Input nodes
	std::vector<std::string> InNODEs;

	// Output nodes
	std::vector<std::string> OutNODEs;
};

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
	std::string model_name_raw;
	std::string primary_input_raw;
	std::string primary_output_raw;
	char delimiter = ' ';
	getline(infile, model_name_raw);
	getline(infile, primary_input_raw);
	getline(infile, primary_output_raw);
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

	// cout << "Model name:\n" << model_name << '\n';

	// cout << "PIs: \n";
	for(std::vector<std::string>::const_iterator it = PIs.begin(); it != PIs.end(); ++it){
		// cout << *it << '\n';

		// NODE initialization
		NODE tmp_node;
		tmp_node.node_name = *it;
		tmp_node.node_type = "0"; // PI

		// PUT new NODE into total NODE list
		Total_NODES.push_back(tmp_node);
		
		// Put new NODE into mapp
		NODES_map[*it] = tmp_node;
	}

	// cout << "POs: \n";
	for(std::vector<std::string>::const_iterator it = POs.begin(); it != POs.end(); ++it){
		// cout << *it << '\n';

		// NODE initialization
		NODE tmp_node;
		tmp_node.node_name = *it;
		tmp_node.node_type = "0"; // PO

		// PUT new NODE into total NODE list
		Total_NODES.push_back(tmp_node);
		
		// Put new NODE into NODE map
		NODES_map[*it] = tmp_node;
	}


	// for(std::map<std::string, NODE>::iterator it = NODES_map.begin(); it != NODES_map.end(); ++it){
	// 	cout << it->first << " -> (" << it->second.node_name << ", " << it->second.node_type << ")\n";
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

		int this_gate_idx = names.size()-1;
		// Create input NODES if needed
		for(int i = 0; i < this_gate_idx; ++i){
			if(NODES_map.count(names[i]) > 0){ // names[i] already exists.
				// Put this gate into its input nodes' output NODE list
				NODES_map[names[i]].OutNODEs.push_back(names[this_gate_idx]);
			}
			else{ // names[i] not exists.
				NODE tmp_node;
				tmp_node.node_name = names[i];
				tmp_node.node_type = "0"; // Initialize its node type to PI/PO
				
				// PUT new NODE into total NODE list
				Total_NODES.push_back(tmp_node);
		
				// Put new NODE into NODE map
				NODES_map[names[i]] = tmp_node;

				tmp_node.OutNODEs.push_back(names[this_gate_idx]);
			}
		}

		// Create gate NODE if needed
		if(NODES_map.count(names[this_gate_idx]) > 0){ // names[-1] already exists.
			NODES_map[names[this_gate_idx]].node_type = node_type; // Update node type
			for(int i = 0; i < this_gate_idx; ++i){ // Update input NODE list
				NODES_map[names[this_gate_idx]].InNODEs.push_back(names[i]);
			}
		}
		else{
			NODE tmp_node;
			tmp_node.node_name = names[this_gate_idx];
			tmp_node.node_type = node_type;
				
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

	// cout << "Node names: \n";
	// for(std::vector<std::vector<std::string> >::const_iterator it = tree_nodes.begin(); it != tree_nodes.end(); ++it){
	// 	for(std::vector<std::string>::const_iterator itr = it->begin(); itr != it->end(); ++itr){
	// 		if(*itr == "1") cout << "AND ";
	// 		else if(*itr == "2") cout << "OR ";
	// 		else if(*itr == "3") cout << "INV ";
	// 		else if(*itr == "4") cout << "CONSTANT_0";
	// 		else if(*itr == "5") cout << "CONSTANT_1";
	// 		else if(*itr == "6") cout << "BUFFER";
	// 		else if(*itr == "7") cout << "UNKNOWN ";
	// 		else cout << *itr << ' ';
	// 	}
	// 	cout << '\n'; 
	// }

	for(std::map<std::string, NODE>::iterator it = NODES_map.begin(); it != NODES_map.end(); ++it){
		cout << it->first << " -> (" << it->second.node_name << ", ";
		if(it->second.node_type == "1") cout << "AND";
		else if(it->second.node_type == "2") cout << "OR";
		else if(it->second.node_type == "3") cout << "INVERTER";
		else if(it->second.node_type == "4") cout << "CONSTANT_0";
		else if(it->second.node_type == "5") cout << "CONSTANT_1";
		else if(it->second.node_type == "6") cout << "BUFFER";
		else if(it->second.node_type == "7") cout << "UNKNOWN";
		else cout << "PI/PO";
		cout << ")\n";
		cout << "Input NODEs (" << it->second.InNODEs.size() << "): \n";
		for(int i = 0; i < it->second.InNODEs.size(); ++i){
			cout << it->second.InNODEs[i] << ' ';
		}
		cout << '\n';

		cout << "Output NODEs (" << it->second.OutNODEs.size() << "): \n";
		for(int i = 0; i < it->second.OutNODEs.size(); ++i){
			cout << it->second.OutNODEs[i] << ' ';
		}
		cout << '\n';
		cout << "-----------------------------------------------\n";
	}



	// Graph construction
	// graph G;
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
	return 0;
}

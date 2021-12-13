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
#include <cmath> // For pow()
#include <bitset> // To generate test pattern
using namespace std;

#include <LEDA/core/list.h>
#include <LEDA/graph/graph.h>
#include <LEDA/graph/min_cut.h>
#include <LEDA/graph/max_flow.h>
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

int dfs(std::string node, std::map<std::string, bool>& visited, std::map<std::string, int>& level_map, std::map<std::string, NODE>& mapping){
  	if(mapping[node].InNODEs.size() == 0){
    	visited[node] = true; // mark this unvisited node as visited

  		level_map[node] = 0;
  		return 0; // The bottom level is set to 0
  	}

  	int max_level = -1;
    for(std::vector<std::string>::iterator it = mapping[node].InNODEs.begin(); it != mapping[node].InNODEs.end(); ++it){
    	int this_level = 0;
        // if the node is not visited
        if (visited.count(*it) == 0) {
            // call the dfs function
            this_level = dfs(*it, visited, level_map, mapping)+1;
        }
        else{
        	this_level = level_map[*it]+1;
        }
        if(this_level > max_level) max_level = this_level;
    }

    visited[node] = true; // mark this unvisited node as visited
    
    level_map[node] = max_level;
    return max_level;
}

// int k2dfs(std::string node, std::map<std::string, bool>& visited, std::map<std::string, NODE>& mapping){
//     // mark every node as visited
//     visited[node] = true;
  	
//   	if(mapping[node].InNODEs.size() == 0) return 0; // The bottom level set to 0

//   	int max_level = 0;
//     for(std::vector<std::string>::iterator it = mapping[node].InNODEs.begin(); it != mapping[node].InNODEs.end(); ++it){
//         // if the node is not visited
//         if (visited.count(*it) == 0) {
//             // call the dfs function
//             int this_level = k2dfs(*it, visited, mapping)+1;
//             if(this_level > max_level) max_level = this_level;
//         }
//     }
//     return max_level;
// }

int klutdfs(std::string node, std::map<std::string, bool>& visited, std::map<std::string, int>& node_output_every_round,
			std::string test_pattern, std::vector<std::string>& input_signals, std::map<std::string, NODE>& mapping){
    
  	int signal = (mapping[node].NODE_type == "1")?1:0; // Set initial value to the non-controlling value of the gate
  	int input_signal_error = 0; // if input signal is conflict with CONSTANT
  	int out_signal;

  	// Check if this node is input signal
  	for(int i = 0; i < input_signals.size(); ++i){
  		if(input_signals[i] == node){
  			istringstream iss( std::string(1,test_pattern[i]) ); // Convert char to string first then to int
			iss >> signal;

			// cout << "Testing input node (" << node << ") with signal: " << signal << '\n';
  			if(mapping[node].NODE_type == "4"){ // CONSTANT 0
  				// signal = 0;
  				if(!signal) return signal;
  				else input_signal_error = 1;
  			}
  			else if(mapping[node].NODE_type == "5"){ // CONSTANT 1
  				// signal = 1;
  				if(signal) return signal;
  				else input_signal_error = 1;
  			}
  			else{
				// istringstream iss( std::string(1,test_pattern[i]) ); // Convert char to string first then to int
				// iss >> signal;

				visited[node] = true; //  mark this unvisited node as visited
				node_output_every_round[node] = signal;
	  			return signal; // Set this input as its test_pattern signal
  			}
  		}
  	}
  	if(input_signal_error) return -1;

  	// cout << "Input nodes for node " << node << ": ";
    for(std::vector<std::string>::iterator it = mapping[node].InNODEs.begin(); it != mapping[node].InNODEs.end(); ++it){
    	// cout << *it << ' ';
        // if the predecessor is not visited
        if (visited.count(*it) == 0) {
        	out_signal = klutdfs(*it, visited, node_output_every_round, test_pattern, input_signals, mapping);
        }
        else{ // Read from previosly computed output of the predecessor
        	out_signal = node_output_every_round[*it];
        }
        if(out_signal == -1) return -1; // Propagate the input signal error
        if(mapping[node].NODE_type == "1"){ // If this node is AND
        	signal = signal & out_signal;
        }
        else if(mapping[node].NODE_type == "2"){ // If this node is OR
        	signal = signal | out_signal;
        }
        else if(mapping[node].NODE_type == "3"){ // If this node is INVERTER
        	signal = 1 ^ out_signal;
        }
        else if(mapping[node].NODE_type == "6"){ // If this node is BUFFER
        	signal = out_signal;
		}
		else ASSERT(true, "UNKNOWN gate type!");
        // node_output_every_round[*it] = signal; // If the predecessor is not visited, add the output value to it
    }
    // cout << '\n';
    // cout << "Output signal for node (" << node << ") is " << signal << '\n';
    node_output_every_round[node] = signal;
    visited[node] = true; //  mark this unvisited node as visited
    return signal;
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
	// cout << "k: " << k << ", input file name: " << input_file_name << ", output file name: " << output_file_name << '\n';

	cout << "\rParsing inputs..." << std::flush;

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
	// cout << "Model name: " << model_name_raw << ", PI: " << primary_input_raw << ", PO: " << primary_output_raw << "\n";

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

	cout << "\rDecomposing..." << std::flush;

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

     // 2-LUT mapping is equivalent to network decomposition
    cout << "\rOutputing..." << std::flush;
    int LUT_count = 0;

	// Output blif file
	ofstream outblif;
	outblif.open(output_file_name);

	outblif << model_name_raw << '\n';
	outblif << primary_input_raw << '\n';
	outblif << primary_output_raw << '\n';

	// cout << model_name_raw << '\n';
	// cout << primary_input_raw << '\n';
	// cout << primary_output_raw << '\n';

	std::queue<std::string> Output_traversal;
	std::map<std::string, bool> Output_visited_NODEs;
	// cout << "PO: ";
	for(int i = 0; i < POs.size(); ++i){
		Output_traversal.push(POs[i]);
		// cout << POs[i] << ' ';
	}
	// cout << '\n';

	while(!Output_traversal.empty()){
		std::string this_node = Output_traversal.front();
		Output_traversal.pop();

		// Skip NODEs that have been added multiple times before it is really traversed
		if(Output_visited_NODEs.count(this_node) != 0) continue;

		Output_visited_NODEs[this_node] = true;

		if(NODES_map[this_node].NODE_type == "0") continue; // PI does not need to be printed as .names

		outblif << ".names ";
		LUT_count++;
		// cout << ".names ";
		for(int i = 0; i < NODES_map[this_node].InNODEs.size(); ++i){
			// Only the unvisited NODEs need to be put in the queue
			if(Output_visited_NODEs.count(NODES_map[this_node].InNODEs[i]) == 0)
				Output_traversal.push(NODES_map[this_node].InNODEs[i]);

			outblif << NODES_map[this_node].InNODEs[i] << ' ';
			// cout << NODES_map[this_node].InNODEs[i] << ' ';
		}
		outblif << this_node << '\n';
		// cout << this_node << '\n';

		if(NODES_map[this_node].NODE_type == "1"){
			// cout << "AND";
			// For input NODEs
			for(int i = 0; i < NODES_map[this_node].InNODEs.size(); ++i){
				outblif << "1";
				// cout << "1";
			}
			outblif << ' ';
			// cout << ' ';

			// For gate NODE
			outblif << "1\n";
			// cout << "1\n";
	    }
		else if(NODES_map[this_node].NODE_type == "2"){
			// cout << "OR";
			for(int i = 0; i < NODES_map[this_node].InNODEs.size(); ++i){
				// For input NODEs
				for(int j = 0; j < NODES_map[this_node].InNODEs.size(); ++j){
					if(i == j){
						outblif << "1";
						// cout << "1";
					}
					else{
						outblif << "-";
						// cout << "-";
					}
				}
				// For output NODEs
				outblif << " 1\n";
				// cout << " 1\n";
			}
		}
		else if(NODES_map[this_node].NODE_type == "3"){
			// cout << "INVERTER";
			outblif << "0 1\n";
			// cout << "0 1\n";
		}
		else if(NODES_map[this_node].NODE_type == "4"){
			// cout << "CONSTANT_0";
			continue;
		}
		else if(NODES_map[this_node].NODE_type == "5"){
			// cout << "CONSTANT_1";
			outblif << "1\n";
			// cout << "1\n";
		}
		else if(NODES_map[this_node].NODE_type == "6"){
			// cout << "BUFFER";
			outblif << "1 1\n";
			// cout << "1 1\n";
		}
		else if(NODES_map[this_node].NODE_type == "7"){
			// cout << "UNKNOWN\n";
		}
	}

	outblif << ".end";
	// cout << ".end\n";
	outblif.close();

	if(k == 2){
	    // DFS from PO to get the level of final graph
    	std::map<std::string, bool> visited;
    	std::map<std::string, int> level_map;
    	int circuit_level = 0;
    	for(int i = 0; i < POs.size(); ++i){
    		int level = dfs(POs[i], visited, level_map, NODES_map);
    		if(level > circuit_level) circuit_level = level;
    	}

	    cout << "\rThe circuit level is " << circuit_level << ".\n";
    	cout << "The number of LUTs is " << LUT_count << ".\n";

    	return 0;
	}

	cout << "\rFlowMap......" << std::flush;

    // Get the topological order for newly decomposed newwork
    std::stack<std::string> FlowMap_cut_Stack;
	topologicalSort(FlowMap_cut_Stack, NODES_map, Total_NODES);

	// Create a map to record l(v) for a given node name
	std::map<std::string, int> FlowMap_NODE_label;
	// Initialize l(v)=0 for v in PI list
	for(std::vector<std::string>::const_iterator it = PIs.begin(); it != PIs.end(); ++it){
		FlowMap_NODE_label[*it] = 0;
	}

	// Save input signals for mincut of N_t'' as a map
	std::map<std::string, std::vector<std::string> > t_KLUT_input_map;

    // FlowMap cut
    // cout << "FlowMap cut sequence:\n";
    while (!FlowMap_cut_Stack.empty()){
		std::string NODE_front = FlowMap_cut_Stack.top();
        FlowMap_cut_Stack.pop();

        // Skip the PI nodes, CONSTANT as well
        // if(NODES_map[NODE_front].InNODEs.size()==0) continue;
        if(std::find(PIs.begin(), PIs.end(), NODE_front) != PIs.end()) continue;

        // cout << "Subgraph ending at " << NODE_front << ": \n";

        // Construct subnetwork N_t backtracking from NODE t (NODE_front)
		graph N_t;
		node tmp_node;
		edge tmp_edge;
		std::stack<std::string> N_t_traversal; // For DFS traversal
		std::map<std::string, bool> visited_N_t_traversal;
		std::map<std::string, node> N_t_node_mapping; // Recording N_t nodes
		std::map<node, std::string> Rev_N_t_node_mapping; // Reverse N_t nodes mapping, which use ndoe as key

		// Create node t and map it with its name
		node t = N_t.new_node();
		N_t_node_mapping[NODE_front] = t;

		// Set the front NODE as DFS traversal root
		N_t_traversal.push(NODE_front);

		// Subnetwork N_t construction
		while(!N_t_traversal.empty()){
			std::string predecessor = N_t_traversal.top();
			N_t_traversal.pop();

			ASSERT(NODES_map[predecessor].InNODEs.size() <= 2, "Input size more than 2!!!");

			// Build edges to connect this node with all its descendant
			for(int i = 0; i < NODES_map[predecessor].InNODEs.size(); ++i){
				if(visited_N_t_traversal.count(NODES_map[predecessor].InNODEs[i]) == 0){
					// If this input NODE has not been visited
					N_t_traversal.push(NODES_map[predecessor].InNODEs[i]);

					if(N_t_node_mapping.count(NODES_map[predecessor].InNODEs[i]) == 0){
						// Create a new node for predecessor's input node that has not been created
						tmp_node = N_t.new_node();
						// Map the new node with its name
						N_t_node_mapping[NODES_map[predecessor].InNODEs[i]] = tmp_node;
					}

					// // Create a edge to connect the new input node and predecessor node
					// tmp_edge = N_t.new_edge(N_t_node_mapping[NODES_map[predecessor].InNODEs[i]], N_t_node_mapping[predecessor]);
					// // cout << NODES_map[predecessor].InNODEs[i] << " -> " << predecessor << '\n';
				}
				// Create a edge to connect the new input node and predecessor node
				tmp_edge = N_t.new_edge(N_t_node_mapping[NODES_map[predecessor].InNODEs[i]], N_t_node_mapping[predecessor]);
				// cout << NODES_map[predecessor].InNODEs[i] << " -> " << predecessor << '\n';
			}
			visited_N_t_traversal[predecessor] = true;
		}

		// Let p = max{l(u): u in input(t)}
		int p = -1;
		for(int i = 0; i < NODES_map[NODE_front].InNODEs.size(); ++i){
			if(FlowMap_NODE_label[NODES_map[NODE_front].InNODEs[i]] > p){
				p = FlowMap_NODE_label[NODES_map[NODE_front].InNODEs[i]];
			}
		}
		FlowMap_NODE_label[NODE_front] = p; // Temporarily assign p as the front NODE's label value

		// Add source node
		// cout << "\nLinking source with PI...\n";
		node print_node;
		edge print_edge;
		node source_node = N_t.new_node();
		N_t_node_mapping["source"] = source_node;
		FlowMap_NODE_label["source"] = 0;
		forall_nodes(print_node, N_t){			
			for(std::vector<std::string>::const_iterator it = PIs.begin(); it != PIs.end(); ++it){
				if(N_t_node_mapping.count(*it) != 0 && print_node == N_t_node_mapping[*it]){ // Use short-circuit evaluation to avoid building null mapping for every PI
					N_t.new_edge(source_node, print_node);
					// N_t.print_node(source_node);
					// cout << " -> ";
					// N_t.print_node(print_node);
					// cout << '\n';
				}
			}
		}

		// cout << "\np(" << NODE_front << "): " << p << '\n';

		// Create a new graph N_t_copied
		graph N_t_copied; // Representing original N_t
		// Copy N_t graph and memory addresses of its nodes have also changed.
		N_t_copied = N_t; // Use N_t as N_t_prime

		// Copplapse all nodes in N_t with label p into t as new node t_prime
		std::vector<std::string> Collapsed_nodes;
		std::vector<std::string> Input_node_list_of_t_prime;
		std::vector<std::string>::iterator iter; // To find the element for vector
		std::string t_prime_name = NODE_front+"_prime";
		node t_prime = N_t.new_node();
		for(std::map<std::string, node>::iterator it = N_t_node_mapping.begin(); it != N_t_node_mapping.end(); ++it){
			if(FlowMap_NODE_label[it->first] == p && it->first != "source"){ // Collapse nodes with label p except for source node
				Collapsed_nodes.push_back(it->first);
			}
		}
		// cout << "Nodes to be collapsed together:\n";
		// for(std::vector<std::string>::iterator it = Collapsed_nodes.begin(); it != Collapsed_nodes.end(); ++it){
		// 	cout << *it << ' ';
		// }
		// cout << '\n';
		for(std::vector<std::string>::iterator it = Collapsed_nodes.begin(); it != Collapsed_nodes.end(); ++it){
			// N_t node with label p
			if(FlowMap_NODE_label[*it] == p){
				// cout << *it << "(";
				// N_t.print_node(N_t_node_mapping[*it]);
				// cout << ") has label: " << FlowMap_NODE_label[*it] << ", and p is: " << p << '\n';

				for(int i = 0; i < NODES_map[*it].InNODEs.size(); ++i){
					// If the input node of collapsed node is also need to be collapsed together with t
					if(FlowMap_NODE_label[NODES_map[*it].InNODEs[i]] == p) continue;

					// Connect input nodes of collapsed node to t_prime
					iter = std::find(Input_node_list_of_t_prime.begin(),
									 Input_node_list_of_t_prime.end(), NODES_map[*it].InNODEs[i]);
					if(iter == Input_node_list_of_t_prime.end()){
						// If input node has not been added into input node list of t_prime
						N_t.new_edge(N_t_node_mapping[NODES_map[*it].InNODEs[i]], t_prime);
						Input_node_list_of_t_prime.push_back(NODES_map[*it].InNODEs[i]);
					}
				}
				if(N_t.indeg(N_t_node_mapping[*it]) == 1 && NODES_map[*it].InNODEs.size() == 0){ // Collapsed node is PIs
					if(std::find(Input_node_list_of_t_prime.begin(),
								 Input_node_list_of_t_prime.end(), "source") == Input_node_list_of_t_prime.end()){
						// If source has not been linked to t_prime
						N_t.new_edge(source_node, t_prime);
					}
				}

				// Delete the collapsed node
				N_t.del_node(N_t_node_mapping[*it]);
				N_t_node_mapping.erase(*it);
			}
		}
		N_t_node_mapping[t_prime_name] = t_prime;

		// node print_node;
		// edge print_edge;

		// Add source node
		// cout << "\nLinking source with PI...\n";
		// node source_node = N_t.new_node();
		// N_t_node_mapping["source"] = source_node;
		// FlowMap_NODE_label["source"] = 0;
		// forall_nodes(print_node, N_t){
		// 	if(N_t.indeg(print_node) == 0 && print_node != source_node){ // Link source node to every PI of N_t
		// 		N_t.new_edge(source_node, print_node);
		// 		// N_t.print_node(source_node);
		// 		// cout << " -> ";
		// 		// N_t.print_node(print_node);
		// 		// cout << '\n';
		// 	}
		// }

		// Build reverse N_t node mapping using current nodes
		for(std::map<std::string, node>::iterator iter = N_t_node_mapping.begin(); iter != N_t_node_mapping.end(); ++iter){
			Rev_N_t_node_mapping[iter->second] = iter->first;
		}

		// cout << "\nCollapsed network:\n";
		// cout << "\tNodes: ";
		// forall_nodes(print_node, N_t)
		// 	N_t.print_node(print_node);
		// cout << "\n\tEdges:\n";
		// forall_edges(print_edge, N_t){
		// 	N_t.print_edge(print_edge);
		// 	cout << '\n';
		// }
		// cout << '\n';

		// cout << "Copied original network:\n";
		// cout << "\tNodes: ";
		// forall_nodes(print_node, N_t_copied)
		// 	N_t_copied.print_node(print_node);
		// cout << "\n\tEdges:\n";
		// forall_edges(print_edge, N_t_copied){
		// 	N_t_copied.print_edge(print_edge);
		// 	cout << '\n';
		// }
		// cout << '\n';

		// Transform N_t_prime to N_t_prime_prime to perform min cut on nodes
		graph N_t_copied_copied; // Representing N_t_prime
		N_t_copied_copied = N_t; // Use N_t as N_t_prime_prime
		std::vector<edge> node_dup_edges; // For edges between node and its duplication

		// To avoid adding new node in N_t and looping each node in N_t at the same time
		std::vector<node> current_total_node_list;
		forall_nodes(print_node, N_t){
			current_total_node_list.push_back(print_node);
		}
		for(int i = 0; i < current_total_node_list.size(); ++i){
			// Duplicate each node except for t_prime and source
			if(current_total_node_list[i] == t_prime || current_total_node_list[i] == source_node) continue;

			// cout << "Splitting node ";
			// N_t.print_node(current_total_node_list[i]);
			// cout << " now...\n";

			node node_dup = N_t.new_node();
			std::string node_dup_name = Rev_N_t_node_mapping[current_total_node_list[i]] + "_dup";
			Rev_N_t_node_mapping[node_dup] = node_dup_name;

			// cout << Rev_N_t_node_mapping[current_total_node_list[i]] << " is now splitting into ";
			// cout << node_dup_name << " and ";
			// cout << Rev_N_t_node_mapping[current_total_node_list[i]] << '\n';


			// Link the input nodes of original node to duplicated node
			list<edge> input_edges = N_t.in_edges(current_total_node_list[i]);
			edge in_edge;
			forall(in_edge, input_edges){
				edge new_in_edge = N_t.new_edge(N_t.source(in_edge), node_dup);
				N_t.del_edge(in_edge);
			}
			// Link between duplicated node and original node
			edge node_edge = N_t.new_edge(node_dup, current_total_node_list[i]); // node_dup -> node_original
			node_dup_edges.push_back(node_edge);
			// N_t.print_node(node_dup);
			// cout << " -> ";
			// N_t.print_node(current_total_node_list[i]);
			// cout << '\n';
		}

		// Assign edge weights
		// 1 for edge between node and its duplication, 10000 for other edges
		edge_array<int> weight(N_t);
		forall_edges(print_edge, N_t){
			weight[print_edge] = 10000; // Can not use INT_MAX here, or it would be overflow
		}
		for(int i = 0; i < node_dup_edges.size(); ++i){
			weight[node_dup_edges[i]] = 1;
		}

		// cout << "Split network:\n";
		// cout << "\tNodes: ";
		// forall_nodes(print_node, N_t)
		// 	N_t.print_node(print_node);
		// cout << "\n\tEdges:\n";
		// forall_edges(print_edge, N_t){
		// 	N_t.print_edge(print_edge);
		// 	cout << ", weight: " << weight[print_edge] << '\n';
		// }
		// cout << '\n';

		// Min cut algorithm (using max flow to achieve min cut)
		// list<node> cut;
		// int cut_value = MIN_CUT(N_t, weight, cut);
		// Since MIN_CUT in LEDA is just a cut with minimum weights and it does not guarantee the cut separates s from t
		edge_array<int> flow(N_t);
		list<node> cut;
		int cut_value = MAX_FLOW(N_t, source_node, t_prime, weight, flow, cut);

		// cout << "****************************************\n";
		// cout << "*  Nodes to Name Mapping Table\n";
		// for(std::map<node, std::string>::iterator iter = Rev_N_t_node_mapping.begin();
		// 											iter != Rev_N_t_node_mapping.end(); ++iter){
		// 	cout << "*     " << iter->second << " -> ";
		// 	N_t.print_node(iter->first);
		// 	cout << '\n';
		// }
		// cout << "****************************************\n";


		// cout << "The minimum cut has value: " << cut_value << '\n';
		// cout << "cut:\n";

		// Build a map for node in the cut for later searching
		node v;
		std::map<node, bool> cut_node_mapping;
		forall(v, cut){
			// N_t.print_node(v);
			cut_node_mapping[v] = true;
		}
		// cout << '\n';

		// Cut nodes are the input signal of KLUT
		std::queue<node> KLUT_input_extract_queue;
		std::vector<std::string> KLUT_inputs;
		std::map<std::string, bool> KLUT_inputs_visited;

		// int print_out = 0;
		// if(Rev_N_t_node_mapping[t_prime] == "i_prime") print_out = 1;

		KLUT_input_extract_queue.push(t_prime);
		while(!KLUT_input_extract_queue.empty()){
			node input_check_node = KLUT_input_extract_queue.front();
			KLUT_input_extract_queue.pop();

			// Find all input nodes of input_check_node
			list<edge> input_edges = N_t.in_edges(input_check_node);
			edge in_edge;
			forall(in_edge, input_edges){
				// if(print_out) cout << "Check node " << Rev_N_t_node_mapping[N_t.source(in_edge)] << " as child of " << Rev_N_t_node_mapping[input_check_node] << '\n';
				if(cut_node_mapping.count(N_t.source(in_edge)) > 0){ // Input node is in cut set
					// if(print_out) cout << "Find a cut between " << Rev_N_t_node_mapping[N_t.source(in_edge)] << " and " << Rev_N_t_node_mapping[input_check_node] << '\n';
					if(KLUT_inputs_visited.count(Rev_N_t_node_mapping[input_check_node]) == 0){ // If the node is already in the list then skip
						KLUT_inputs.push_back(Rev_N_t_node_mapping[input_check_node]);
						KLUT_inputs_visited[Rev_N_t_node_mapping[input_check_node]] = true;
					}
				}
				else{ // Input node is NOT in cut set
					// Put the input node into the KLUT_input_extract_queue for further search
					KLUT_input_extract_queue.push(N_t.source(in_edge));
				}
			}	
		}
		
		if(cut_value <= k){ // node t can be packed with the nodes in X_t_bar
			FlowMap_NODE_label[NODE_front] = p;

			// cout << "Input signals for " << Rev_N_t_node_mapping[t_prime] << ": ";
			// for(int i = 0; i < KLUT_inputs.size(); ++i){
			// 	cout << KLUT_inputs[i] << ' ';
			// }
			// cout << '\n';
		}
		else{ // node t is self cut out
			FlowMap_NODE_label[NODE_front] = p+1;

			// Put the node t itself to the input signal set for no cut found
			KLUT_inputs.clear(); // Clear out the unfeasible inputs
			KLUT_inputs.push_back(NODE_front);
		}

		ASSERT(KLUT_inputs.size() <= k, "Input signal number are more than k!");
		// Map the vector of input signals of KLUT with the output node t
		t_KLUT_input_map[NODE_front] = KLUT_inputs;

		// cout << "l(" << NODE_front << "): " << FlowMap_NODE_label[NODE_front] << '\n';
		// cout << "-----------------------------------------------\n";
    }

    cout << "\rK-LUT......." << std::flush;
    // cout << "\rK-LUT.......\n";

    // for(std::map<std::string, std::vector<std::string> >::iterator it = t_KLUT_input_map.begin();
    // 																it != t_KLUT_input_map.end(); ++it){
    // 	cout << "Input signals for " << it->first << ": ";
    // 	for(int i = 0; i < it->second.size(); ++i){
    // 		cout << it->second[i] << ' ';
    // 	}
    // 	cout << "\n\n";
    // }

    // FlowMap Mapping Phase
    // Build a map for searching PI NODEs
    std::map<std::string, bool> PIs_mapping;
    for(int i = 0; i < PIs.size(); ++i){
    	PIs_mapping[PIs[i]] = true;
    }

    std::stack<std::string> mapping_phase_traversal_list;
    for(int i = 0; i < POs.size(); ++i){
    	mapping_phase_traversal_list.push(POs[i]);
    }

    // New graph for final K-LUT covered graph
    std::map<std::string, NODE> final_graph_mapping;
    std::map<std::string, bool> visited_input_nodes_mapping;

	outblif.open(output_file_name);

	outblif << model_name_raw << '\n';
	outblif << primary_input_raw << '\n';
	outblif << primary_output_raw << '\n';

    // cout << model_name_raw << '\n';
	// cout << primary_input_raw << '\n';
	// cout << primary_output_raw << '\n';

    // cout << "LUT output signals: ";
    LUT_count = 0; // Count LUT number
    while(!mapping_phase_traversal_list.empty()){
    	std::string node_v = mapping_phase_traversal_list.top();
    	mapping_phase_traversal_list.pop();

    	if(PIs_mapping.count(node_v) == 1) continue; // If this node is PI then skip
    	if(visited_input_nodes_mapping.count(node_v) == 1 && visited_input_nodes_mapping[node_v] == true) continue; // If this node is visited

    	// LUT_count++;
    	visited_input_nodes_mapping[node_v] = true;
    	// cout << node_v << ' ';

    	// New NODE
    	NODE NODE_v;
    	NODE_v.NODE_name = node_v;
    	NODE_v.NODE_level = 0;

    	// outblif << ".names ";
    	// // cout << ".names ";
    	if(t_KLUT_input_map[node_v].size() == 1 && t_KLUT_input_map[node_v][0] == node_v){
    		LUT_count++;
    		outblif << ".names ";
    		// cout << ".names ";

    		// If the cut contain only one node, add that node's input NODEs into the list
    		for(int i = 0; i < NODES_map[node_v].InNODEs.size(); ++i){
    			// if(PIs_mapping.count(NODES_map[node_v].InNODEs[i]) == 0){ // Input NODE is not in PIs
    				// if(visited_input_nodes_mapping.count(NODES_map[node_v].InNODEs[i]) == 0){ // Input NODE is not visited
    					mapping_phase_traversal_list.push(NODES_map[node_v].InNODEs[i]);
    					NODE_v.InNODEs.push_back(NODES_map[node_v].InNODEs[i]);
    				// }
    			// }
    			outblif << NODES_map[node_v].InNODEs[i] << ' ';
    			// cout << NODES_map[node_v].InNODEs[i] << ' ';
    		}
    		outblif << node_v << '\n';
    		// cout << node_v << '\n';
    		if(NODES_map[node_v].NODE_type == "1"){
				// cout << "AND";
			
				outblif << "11 1\n";
				// cout << "11 1\n";
	    	}
			else if(NODES_map[node_v].NODE_type == "2"){
				// cout << "OR";
				
				outblif << "1- 1\n";
				outblif << "-1 1\n";
				// cout << "1- 1\n";
				// cout << "-1 1\n";
			}
			else if(NODES_map[node_v].NODE_type == "3"){
				// cout << "INVERTER";
				outblif << "0 1\n";
				// cout << "0 1\n";
			}
			else if(NODES_map[node_v].NODE_type == "4"){
				// cout << "CONSTANT_0";
				continue;
			}
			else if(NODES_map[node_v].NODE_type == "5"){
				// cout << "CONSTANT_1";
				outblif << "1\n";
				// cout << "1\n";
			}
			else if(NODES_map[node_v].NODE_type == "6"){
				// cout << "BUFFER";
				outblif << "1 1\n";
				// cout << "1 1\n";
			}
			else if(NODES_map[node_v].NODE_type == "7"){
				// cout << "UNKNOWN\n";
				ASSERT(true, "UNKNOWN gate type!");
			}
    	}
    	else{
    		// // t_KLUT_input_map[node_v] has the input nodes for the cut
    		// for(int i = 0; i < t_KLUT_input_map[node_v].size(); ++i){
    		// 	// if(PIs_mapping.count(t_KLUT_input_map[node_v][i]) == 0){ // Input NODE is not in PIs
    		// 		// if(visited_input_nodes_mapping.count(t_KLUT_input_map[node_v][i]) == 0){ // Input NODE is not visited
    		// 			mapping_phase_traversal_list.push(t_KLUT_input_map[node_v][i]);
    		// 			NODE_v.InNODEs.push_back(t_KLUT_input_map[node_v][i]);
    		// 		// }
    		// 	// }
    		// 	outblif << t_KLUT_input_map[node_v][i] << ' ';
    		// 	// cout << t_KLUT_input_map[node_v][i] << ' ';
    		// }
    		// outblif << node_v << '\n';
    		// // cout << node_v << '\n';

    		std::map<std::string, bool> klut_visited;
    		std::map<std::string, int> node_output_every_round;
    		std::bitset<6> test_pattern_bit;
    		std::string test_pattern;
    		int output_signal;
    		std::vector<std::string> test_pattern_to_1;
    		for(int i = 0; i < std::pow(2,t_KLUT_input_map[node_v].size()); ++i){
    			klut_visited.clear(); // Clear the visited map every round
    			node_output_every_round.clear(); // Clear the node output map every round

    			// test_pattern generation
    			test_pattern_bit = i;
    			test_pattern = test_pattern_bit.to_string();
    			test_pattern = test_pattern.substr(test_pattern.size() - t_KLUT_input_map[node_v].size());
    			// cout << "Testing: " << test_pattern << '\n';

    			output_signal = klutdfs(node_v, klut_visited, node_output_every_round,
    									test_pattern, t_KLUT_input_map[node_v], NODES_map);
    			// cout << output_signal << '\n';
    			if(output_signal == 1) test_pattern_to_1.push_back(test_pattern);
    		}
    		// for(int i = 0; i < test_pattern_to_1.size(); ++i){
    		// 	outblif << test_pattern_to_1[i] << " 1\n";
    		// 	// cout << test_pattern_to_1[i] << " 1\n";
    		// }
    		if(test_pattern_to_1.size() > 0){
    			LUT_count++;
    			outblif << ".names ";
    			// cout << ".names ";

    			// t_KLUT_input_map[node_v] has the input nodes for the cut
    			for(int i = 0; i < t_KLUT_input_map[node_v].size(); ++i){
    				// if(PIs_mapping.count(t_KLUT_input_map[node_v][i]) == 0){ // Input NODE is not in PIs
    					// if(visited_input_nodes_mapping.count(t_KLUT_input_map[node_v][i]) == 0){ // Input NODE is not visited
    						mapping_phase_traversal_list.push(t_KLUT_input_map[node_v][i]);
    						NODE_v.InNODEs.push_back(t_KLUT_input_map[node_v][i]);
    					// }
    				// }
    				outblif << t_KLUT_input_map[node_v][i] << ' ';
    				// cout << t_KLUT_input_map[node_v][i] << ' ';
    			}
    			outblif << node_v << '\n';
    			// cout << node_v << '\n';
    			for(int i = 0; i < test_pattern_to_1.size(); ++i){
    				outblif << test_pattern_to_1[i] << " 1\n";
    				// cout << test_pattern_to_1[i] << " 1\n";
    			}
    		}
    		else continue;
    		ASSERT(test_pattern_to_1.size() > 0, "Truth table not found!");
    	}

    	final_graph_mapping[node_v] = NODE_v;
    }
    outblif << ".end";
	// cout << ".end\n";
	outblif.close();

    // DFS from PO to get the level of final graph
    std::map<std::string, bool> visited;
    std::map<std::string, int> level_map;
    int circuit_level = 0;
    for(int i = 0; i < POs.size(); ++i){
    	int level = dfs(POs[i], visited, level_map, final_graph_mapping);
    	if(level > circuit_level) circuit_level = level;
    }
    
 //    for(std::map<std::string, NODE>::iterator it = final_graph_mapping.begin(); it != final_graph_mapping.end(); ++it){
	// 	cout << "Index: " << it->first << " -> NODE: " << it->second.NODE_name << '\n';
	// }

    cout << "\rThe circuit level is " << circuit_level << ".\n";
    cout << "The number of LUTs is " << LUT_count << ".\n";

	return 0;
}
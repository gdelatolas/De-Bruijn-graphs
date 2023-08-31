#ifndef DEBRUIJNGRAPH_H
#define DEBRUIJNGRAPH_H

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>



enum class Amino : int32_t
{
    A = 1,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,
};


class De_Bruijn_Graph {
    private:

        // Returns 1 if the argument element exist in vector v.
        template <typename T>
        bool element_exists_in_vector(const std::vector<T>& v, const T& element, size_t start = 0);
        
        // Returns the index of an element at a vector.
        template <typename TT>
        int index_of_element_in_vector(const std::vector<TT>& v, int element);
        


        // Adds an edge from vertex s to vertex d.
        // We will have one edge per kmer.
        void add_edge(int s, int d);     

        // Simple DFS implementation.
        void dfs(int curr, std::vector<bool>& visited, std::vector<int>& path);
        
        // Returns 1 if the graph is fully-connected, otherwise 0.
        bool strongly_connected_graph();
        


        // Function that counts the in_degree in edges at one vertex v.
        int in_degree(int v);
        // Function that counts the out_degree in edges at one vertex v.
        int out_degree(int v);

    public:
        De_Bruijn_Graph(int lines, int line_size, int kmer, std::vector<std::vector<Amino>> aminos);
        ~De_Bruijn_Graph();
        
        // Creates the Graph.
        void create_the_graph();
        // Print the graph using: m_adj.
        void print_graph();

        // Just calls the function: print_eulerian_path_cycle.
        void create_euler_path_cycle();

        // Creates the eulerian path cycle, fill the m_final_path with the correct nodes.
        void print_eulerian_path_cycle(int start_node);
        
        // Prints the eulerian path cycle using m_final_path.
        void print_euler_path();
        
        // Checks if the graph has an eulerian path or not.
        int find_euler(int start_node);
        
        //==============================================================================================
        int m_lines, m_line_size, m_kmer, m_no_vertices;
        
        std::vector<int> m_k_1_mers_int;    // vector of integers.
                                            // each line has a k-1 mer at int form.
                                            // We have created this vector in order
                                            // to have acces to k_1_mers_int 
                                            // from every function of the class.
        
                                
        std::vector<std::string> m_k_1_mers_string;     // string vector it is used to print the graph correctly
                                                        // k-1 mers at string form.
                                                        // We have created this vector in order
                                                        // to have acces to k_1_mers_string 
                                                        // from every function of the class.
                                        
       
        std::vector<int> m_final_path;        
        

        std::vector<int> m_vec1; // We store the nodes at integer form

        std::vector<std::string> m_node;// Vector of strings. 
                                        // We will store the nodes at string form.
        
        std::vector<std::vector<int>> m_adj; // Array of std::vectors. 
                                 // This is our graph.

        std::vector<int> path; // It is used for DFS.
};

#endif // DEBRUIJNGRAPH_H

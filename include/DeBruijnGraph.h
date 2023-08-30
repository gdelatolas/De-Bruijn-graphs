#ifndef DEBRUIJNGRAPH_H
#define DEBRUIJNGRAPH_H

#include <vector>
#include <string>
#include <unordered_map>

enum class Amino : int32_t
{
    A = 1,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,
};

class De_Bruijn_Graph {
    private:
        template <typename T>
        bool elementExistsInVector(const std::vector<T>& v, const T& element, size_t start = 0);
        template <typename TT>
        int indexOfElementInVector(const std::vector<TT>& v, int element) ;
        
        void addEdge(int s, int d);     // We will have one edge per kmer.
        void DFS(int curr, std::vector<bool>& visited, std::vector<int>& path);
        
        bool Strongly_Connected_Graph();
        
        int in_degree(int v);
        int out_degree(int v);

    public:
        De_Bruijn_Graph();
        ~De_Bruijn_Graph();

        // Public member functions
        void create_the_graph();
        void printGraph();
        void printEulerianPathCycle(int start_node);
        int find_Euler(int start_node);
        void create_Euler_Path_Cycle();
        void print_Euler_path();

        // ... (other public member functions)
        int m_lines, m_line_size, m_kmer, m_no_vertices;
        int *m_k_1_mers_int_11; // array of integers.
                                // each line has a k-1 mer at int form.
        std::string *m_k_1_mers_string_11;  // string array it is used to print the graph correctly
                                            // k-1 mers at string form.
                                        
        Amino **m_k_1_mers11;   // 2-d Array. k-1 mers at Amino form

        std::vector<Amino> m_aminos;
        std::vector<int> m_final_path;        
        

        std::vector<int> m_vec1; // the two mers values.

        std::vector<std::string> m_node; // Vector of strings. 
                                        // We will store the nodes at string form.
        
        // Add edge, including self-loop
        std::vector<int> *m_adj; // Change to pointer
        std::vector<int> path; // we will store the eulerian path.s
};

#endif // DEBRUIJNGRAPH_H
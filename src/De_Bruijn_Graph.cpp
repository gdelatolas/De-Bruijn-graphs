#include "../include/DeBruijnGraph.h"
#include <iostream>

std::unordered_map<char, Amino> aminoMap = {
        {'A', Amino::A},
        {'B', Amino::B},
        {'C', Amino::C},
        {'D', Amino::D},
        {'E', Amino::E},
        {'F', Amino::F},
        {'G', Amino::G},
        {'H', Amino::H},
        {'I', Amino::I},
        {'J', Amino::J},
        {'K', Amino::K},
        {'L', Amino::L},
        {'M', Amino::M},
        {'N', Amino::N},
        {'O', Amino::O},
        {'P', Amino::P},
        {'Q', Amino::Q},
        {'R', Amino::R},
        {'S', Amino::S},
        {'T', Amino::T},
};

// Constructor
De_Bruijn_Graph::De_Bruijn_Graph() {

    std::cout << "Input the number of lines: ";
    std::cin >> m_lines;

    std::cout << "Input the size of each line: ";
    std::cin >> m_line_size;

    //The vector aminos stores our input lines.
    std::vector<std::vector<Amino>> aminos(m_lines, std::vector<Amino>());

    for (int i = 0; i < m_lines; i++){
        
        std::cout << "Input a sequence of "<<m_line_size<<" amino acids as letters (A, B, ..., T): ";
        std::string userInput;
        std::cin >> userInput;
        for (char c : userInput) {
            if (aminoMap.find(c) != aminoMap.end()) 
                aminos[i].push_back(aminoMap[c]);
            
            else 
                std::cout << "Invalid amino acid: " << c << "\n";
        } 
    }
    std::cout << "size = " << aminos.size() << "\n";
    //PRINTING THE aminos VECTOR.  
    std::cout << "\nAminos\n";
    for (int i = 0; i < m_lines; i++) {
        for (Amino amino : aminos[i]) {
            char aminoChar = static_cast<char>(static_cast<int>(amino) + 'A' - 1);
            std::cout << aminoChar;
        }
        std::cout << "\n";
    } 
        
    std::cout << "Input the k-mer you want: ";
    std::cin >> m_kmer;
    
    // X dimension of k_mer array must have the size (m_line_size-(m_kmer-1))*m_lines.
    // (m_line_size-(m_kmer-1)) is the number of k_mers we can extract from an input line,
    // For example, first line = ABCBBD (m_line_size = 6) and k_mer = 3 the 3-mers are: 
    // ABC, BCB , CBB, BDD 
    // So, we have four 3-mers, 6-(3-1) = 4.
    // We extract the same amount of k-mers from each input line, due to that
    // we multiply with m_line (number of the input lines).
    Amino k_mer[(m_line_size-(m_kmer-1))*m_lines][m_kmer];

    int counter = 0 ; // counter is used to identify the input line.
    int index = 0;    // index indicates the line of the k_mer array.
                        // index = numbers of total k_mers = (m_line_size - (m_kmer -1))* m_lines

    while (counter < m_lines){
        
        for (int i = 0; i < m_line_size - (m_kmer -1); i++){

            for(int j = 0; j < m_kmer; j++){
                k_mer[index][j] = aminos[counter][i+j];
            }
            index++;
        }
        counter++;
    }

    //FOR PRINTING K_MERS.
    std::cout << "\nKmers\n";
    for (int i = 0; i < (m_line_size-(m_kmer-1))*m_lines; i++){

        for(int j = 0; j < m_kmer; j++){

            char h = static_cast<char>(static_cast<int>(k_mer[i][j]) + 'A' - 1);
            std::cout << h;
        }
        std::cout << "\n";
    } 

    // We create K-1 mers.
    // X-dimension must have the size of (m_line_size-(m_kmer-1))*m_lines*2, 
    // because for each k_mer we create 2 k_1_mers.
    counter = 0;
    Amino k_1_mer [(m_line_size-(m_kmer-1))*m_lines*2][m_kmer-1];
    for (int i = 0; i < (m_line_size-(m_kmer-1))*m_lines*2; i+=2){
        
        for (int j = 0, c = 0; j <= m_kmer-2; j++, c++){
            k_1_mer[i][j] = k_mer[counter][j];
            k_1_mer[i+1][j] = k_mer[counter][j+1];
        }
        counter++;
    }
    
    // =====================
    // We create those two arrays k_1_mers_int and k_1_mers_string in order 
    // to have acces to those values by m_k_1_mers11 and m_k_1_mers_string_11
    // by the other functions of the class.
    int k_1_mers_int[(m_line_size-(m_kmer-1))*m_lines*2];
    std::string k_1_mers_string[(m_line_size-(m_kmer-1))*m_lines*2];
    
    std::cout<< "\nk-1 mers\n";
    for(int i=0; i < (m_line_size-(m_kmer-1))*m_lines*2; i++){
        std::string s="";
        std::string s2="";
        int helper;
        char c;
        for(int j =0 ; j < m_kmer-1; j++){
            //std::cout << static_cast<int>(k_1_mers[i][j]) << " ";
            helper = static_cast<int>(k_1_mer[i][j]);  
            c = static_cast<char>(static_cast<int>(k_1_mer[i][j]) + 'A' - 1);
            
            s = s +  std::to_string(helper);
            s2 += c;
        }
        k_1_mers_int[i] = std::stoi(s);
        k_1_mers_string[i] = s2;

    }
    // Dynamically allocate memory for two_mers11
    m_k_1_mers11 = new Amino*[(m_line_size-(m_kmer-1))*m_lines*2];
    for (int i = 0; i < ((m_line_size-(m_kmer-1))*m_lines*2); i++) {
        m_k_1_mers11[i] = new Amino[m_kmer-1];
        for (int j = 0; j<m_kmer-1; j++){
            m_k_1_mers11[i][j] = k_1_mer[i][j];
        }
    }
    m_k_1_mers_int_11 = new int [(m_line_size-(m_kmer-1))*m_lines*2];
    for (int i = 0; i < (m_line_size-(m_kmer-1))*m_lines*2; i++){
        m_k_1_mers_int_11[i] = k_1_mers_int[i];
    }

    
    m_k_1_mers_string_11 = new std::string [(m_line_size-(m_kmer-1))*m_lines*2];
    for (int i = 0; i < (m_line_size-(m_kmer-1))*m_lines*2; i++){
        m_k_1_mers_string_11[i] = k_1_mers_string[i];
    }
    
}

// Destructor
De_Bruijn_Graph::~De_Bruijn_Graph() {
    // Deallocate m_k_1_mers11
    for (int i = 0; i < (m_line_size-(m_kmer-1))*m_lines*2; i++) {
        delete[] m_k_1_mers11[i];
    
    }
    delete[] m_k_1_mers11;

    // Deallocate m_k_1_mers_int_11
    delete[] m_k_1_mers_int_11;

    // Deallocate m_k_1_mers_string_11
    delete[] m_k_1_mers_string_11;

    // Deallocate m_adj (the dynamically allocated array of vectors)
    delete[] m_adj;
}


//=========================================================================================================//
//          Function that checks if an element already exists in a vector.
template <typename T>
bool De_Bruijn_Graph::elementExistsInVector(const std::vector<T>& v, const T& element, size_t start) {
    std::vector<T> temp = v;
    for (size_t i = start; i < temp.size(); i++) {
        if (temp[i] == element)
            return true;
    }
    return false;
}


//=========================================================================================================//
//          Function that returns the index of an element at a vector.
template <typename TT>
int De_Bruijn_Graph::indexOfElementInVector(const std::vector<TT>& v, int element) {
    std::vector<TT> temp = v;
    for (int i = 0; i < temp.size(); i++){
        if (temp[i] == element)
            return i;
    }
    return -1;
}


//=========================================================================================================//
//          Function that adds an edge
void De_Bruijn_Graph::addEdge(int s, int d) {
    m_adj[s].push_back(d);
}


void De_Bruijn_Graph::create_the_graph(){
    int counter = 0;
    for (int i = 0; i < (m_line_size-(m_kmer-1))*m_lines*2; i++){
        if(! elementExistsInVector(m_vec1, m_k_1_mers_int_11[i])){
            m_vec1.push_back(m_k_1_mers_int_11[i]);  
             m_node.push_back(m_k_1_mers_string_11[i]); 
            counter++;
        }
        
    }
    // We will have one node per distinct k-1 mer.
    m_no_vertices = counter;
    m_adj = new std::vector<int>[m_no_vertices];

    // We use this vector in order to store the edges as strings.
    std::vector<std::string> vec_edges;

    std::cout << "\nVECTOR\n";
    for (int i = 0; i < m_vec1.size(); i ++){
        std::cout << m_vec1[i] << "\n";
    }
    int i_s,i_d;
    for (int i = 0; i < (m_line_size-(m_kmer-1))*m_lines*2; i+=2){
        i_s = indexOfElementInVector(m_vec1, m_k_1_mers_int_11[i]);
        i_d = indexOfElementInVector(m_vec1, m_k_1_mers_int_11[i+1]);
        std::string i_s_String = std::to_string(i_s);
        std::string i_d_String = std::to_string(i_d);


        // We check if the edge already exist.
        if(!elementExistsInVector(vec_edges, i_s_String + i_d_String)){
            vec_edges.push_back(i_s_String + i_d_String);
            addEdge(i_s, i_d);
        }

    }

}
// Print the graph
// Keep in mind that the walk that allows us to reconstruct the genome 
// is the walk that crosses each edge exactly once, wich means that, 
// uses each k-mer exactly once.
// This is the eulerian walk.
void De_Bruijn_Graph::printGraph() {
    for (int d = 0; d < m_no_vertices; ++d) {
        std::cout << "\n Vertex " << m_node[d] << ":";
        for (auto x : m_adj[d]){
            std::cout << "-> " << m_node[x];
        }
        std::cout<<"\n";
    }
}
void De_Bruijn_Graph::printEulerianPathCycle(int start_node)
{
    std::vector<int> circuit;

    while(m_adj[start_node].size())
    {
        int next_node = m_adj[start_node].back();
        m_adj[start_node].pop_back();
        printEulerianPathCycle(next_node);
    }

    circuit.push_back(start_node);
    
    for(auto node: circuit)
        m_final_path.push_back(node);
}
void De_Bruijn_Graph::print_Euler_path(){
    for(int i = m_final_path.size() - 1; i >= 0; i--){
        int node = m_final_path[i];
        if(i != 0)
            std::cout << m_node[node] << " -> ";
        else
            std::cout << m_node[node];
    }
}

int De_Bruijn_Graph::find_Euler(int start_node)
{
    if(!Strongly_Connected_Graph())	{//multi-component edged graph
        std::cout << "The graph is not fully-connected.";
        return 0;		//All non-zero degree vertices should be connected

    } 
    // Counting in-degrees and out-degrees
    int start_nodes = 0, end_nodes = 0;
    for(int i = 0; i < m_no_vertices; ++i)
    {
        int in = in_degree(i);
        int out = out_degree(i);
        if(out - in > 1 or in - out > 1)
            return 0;
        else if(out - in == 1)
        {
            start_nodes++;
            start_node = i;
            std::cout << "m_node[i] = " << m_node[i] << "\n";
        }
        else if(in - out == 1)
            end_nodes++;
    }
    std::cout << start_nodes << "  " << end_nodes << "\n";
    // Only start and end nodes can have in-degree and out-degree difference
    if(start_nodes == end_nodes and (start_nodes == 0 or start_nodes == 1))
        return (start_nodes == 0) ? 2 : 1;

    return 0;
}

void De_Bruijn_Graph::create_Euler_Path_Cycle()
{
    int start_node = 0;
    int ret = find_Euler(start_node);
    if(ret == 0)
        std::cout << "Graph is NOT an Euler graph\n";
    else if(ret == 1)
    {
        std::cout << "Graph is Semi-Eulerian\n";
        printEulerianPathCycle(start_node);
    }
    else
    {
        std::cout << "Graph is Eulerian (Euler circuit)\n";
        printEulerianPathCycle(start_node);
    }
    std::cout << "\n";
}
void De_Bruijn_Graph::DFS(int curr, std::vector<bool>& visited, std::vector<int>& path)
{
    visited[curr] = true;
    path.push_back(curr);
    for(auto it: m_adj[curr])
    {
        if(!visited[it])
            DFS(it, visited, path);
    }
}
bool De_Bruijn_Graph::Strongly_Connected_Graph()
{
    std::vector<bool> visited(m_no_vertices, false);
    int node = -1;	//Node to start DFS
    for(int i = 0; i < m_no_vertices; ++i)
        if(m_adj[i].size() > 0)
        {
            node = i;	//Found a node to start DFS
            break;
        }
    if(node == -1)	//No start node was found-->No edges are present in graph
        return true; //It's always Eulerian

    DFS(node, visited, path);
    //Check if all the non-zero vertices are visited
    for(int i = 0; i < m_no_vertices; ++i)
        if(visited[i] == false and m_adj[i].size() > 0)
            return false;	//We have edges in multi-component
    return true;
}
int De_Bruijn_Graph::in_degree(int v)
{
    int degree = 0;
    for(int i = 0; i < m_no_vertices; i++)
    {
        for(auto it: m_adj[i])
        {
            if(it == v)
                degree++;
        }
    }
    return degree;
}

int De_Bruijn_Graph::out_degree(int v)
{
    return m_adj[v].size();
}
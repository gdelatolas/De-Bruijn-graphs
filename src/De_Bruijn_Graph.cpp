#include "../include/DeBruijnGraph.h"


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
De_Bruijn_Graph::De_Bruijn_Graph(int lines, int kmer, std::vector<std::vector <Amino>> k_mer_vec){

    m_lines = lines;

   
    m_kmer = kmer;
    
    // X dimension of k_mer array must have the size (m_line_size-(m_kmer-1))*m_lines.
    // (m_line_size-(m_kmer-1)) is the number of k_mers we can extract from an input line,
    // For example, first line = ABCBBD (m_line_size = 6) and k_mer = 3 the 3-mers are: 
    // ABC, BCB , CBB, BDD 
    // So, we have four 3-mers, 6-(3-1) = 4.
    // We extract the same amount of k-mers from each input line, due to that
    // we multiply with m_line (number of the input lines).


    m_k_mer_size = k_mer_vec.size();

    //FOR PRINTING K_MERS.
    std::cout << "\nKmers\n";
    for (int i = 0; i < m_k_mer_size; i++){

        for(int j = 0; j < m_kmer; j++){

            char h = static_cast<char>(static_cast<int>(k_mer_vec[i][j]) + 'A' - 1);
            std::cout << h;
        }
        std::cout << "\n";
    } 

    // We create K-1 mers.
    // X-dimension must have the size of (m_line_size-(m_kmer-1))*m_lines*2, 
    // because for each k_mer we create 2 k_1_mers.
    int counter = 0;
    Amino k_1_mer [m_k_mer_size*2][m_kmer-1];   //2-D Amino array.
    
    for (int i = 0; i < m_k_mer_size*2; i+=2){
        
        for (int j = 0, c = 0; j <= m_kmer-2; j++, c++){
            k_1_mer[i][j] = k_mer_vec[counter][j];
            k_1_mer[i+1][j] = k_mer_vec[counter][j+1];
        }
        counter++;
    }
    
    std::cout<< "\nk-1 mers\n";
    for(int i=0; i < m_k_mer_size*2; i++){
        std::string s="";
        std::string s2="";
        int helper;
        char c;
        for(int j =0 ; j < m_kmer-1; j++){
            // For example if k_1_mer[i][j] =  A then helper =  1
            helper = static_cast<int>(k_1_mer[i][j]);
            
            // For example if k_1_mer[i][j] =  A then c =  A  
            c = static_cast<char>(static_cast<int>(k_1_mer[i][j]) + 'A' - 1);
            
            s = s +  std::to_string(helper);
            s2 += c;
        }
        m_k_1_mers_int.push_back(std::stoi(s));
        m_k_1_mers_string.push_back(s2);
    }
    create_the_graph();
}


//=========================================================================================================//
//          Function that creates the graph.
void De_Bruijn_Graph::create_the_graph(){
    int counter = 0;
    for (int i = 0; i < m_k_mer_size*2; i++){
        if(! element_exists_in_vector(m_vec1, m_k_1_mers_int[i])){
            m_vec1.push_back(m_k_1_mers_int[i]);  
            m_node.push_back(m_k_1_mers_string[i]); 
            counter++;
        }
        
    }
    // We will have one node per distinct k-1 mer.
    m_no_vertices = counter;
    
    m_adj.resize(m_no_vertices);
    
    // We use this vector in order to store the edges as strings.
    std::vector<std::string> vec_edges;
    
    std::cout << "\nVECTOR\n";
    for (int i = 0; i < m_vec1.size(); i ++){
        std::cout << m_vec1[i] << "\n";
    }
    int i_s, i_d;
    for (int i = 0; i < m_k_mer_size*2; i+=2){
        i_s = index_of_element_in_vector(m_vec1, m_k_1_mers_int[i]);
        i_d = index_of_element_in_vector(m_vec1, m_k_1_mers_int[i+1]);

        // index_of_element_in_vector always return a value because every element of
        // the arrau m_k_1_int exist exactly once at the m_vec1 vector.
        std::string i_s_String = std::to_string(i_s);
        std::string i_d_String = std::to_string(i_d);


        // We check if the edge already exist.
        if(!element_exists_in_vector(vec_edges, i_s_String + i_d_String)){
            vec_edges.push_back(i_s_String + i_d_String);
            add_edge(i_s, i_d);
        }

    }

    // Compute in_degree and out_degree for each node.
    for(int i = 0; i < m_no_vertices; ++i)
    {
        in__degree.push_back( in_degree(i));
        out__degree.push_back(out_degree(i));
    }

}




//=========================================================================================================//
//          Destructor.
De_Bruijn_Graph::~De_Bruijn_Graph() {}



//auto it_9    = std::find(seq.begin(), seq.end(), 9);

//auto index_9 = (index_t)std::distance(seq.begin(), it_9);


//=========================================================================================================//
//          Function that checks if an element already exists in a vector.
template <typename T>
bool De_Bruijn_Graph::element_exists_in_vector(const std::vector<T>& v, const T& element) {
    std::vector<T> temp = v;
    auto finder = std::find(v.begin(), v.end(), element);
    if (finder == v.end())
        return false;
    return true;
}




//=========================================================================================================//
//          Function that returns the index of an element at a vector.
template <typename TT>
int De_Bruijn_Graph::index_of_element_in_vector(const std::vector<TT>& v, int element) {
    auto finder = std::find(v.begin(), v.end(), element);
    
    if (finder == v.end())
        return -1; // Element not found, return -1
    else
        return std::distance(v.begin(), finder); // Return the index of the found element
}






//=========================================================================================================//
//          Function that adds an edge
void De_Bruijn_Graph::add_edge(int s, int d) {
    m_adj[s].push_back(d);
}





// Print the graph
// Keep in mind that the walk that allows us to reconstruct the genome 
// is the walk that crosses each edge exactly once, wich means that, 
// uses each k-mer exactly once.
// This is the eulerian walk.
void De_Bruijn_Graph::print_graph() {
    for (int d = 0; d < m_no_vertices; ++d) {
        std::cout << "\n Vertex " << m_node[d] << ":";
        for (auto x : m_adj[d]){
            std::cout << "-> " << m_node[x];
        }
        std::cout<<"\n";
    }
}







void De_Bruijn_Graph::print_eulerian_path_cycle(int start_node)
{
    std::vector<int> circuit;

    while(m_adj[start_node].size())
    {
        int next_node = m_adj[start_node].back();
        m_adj[start_node].pop_back();
        print_eulerian_path_cycle(next_node);
    }

    circuit.push_back(start_node);
    
    for(auto node: circuit)
        m_final_path.push_back(node);
}








void De_Bruijn_Graph::print_euler_path(){
    for(int i = m_final_path.size() - 1; i >= 0; i--){
        int node = m_final_path[i];
        if(i != 0)
            std::cout << m_node[node] << " -> ";
        else
            std::cout << m_node[node];
    }
}








int De_Bruijn_Graph::find_euler(int start_node)
{
    if(!strongly_connected_graph())	{   //multi-component edged graph
        std::cout << "The graph is not fully-connected.";
        return 0;		//All non-zero degree vertices should be connected

    } 
    // Counting in-degrees and out-degrees
    int start_nodes = 0, end_nodes = 0;
    for(int i = 0; i < m_no_vertices; ++i)
    {
        int in = in__degree[i];
        int out = out__degree[i];
        if(out - in > 1 or in - out > 1)
            return 0;
        else if(out - in == 1)
        {
            start_nodes++;
            start_node = i;
        }
        else if(in - out == 1)
            end_nodes++;
    }
    
    // Only start and end nodes can have in-degree and out-degree difference
    if(start_nodes == end_nodes and (start_nodes == 0 or start_nodes == 1))
        return (start_nodes == 0) ? 2 : 1;

    return 0;
}










void De_Bruijn_Graph::create_euler_path_cycle()
{
    int start_node = 0;
    int ret = find_euler(start_node);
    if(ret == 0)
        std::cout << "Graph is NOT an Euler graph\n";
    else if(ret == 1)
    {
        std::cout << "Graph is Semi-Eulerian\n";
        print_eulerian_path_cycle(start_node);
    }
    else
    {
        std::cout << "Graph is Eulerian (Euler circuit)\n";
        print_eulerian_path_cycle(start_node);
    }
    std::cout << "\n";
}






void De_Bruijn_Graph::dfs(int curr, std::vector<bool>& visited, std::vector<int>& path)
{
    visited[curr] = true;
    path.push_back(curr);
    for(auto it: m_adj[curr])
    {
        if(!visited[it])
            dfs(it, visited, path);
    }
}






bool De_Bruijn_Graph::strongly_connected_graph()
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

    dfs(node, visited, path);
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
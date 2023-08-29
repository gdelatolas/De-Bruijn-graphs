#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

enum class Amino : int32_t
{
    A = 1,
    B,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    J,
    K,
    L,
    M,
    N,
    O,
    P,
    Q,
    R,
    S,
    T,
};
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

/* 
Input a sequence of 5 amino acids as letters (A, B, ..., T): ABCBB
Input a sequence of 5 amino acids as letters (A, B, ..., T): KLMNO
Input a sequence of 5 amino acids as letters (A, B, ..., T): FGEHJ
 */
class De_Bruijn_Graph {

    public:
        int m_lines, m_line_size, m_kmer, m_no_vertices;
        int *m_k_1_mers_int_11; // k-1 mers at int form
        std::string *m_k_1_mers_string_11; // k-1 mers at string form.
                                           // This pointer-array it is used to print the graph correctly.
        Amino **m_k_1_mers11;   // k-1 mers at Amino form

        std::vector<Amino> m_aminos;

        std::vector<int> m_vec1; // the two mers values.
        std::vector<int> m_vec2; // two mers indicies at the two merss array.
        std::vector<int> m_vec3; // values 0 1 2 .. no_vertices-1.
        std::vector<std::string> m_node; // We will store the nodes at string form.
        
        // Add edge, including self-loop
        std::vector<int> *m_adj; // Change to pointer
        

        De_Bruijn_Graph(){

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

        template <typename T>
        bool elementExistsInVector(const std::vector<T>& v, const T& element, size_t start = 0) {
            std::vector<T> temp = v;
            for (size_t i = start; i < temp.size(); i++) {
                if (temp[i] == element)
                    return true;
            }
            return false;
        }
        template <typename TT>
        int indexOfElementInVector(const std::vector<TT>& v, int element) {
            std::vector<TT> temp = v;
            for (int i = 0; i < temp.size(); i++){
                if (temp[i] == element)
                    return i;
            }
            return -1;
        }
        // Add edge, including self-loop
        // We will have one edge per kmer.
        void addEdge(int s, int d) {
            m_adj[s].push_back(d);
        }
        

        void create_the_graph(){
            int counter = 0;
            for (int i = 0; i < (m_line_size-(m_kmer-1))*m_lines*2; i++){
                if(! elementExistsInVector(m_vec1, m_k_1_mers_int_11[i])){
                    m_vec1.push_back(m_k_1_mers_int_11[i]);  
                    m_vec2.push_back(i);             
                    m_vec3.push_back(counter);
                    m_node.push_back(m_k_1_mers_string_11[i]); 
                    counter++;
                }
                
            }
            // We will have one node per distinct k-1 mer.
            m_no_vertices = counter;
            m_adj = new std::vector<int>[m_no_vertices];
            
            // We use this vector in order to store the edges as strings.
            std::vector<std::string> vec_edges;

            std::cout << "\nVECTORS\n";
            for (int i = 0; i < m_vec1.size(); i ++){
                std::cout << m_vec1[i] << " " << m_vec2[i] << "  "<< m_vec3[i] << "\n";
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
        void printGraph() {
            for (int d = 0; d < m_no_vertices; ++d) {
                std::cout << "\n Vertex " << m_node[d] << ":";
                for (auto x : m_adj[d]){
                    std::cout << "-> " << m_node[x];
                }
                std::cout<<"\n";
            }
        }
        // Destructor
        ~De_Bruijn_Graph() {
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
        
};

int main(){

    De_Bruijn_Graph g1;
    g1.create_the_graph();
    g1.printGraph();

}

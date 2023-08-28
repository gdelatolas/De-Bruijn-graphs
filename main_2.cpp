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

class De_Bruijn_Graph {
    public:
        int m_size, m_no_vertices, m_size_reads, m_kmer;
        int *m_k_1_mers_int_11;
        Amino **m_k_1_mers11;
        
        std::vector<Amino> m_aminos;

        std::vector<int> m_vec1; // the two mers values.
        std::vector<int> m_vec2; // two mers indicies at the two merss array.
        std::vector<int> m_vec3; // values 0 1 2 .. no_vertices-1.
        // Add edge, including self-loop
        std::vector<int> *m_adj; // Change to pointer
        

        De_Bruijn_Graph (int n, int size_reads, int k){  //n = size of sequence 
            m_size = n;
            m_size_reads = size_reads;
            m_kmer = k;
            std::cout << "Input a sequence of "<<n<<" amino acids as letters (A, B, ..., T): ";
            std::string userInput;
            std::cin >> userInput;

        
            for (char c : userInput) {
                if (aminoMap.find(c) != aminoMap.end()) {
                    m_aminos.push_back(aminoMap[c]);
                } 
                else {
                    std::cout << "Invalid amino acid: " << c << "\n";
                }
            }
            Amino k_mers[m_size-(m_kmer-1)][m_kmer];
            for(int i = 0; i < m_size-(m_kmer-1); i++){
                for(int j = 0; j < m_kmer; j++){
                    k_mers[i][j] = m_aminos[i+j];
                } 
            }

            Amino k_1_mers[(m_size-(m_kmer-1))*2][m_kmer-1]; // k-1 mers
            int counter = 0;
            for (int i = 0; i < (m_size - (m_kmer-1))* 2; i+=2){
                for(int j = 0 ,c = 0; j  <= m_kmer; j++,c++){
                    if (j <= m_kmer-2){
                        k_1_mers[i][j] = k_mers[counter][j];
                    }
                    if (j+1 <= m_kmer -1){
                        k_1_mers[i+1][j] = k_mers[counter][j+1];
                    }
                }
                counter++;
            }
                      
            // =====================
            // FOR TESTING
            int k_1_mers_int[(m_size-(m_kmer-1))*2];

            std::cout<< "\nk-1 mers\n";
            for(int i=0; i < (m_size-(m_kmer-1))*2; i++){
                std::string s="";
                int helper;
                for(int j =0 ; j < m_kmer-1; j++){
                    //std::cout << static_cast<int>(k_1_mers[i][j]) << " ";
                    helper = static_cast<int>(k_1_mers[i][j]);  
                    s = s +  std::to_string(helper);
                }
                //std::cout << s << "\n";
                k_1_mers_int[i] = std::stoi(s);
                std::cout << "\n";
            }
            // =====================
            // Dynamically allocate memory for two_mers11
            m_k_1_mers11 = new Amino*[(m_size - (m_kmer-1))*2];
            for (int i = 0; i < ((m_size - (m_kmer-1))*2); i++) {
                m_k_1_mers11[i] = new Amino[m_kmer-1];
                for (int j = 0; j<m_kmer-1; j++){
                    m_k_1_mers11[i][j] = k_1_mers[i][j];
                }
            }
            m_k_1_mers_int_11 = new int [(m_size-(m_kmer-1))*2];
            for (int i = 0; i < (m_size-(m_kmer-1))*2; i++){
                m_k_1_mers_int_11[i] = k_1_mers_int[i];
            }
        }
        bool elementExistsInVector(const std::vector<int>& v, int element, int start = 0) {
            std::vector<int> temp = v;
            // start is the element we search from this an beyond.
            for (int i = start; i < temp.size(); i++){
                if (temp[i] == element)
                    return true;
            }
            return false;
        }

        
        // Add edge, including self-loop
        void addEdge(int s, int d) {
            m_adj[s].push_back(d);
        }

        
        void create_the_graph(){
            
            int counter = 0;

            //int kmerss[(size-(kmer-1))*(kmer-1)];
            
            for(int i = 0; i < (m_size-(m_kmer-1))*2;  i+=2){

                //kmerss[i] = static_cast<int>(m_k_1_mers11[i][0]) *10 + static_cast<int>(m_k_1_mers11[i][1]);
                //kmerss[i+1] = static_cast<int>(m_k_1_mers11[i+1][0]) *10 + static_cast<int>(m_k_1_mers11[i+1][1]);
            
                if(! elementExistsInVector(m_vec1, m_k_1_mers_int_11[i])){
                    m_vec1.push_back(m_k_1_mers_int_11[i]);  
                    m_vec2.push_back(i);             
                    m_vec3.push_back(counter);
                    counter++;
                }
                if (i + (m_kmer-1) >= (2)*(m_size-(m_kmer-1))){
                    if(! elementExistsInVector(m_vec1, m_k_1_mers_int_11[i+1])){
                        m_vec1.push_back(m_k_1_mers_int_11[i+1]);  
                        m_vec2.push_back(i+1);             
                        m_vec3.push_back(counter);
                        counter++;
                    }
                }
            }
            
            m_no_vertices = counter;
            m_adj = new std::vector<int>[m_no_vertices];
            
            int i = 0;
            int search = 0;
            for(i = 0; i < (m_size-(m_kmer-1))*2; i+=2){

                if (elementExistsInVector(m_vec2,i,search)){
                    search++;
                    std::cout <<"i = "<< i << "\n";
                    if(m_k_1_mers_int_11[i] == m_k_1_mers_int_11[i+1]){
                        addEdge(m_vec3[search-1], m_vec3[search-1]);
                        continue;    
                    }
                    addEdge(m_vec3[search-1], m_vec3[search]);
                    std::cout << "Check:  "<<m_k_1_mers_int_11[i] << "     " << m_k_1_mers_int_11[i+1] << "\n";
                    if (m_k_1_mers_int_11[i] != m_k_1_mers_int_11[i+1]){
                        std::cout << "double way : \n";
                        addEdge(m_vec3[search], m_vec3[search-1]);
                    }
                    
                }
            }
        
        }
        // Print the graph
        void printGraph() {
            for (int d = 0; d < m_no_vertices; ++d) {
                std::cout << "\n Vertex " << m_vec1[d] << ":";
                for (auto x : m_adj[d]){
                    std::cout << "-> " << m_vec1[x];
                }
                std::cout<<"\n";
            }
        }
        // Destructor
        ~De_Bruijn_Graph() {
            // Deallocate m_k_1_mers11
            for (int i = 0; i < (m_size - (m_kmer - 1)) * 2; i++) {
                delete[] m_k_1_mers11[i];
            }
            delete[] m_k_1_mers11;

            // Deallocate m_k_1_mers_int_11
            delete[] m_k_1_mers_int_11;

            // Deallocate m_adj (the dynamically allocated array of vectors)
            delete[] m_adj;
        }
    
};

int main(){
    int length, size_reads, k;
    std::cout << "Input the length of the sequence: ";
    std::cin >> length;
    std:: cout << "Input the reads size: ";
    std::cin >> size_reads;
    std::cout << "Input the k for k-mers you want: ";
    std::cin >> k;
    
    De_Bruijn_Graph g1(length, size_reads, k);

    g1.create_the_graph();
    g1.printGraph();
    
}


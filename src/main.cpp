#include "../include/DeBruijnGraph.h"

std::unordered_map<char, Amino> aminoMap_1 = {
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

int main() {
    
    int lines, line_size, kmer;
    std::cout << "Input the number of lines: ";
    std::cin >> lines;

   
    //The vector aminos stores our input lines.
    std::vector<std::vector<Amino>> aminos(lines, std::vector<Amino>());

    for (int i = 0; i < lines; i++){
        
        std::cout << "Input a sequence of "<<line_size<<" amino acids as letters (A, B, ..., T): ";
        std::string userInput;
        std::cin >> userInput;
        for (char c : userInput) {
            if (aminoMap_1.find(c) != aminoMap_1.end()) 
                aminos[i].push_back(aminoMap_1[c]);
            
            else 
                std::cout << "Invalid amino acid: " << c << "\n";
        } 
    }
    std::cout << "size = " << aminos.size() << "\n";
    
    int all_aminos=0;
    //PRINTING THE aminos VECTOR.  
    std::cout << "\nAminos\n";
    for (int i = 0; i < lines; i++) {
        for (Amino amino : aminos[i]) {
            char aminoChar = static_cast<char>(static_cast<int>(amino) + 'A' - 1);
            std::cout << aminoChar;
            all_aminos++;
        }
        
        std::cout << "\n";
    } 
    std::cout << "All aminos : " << all_aminos << "\n";

    std::cout << "Input the k-mer you want: ";
    std::cin >> kmer;
    
    // I HAVE TO CORRECT THE K_MER_VEC_SIZE. SOSOSOSOSO
    int k_mer_vec_size = all_aminos  - (kmer-1) * lines; 
    std::vector<std::vector<Amino>> k_mer_vec(k_mer_vec_size, std::vector<Amino>(kmer));

    int counter = 0 ; // counter is used to identify the input line.
    int index = 0;    // index indicates the line of the k_mer array.
                        // index = numbers of total k_mers = (m_line_size - (m_kmer -1))* m_lines
    std::cout << "HERE \n";
    while (counter < lines){
        //k_mer_vec[lines].reserve(aminos[counter].size() - (kmer -1));
        for (int i = 0; i < aminos[counter].size() - (kmer -1); i++){

            for(int j = 0; j < kmer; j++){
                k_mer_vec[index][j] = aminos[counter][i+j];
            }
            index++;
        }
        counter++;
    }
    std::cout << "End of checkkk \n";
    
    
    De_Bruijn_Graph g1(lines, kmer, k_mer_vec);
    g1.print_graph();

    if (g1.find_euler(g1.m_vec1[0]) > 0) {
        g1.create_euler_path_cycle();
        g1.print_euler_path();
    }
    return 0;
}

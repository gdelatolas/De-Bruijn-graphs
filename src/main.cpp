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

    std::cout << "Input the size of each line: ";
    std::cin >> line_size;

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
    
    
    //PRINTING THE aminos VECTOR.  
    std::cout << "\nAminos\n";
    for (int i = 0; i < lines; i++) {
        for (Amino amino : aminos[i]) {
            char aminoChar = static_cast<char>(static_cast<int>(amino) + 'A' - 1);
            std::cout << aminoChar;
        }
        std::cout << "\n";
    } 
        
    std::cout << "Input the k-mer you want: ";
    std::cin >> kmer;
    



    //=============================================================================================================

    
    De_Bruijn_Graph g1(lines, line_size, kmer, aminos);
    g1.create_the_graph();
    g1.print_graph();

    if (g1.find_euler(g1.m_vec1[0]) > 0) {
        g1.create_euler_path_cycle();
        g1.print_euler_path();
    }
    return 0;
}

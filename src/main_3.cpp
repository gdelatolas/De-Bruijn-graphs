#include "../include/DeBruijnGraph.h"

int main() {
    De_Bruijn_Graph g1;
    g1.create_the_graph();
    g1.printGraph();

    if (g1.find_Euler(g1.m_vec1[0]) > 0) {
        g1.create_Euler_Path_Cycle();
        g1.print_Euler_path();
    }

    return 0;
}

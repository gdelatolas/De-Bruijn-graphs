#include "../include/DeBruijnGraph.h"

int main() {
    De_Bruijn_Graph g1;
    g1.create_the_graph();
    g1.print_graph();

    if (g1.find_euler(g1.m_vec1[0]) > 0) {
        g1.create_euler_path_cycle();
        g1.print_euler_path();
    }

    return 0;
}

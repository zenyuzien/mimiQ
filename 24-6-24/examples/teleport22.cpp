#include <bits/stdc++.h>
// please replace this path with yours
std::string ABSPATH = "C:/Users/rushi/Downloads/mimiq/" ;
#include "../mimiq.h"

int main()
{   
    Qcircuit qc("q_teleport_ptcl_cycle");
    
    qc.printVector();
    auto res = qc.simulate(1024);
    
    res.get_counts();
    qc.get_pdf();

    return 0 ;
}

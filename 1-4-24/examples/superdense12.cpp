#include <bits/stdc++.h>
// please replace this path with yours
std::string ABSPATH = "C:/Users/rushi/Downloads/mimiq/" ;
#include "..\mimiq.h"


int main()
{   
    Qcircuit qc("q_superdense_ptcl"); 
    auto res = qc.simulate(1);
    
    res.print_state();
    res.get_probs();
    res.get_counts();
    get_pdf();
    return 0 ;
}
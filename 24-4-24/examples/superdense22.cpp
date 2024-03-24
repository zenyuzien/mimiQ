#include <bits/stdc++.h>
#include "..\mimiq.h"

int main()
{   
    Qcircuit qc("q_superdense_ptcl_rand"); 
    auto res = qc.simulate(1024);
    
    res.print_state();
    res.get_counts();
    return 0 ;
}
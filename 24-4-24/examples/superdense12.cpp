#include <bits/stdc++.h>
#include "..\mimiq.h"


int main()
{   
    Qcircuit qc("q_superdense_ptcl"); 
    auto res = qc.simulate(1);
    
    res.print_state();
    res.get_counts();
    return 0 ;
}
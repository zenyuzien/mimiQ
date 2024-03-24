#include <bits/stdc++.h>
#include "../mimiq.h"

int main()
{   
    Qcircuit qc("q_teleport_ptcl");
    auto res = qc.simulate(1024);
    
    res.print_state();
    res.get_counts();

    return 0 ;
}

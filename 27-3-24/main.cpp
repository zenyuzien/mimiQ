#include <bits/stdc++.h>

// please replace this path with yours
std::string ABSPATH = "C:/Users/rushi/Downloads/mimiq/" ;
#include "mimiq.h"

int main()
{   
    Qcircuit qc(3,3);
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2);
    qc.applyControlledGate("z","q0",1);
    qc.applyControlledGate("z","q0",2);
    qc.applyGate("h",1);
    qc.applyGate("h",2);
    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);
    auto res = qc.simulate();
    res.print_state();
    res.get_probs();
    res.get_counts();
    get_pdf();
    return 0 ;
}
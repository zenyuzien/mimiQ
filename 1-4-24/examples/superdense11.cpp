#include <bits/stdc++.h>
// please replace this path with yours
std::string ABSPATH = "C:/Users/rushi/Downloads/mimiq/" ;
#include "..\mimiq.h"

int main()
{   
    Qcircuit qc(2,4); 
    qc.applyGate("h",0);
    qc.applyControlledGate("x","q0",1);
    
    // c0-> c , c1 ->d
    qc.setCbit0(0); 
    qc.setCbit1(1); 

    qc.applyControlledGate("z","c1",0); 
    qc.applyControlledGate("x","c0",0);

    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);

    qc.measure(0,3);
    qc.measure(1,2);
    
    qc.printVector();
    auto res = qc.simulate(1,"q_superdense_ptcl");
    
    res.print_state();
    res.get_counts();
    res.get_probs();
    get_pdf();
    return 0 ;
}
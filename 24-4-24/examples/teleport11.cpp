#include <bits/stdc++.h>
#include "..\mimiq.h"

int main()
{   
    Qcircuit qc(3,3); 
    qc.applyGate("u", 0, {0.3,0.2,0.1}); // U(0.3,0.2,0.1) on qbit 0 
    qc.applyGate("h",1);
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);  
    qc.measure(0,0); // measure 0th qbit on 0th classical bit // alice measurement
    qc.measure(1,1); // alice measurement
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("z","q0",2);
    qc.measure(2,2,0); // bob measurement 
     
    qc.printVector();
    auto res = qc.simulate(1024,"q_teleport_ptcl");
    
    res.print_state();
    res.get_counts();

 

    return 0 ;
}
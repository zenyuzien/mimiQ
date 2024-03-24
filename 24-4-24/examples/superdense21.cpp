#include <bits/stdc++.h>
#include "..\mimiq.h"


int main()
{   
    Qcircuit qc(3,4); //3rd qbit for randomising the classical bits initially

     // c0-> c , c1 ->d
    qc.applyGate("h",2);
    qc.measure(2,0);
    qc.measure(2,1);

    qc.applyGate("h",0);
    qc.applyControlledGate("x","q0",1);
    
    qc.applyControlledGate("z","c1",0); 
    qc.applyControlledGate("x","c0",0);

    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);

    qc.measure(0,3);
    qc.measure(1,2);
    
    qc.printVector();
    auto res = qc.simulate(1024,"q_superdense_ptcl_rand");
    
    res.print_state();
    res.get_counts();
    return 0 ;
}
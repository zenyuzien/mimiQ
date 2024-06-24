#include <bits/stdc++.h>
// please replace this path with yours
std::string ABSPATH = "C:/Users/rushi/Downloads/mimiq/" ;
#include "..\mimiq.h"

//g++.exe (GCC) 11.2.0

int main()
{   
    Qcircuit qc(3,3); 
    qc.applyGate("u", 0, {0.3,0.2,0.1}); // U(0.3,0.2,0.1) on qbit 0 
    qc.applyGate("h",1);
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);  
    qc.print_state();
    qc.measure(0,0,1); // measure 0th qbit on 0th classical bit // alice measurement
    qc.print_state();
    qc.measure(1,1,1); // alice measurement
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("z","q0",2);
    qc.measure(2,2,1,0); // bob measurement 
     
    qc.printVector();
    //qc.printQ();
    auto res = qc.simulate(10);//,"q_teleport_ptcl");
    
    qc.print_state();
    //res.print_state();
    res.get_counts();
    //res.get_probs();
    
    //get_pdf();

 

    return 0 ;
}
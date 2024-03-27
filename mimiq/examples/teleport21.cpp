#include <bits/stdc++.h>
// please replace this path with yours
std::string ABSPATH = "C:/Users/rushi/Downloads/mimiq/" ;
#include "..\mimiq.h"

int main()
{   
    Qcircuit qc(3,3); 
    qc.applyGate("u", 0, {6.28,0.408,2.22}); // U(0.3,0.2,0.1) on qbit 0 
    qc.applyGate("h",1);
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);  
    qc.measure(0,0); // measure 0th qbit on 0th classical bit // alice measurement
    qc.measure(1,1); // alice measurement
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("z","q0",2);

    qc.applyGate("u3",2,{-6.28,-2.22,-0.408});
    /* or you can also do 
    qc.applyGate("iu", 2, {6.28,0.408,2.22}); 
    */

    qc.measure(2,2,0); // bob measurement //qc.applyGate("iu", 2, {0.3,0.2,0.1}); // inverse U(0.3,0.2,0.1) on qbit 2 //if uncommented this line, 000 , 010, 100, 110 equal distributions, else, as per U gate 
    qc.printVector();
    auto res = qc.simulate(1024,"q_teleport_ptcl_cycle");
    
    res.print_state();
    res.get_counts();
    res.get_probs();
    get_pdf();
    return 0 ;
}
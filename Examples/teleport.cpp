

#include "../src/structures/structures.h"
#include "../src/gates/gates.h"

Experiment quantum_teleportation_prtcl_var1(mimiqHandler* handler) 
{
    Qcircuit qc(handler, 3,3, "quantum teleportation 1"); 
    
    qc.applyGate("u", 0, {0.3,0.2,0.1}); // U(0.3,0.2,0.1) on qbit 0 
    qc.applyGate("h",1);
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);  
    qc.measure(0,0); // measure 0th qbit on 0th classical bit // alice measurement
    qc.measure(1,1); // alice measurement
    qc.applyControlledGate("x","c1",2);
    qc.applyControlledGate("z","c0",2);
    qc.measure(2,2); // bob measurement 
     
    return qc.simulation();
}

Experiment quantum_teleportation_prtcl_var2(mimiqHandler* handler)
{
    Qcircuit qc(handler, 3,3, "quantum teleportation ibm"); 
    
    qc.applyGate("u", 0, {6.28,0.408,2.22}); // U(0.3,0.2,0.1) on qbit 0 
    qc.applyGate("h",1);
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);  
    qc.measure(0,0); // measure 0th qbit on 0th classical bit // alice measurement
    qc.measure(1,1); // alice measurement
    qc.applyControlledGate("x","c1",2);
    qc.applyControlledGate("z","c0",2);
    qc.applyGate("u3",2,{-6.28,-2.22,-0.408}); // inverse the initial gate
    qc.measure(2,2); // bob measurement 
     
    return qc.simulation();   
}

int main()
{  
    mimiqHandler* handler= new mimiqHandler("/home/zenyuzien/Downloads/"); // path to store report, default is current directory
    result res; 

    res = simulate(handler, quantum_teleportation_prtcl_var1, 1024 );
    res.print_counts();

    res = simulate(handler, quantum_teleportation_prtcl_var2, 1024 );
    res.print_counts();
    handler->clean();

    
    handler->generateReport();
    delete handler;
    return 0 ;
}

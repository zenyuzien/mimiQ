

#include "../src/structures/structures.h"
#include "../src/gates/gates.h"

Experiment quantum_superdense_random(mimiqHandler* handler)
{
    Qcircuit qc(handler, 3,4, "superdense coding - random"); //3rd qbit for randomising the classical bits initially
     // c0-> c , c1 ->d
    qc.applyGate("h",2);                    //init randomiser
    qc.measure(2,0);                        // c0 has random value
    qc.applyControlledGate("x","c0",2);     // reset wire
    qc.applyGate("h",2);                    //init randomiser
    qc.measure(2,1);                        // c1 has random value

    // algorithm
    qc.applyGate("h",0);
    qc.applyControlledGate("x","q0",1);
    qc.applyControlledGate("z","c1",0); 
    qc.applyControlledGate("x","c0",0);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);

    // final meaurements
    qc.measure(0,3);
    qc.measure(1,2);
    return qc.simulation();
}

Experiment quantum_superdense(mimiqHandler* handler)
{
    Qcircuit qc(handler, 2,4, "superdense coding"); 

    // c0-> c , c1 ->d
    qc.applyGate("x",0);
    qc.measure(0,1);
    qc.applyGate("x",0);
    qc.applyGate("h",0);
    qc.applyControlledGate("x","q0",1);
    qc.applyControlledGate("z","c1",0); 
    qc.applyControlledGate("x","c0",0);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);
    qc.measure(0,3);
    qc.measure(1,2);
    
    return qc.simulation();
}


int main()
{  
    mimiqHandler* handler= new mimiqHandler("/home/zenyuzien/Downloads/"); // path to store report, default is current directory
    result res; 

    res = simulate(handler, quantum_superdense, 1024 );
    res.print_counts();

    res = simulate(handler, quantum_superdense_random, 1024 );
    res.print_counts();
    handler->clean();

    
    handler->generateReport();
    delete handler;
    return 0 ;
}

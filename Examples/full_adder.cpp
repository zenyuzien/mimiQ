

#include "../src/structures/structures.h"
#include "../src/gates/gates.h"

Experiment Quantum_full_adder(mimiqHandler* handler)
{
    Qcircuit qc(handler, 4,4, "full adder");
    /*
        q0 -> A 
        q1 -> B
        q2 -> carryin_sumout
        q3 -> carryout
    */    

    // initialize inputs to some values
    // initialize inputs A=1, B=0 and carry_in=1
    qc.applyGate("x",0);
    qc.applyGate("x",2);

    // perform addition
    // decomposition of toffoli q[0], q[1], q[2]
    qc.applyGate("h",2);
    qc.applyControlledGate("x","q1",2);
    qc.applyGate("td",2);
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("t",2);
    qc.applyControlledGate("x","q1",2);
    qc.applyGate("td",2);
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("t",1);
    qc.applyGate("t",2);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",2);
    qc.applyGate("t",0);
    qc.applyGate("td",1);
    qc.applyControlledGate("x","q0",1);

    qc.applyControlledGate("x","q0",1);

    // decomposition of toffoli q[1], q[2], q[3]
    qc.applyGate("h",3);
    qc.applyControlledGate("x","q2",3);
    qc.applyGate("td",3);
    qc.applyControlledGate("x","q1",3);
    qc.applyGate("t",3);
    qc.applyControlledGate("x","q2",3);
    qc.applyGate("td",3);
    qc.applyControlledGate("x","q1",3);
    qc.applyGate("t",2);
    qc.applyGate("t",3);
    qc.applyControlledGate("x","q1",2);
    qc.applyGate("h",3);
    qc.applyGate("t",1);
    qc.applyGate("td",2);
    qc.applyControlledGate("x","q1",2); 

    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("x","q0",1);

    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);
    qc.measure(3,3);

    return qc.simulation();
}

int main()
{  
    mimiqHandler* handler= new mimiqHandler("/home/zenyuzien/Downloads/"); // path to store report, default is current directory
    result res; 

    res = simulate(handler, Quantum_full_adder, 1000 );
    res.print_counts();

    handler->generateReport();
    delete handler;
    return 0 ;
}

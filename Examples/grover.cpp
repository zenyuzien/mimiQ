

#include "../src/structures/structures.h"
#include "../src/gates/gates.h"

Experiment Grover(mimiqHandler* handler)
{
    // grover's algorithm for searching the decimal number 6 in a database of size 2^3
    Qcircuit qc(handler,3,3,"grover");
    
    
    qc.applyGate("h",2);
    qc.applyGate("h",1);
    qc.applyGate("h",0);

    // oracle first run
    qc.applyGate("x",0);
    qc.applyGate("h",2);


    //qc.applyControlledGate("x","q0,q1",2); // circuit can't be drawn for now, so going with decomp
    //decomposition of toffoli q0,q1,q2
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

    qc.applyGate("h",2);
    qc.applyGate("x",0);

    // diffusion operator
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2);    
    qc.applyGate("x",0);
    qc.applyGate("x",1);
    qc.applyGate("x",2);  
    qc.applyGate("h",2);    

    // decomposition of toffoli q0,q1,q2
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

    qc.applyGate("h",2);
    qc.applyGate("x",1);
    qc.applyGate("x",0);
    qc.applyGate("x",2);     
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2); 

    // oracle second run
    qc.applyGate("x",0);
    qc.applyGate("h",2);

    // decomposition of toffoli q0,q1,q2
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

    qc.applyGate("h",2);
    qc.applyGate("x",0);

    // diffusion operator
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2);    
    qc.applyGate("x",0);
    qc.applyGate("x",1);
    qc.applyGate("x",2);  
    qc.applyGate("h",2);  

    // decomposition of toffoli q0,q1,q2
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

    qc.applyGate("h",2);
    qc.applyGate("x",1);
    qc.applyGate("x",0);
    qc.applyGate("x",2);     
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2); 

    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);

    return qc.simulation();
}


int main()
{  
    mimiqHandler* handler= new mimiqHandler("/home/zenyuzien/Downloads/"); // path to store report, default is current directory
    result res; 

    res = simulate(handler, Grover, 1000 );
    res.print_counts();

    handler->generateReport();
    delete handler;
    return 0 ;
}

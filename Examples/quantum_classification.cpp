

#include "../src/structures/structures.h"
#include "../src/gates/gates.h"


Experiment quantum_classification(mimiqHandler* handler)
{
    Qcircuit qc(handler,4,4,"quantum classification");
    qc.applyGate("h",0);
    qc.applyGate("h",1);

    /* Encode new data point */
    // encode the new data point, this implements a cRy(omega)
    qc.applyControlledGate("x","q1",2);
    qc.applyGate("ry",2,{0.1105});// (Ry q[2], -omega/2)
    qc.applyControlledGate("x","q1",2);
    qc.applyGate("ry",2,{-0.1105});// (Ry q[2], omega/2)
    qc.applyGate("x",1);


    /* Encode first training point */
    // encode the first data point, this implements a ccRy(theta)

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

    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{0});// (Ry q[2], theta/4)
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{0});// (Ry q[2], -theta/4)

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

    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{0});// (Ry q[2], theta/4)
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{0});// (Ry q[2], -theta/4)
    qc.applyGate("x",0);

    /* .Encode second training point */
    // encode the second data point, this implements a ccRy(phi)

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

    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{-1.511125});// (Ry q[2], phi/4)
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{1.511125});// (Ry q[2], -phi/4)

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

    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{1.511125});// (Ry q[2], -phi/4)
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{-1.511125});// (Ry q[2], phi/4)

    /* Labels */
    // encode the labels
    qc.applyControlledGate("x","q0",3);

    /* Algorithm */
    // The actual algorithm
    qc.applyGate("h",1);

    qc.measure(1,1);
    qc.measure(3,3);
    return qc.simulation();
}
int main()
{  
    mimiqHandler* handler= new mimiqHandler("/home/zenyuzien/Downloads/"); // path to store report, default is current directory
    result res; 

    res = simulate(handler, quantum_classification, 1000 );
    res.print_counts();

    handler->generateReport();
    delete handler;
    return 0 ;
}

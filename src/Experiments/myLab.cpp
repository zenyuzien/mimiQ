
#include "myLab.h"
#include <array>

#include "../structures/structures.h"
#include "../gates/gates.h"


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

    // decomposition of toffoli q0,q1,q2
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

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
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

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
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

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
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

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

Experiment Grover2(mimiqHandler* handler)
{
    Qcircuit qc(handler,3,3,"grover by circuit");
    
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2); 
    qc.applyGate("x",0);

    qc.applyGate("h",2); 
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

    qc.applyGate("x",0);
    qc.applyGate("h",2); 
   
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2); 
    qc.applyGate("x",1);
    qc.applyGate("x",0);
    qc.applyGate("x",2);  

    qc.applyGate("h",2); 
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

    qc.applyGate("h",2); 

    qc.applyGate("x",1);
    qc.applyGate("x",0);
    qc.applyGate("x",2);     
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2);

    qc.applyGate("x",0);

    qc.applyGate("h",2); 
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */


    qc.applyGate("x",0);
    qc.applyGate("h",2); 
 
    qc.applyGate("h",0);
    qc.applyGate("h",1);
    qc.applyGate("h",2);
    qc.applyGate("x",1);
    qc.applyGate("x",0);
    qc.applyGate("x",2);    

    qc.applyGate("h",2); 
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

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

Experiment BB84_QKD(mimiqHandler* handler)
{
    Qcircuit qc(handler,2,6, "BB84-ptcl QKD");
    // alice
    qc.applyGate("h",1);                                                // init randomiser
    qc.measure(1,0);                                         // c0 has a random -> alice bit 
    qc.applyControlledGate("x","c0",1);                                 // reset randomiser
    qc.applyControlledGate("x","c0",0);                                 // alice bit applied
    qc.applyGate("h",1);                                                // init randomiser
    qc.measure(1,1);                                         // c1 has a random -> alice basis
    qc.applyControlledGate("x","c1",1);                                 // reset randomiser
    qc.applyControlledGate("h","c1",0);                                 // alice basis applied
    // eve 
    qc.applyGate("h",1);                                                // init randomiser
    qc.measure(1,4);                                         // c4 has a random value -> eve presence
    qc.applyControlledGate("x","c4",1);                                 // reset randomiser
    if( qc.accessCreg(4) )
    {
        qc.applyGate("h",1);                                            // init randomiser
        qc.measure(1,5);                                     // c5 has random value -> eve bit
        qc.applyControlledGate("x","c5",1);                             // reset randomiser
        qc.applyControlledGate("x","c5",0);                             // eve bit applied
        qc.applyGate("h",1);                                            // init randomiser
        qc.measure(1,5);                                     // c5 has random value -> eve basis
        qc.applyControlledGate("x","c5",1);                             // reset randomiser
        qc.applyControlledGate("h","c5",0);                             // eve basis applied
    }
    // bob 
    qc.applyGate("h",1);                                                // init randomiser
    qc.measure(1,2);                                         // c2 has randomvalue -> bob basis
    qc.measure(0,3,0,qc.accessCreg(2)); // c3 has bob's measurement -> bob bit
    return qc.simulation();
}

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
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

    qc.applyControlledGate("x","q0",1);

    // decomposition of toffoli q[1], q[2], q[3]
    qc.applyControlledGate("x","q1,q2",3);
    /*
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
    qc.applyControlledGate("x","q1",2); */

    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("x","q0",1);

    qc.measure(0,0);
    qc.measure(1,1);
    qc.measure(2,2);
    qc.measure(3,3);

    return qc.simulation();
}

Experiment trial(mimiqHandler* handler)
{
    Qcircuit qc(handler,3,3,"trial");

    //1 
    qc.applyGate("u3",0,{1.111111,1.111111,1.1111111});//5
    qc.applyGate("u3",0,{1.111111,1.111111,1.1111111});//9
    qc.applyGate("h",0);//10
    qc.applyGate("h",0);    //11
    qc.applyGate("h",0);//12
    qc.applyGate("h",0);//13
    qc.applyGate("u3",0,{1.111111,1.111111,1.1111111});//17
    qc.applyGate("u3",0,{1.111111,1.111111,1.1111111});//21

    qc.applyGate("h",1);
    qc.applyGate("h",1);
    qc.applyGate("u1",1,{1.11111});
    qc.applyGate("u2",1,{1.23456,7.89012});
    qc.applyGate("u3",1,{1.111111,1.111111,1.1111111});
    qc.applyGate("h",1);
    qc.applyGate("h",1);    
    qc.applyGate("h",1);
    qc.applyGate("u1",1,{1.11111});
    qc.applyGate("u2",1,{1.23456,7.89012});
    qc.applyGate("u3",1,{1.111111,1.111111,1.1111111});
    qc.applyGate("h",1);
    qc.applyGate("u1",1,{1.11111});
    qc.applyGate("u2",1,{1.23456,7.89012});
    qc.applyGate("u3",1,{1.111111,1.111111,1.1111111});
    qc.applyGate("h",1);    
    qc.applyGate("h",1);
    qc.applyGate("u1",1,{1.11111});
    qc.applyGate("u2",1,{1.23456,7.89012});
    qc.applyGate("u3",1,{1.111111,1.111111,1.1111111});

    qc.applyGate("u1",2,{1.11111});//3
    qc.applyGate("u2",2,{1.23456,7.89012});//6
    qc.applyGate("u1",2,{1.11111});//8
    qc.applyGate("u2",2,{1.23456,7.89012});//11
    qc.applyGate("u1",2,{1.11111});//13
    qc.applyGate("u2",2,{1.23456,7.89012});//16
    qc.applyGate("u1",2,{1.11111});
    qc.applyGate("u2",2,{1.23456,7.89012});

    return qc.simulation();
}

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
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{0});// (Ry q[2], theta/4)
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{0});// (Ry q[2], -theta/4)

    // decomposition of toffoli q[0], q[1], q[2]
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{0});// (Ry q[2], theta/4)
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{0});// (Ry q[2], -theta/4)
    qc.applyGate("x",0);

    /* .Encode second training point */
    // encode the second data point, this implements a ccRy(phi)

    // decomposition of toffoli q[0], q[1], q[2]
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{-1.511125});// (Ry q[2], phi/4)
    qc.applyControlledGate("x","q0",2);
    qc.applyGate("ry",2,{1.511125});// (Ry q[2], -phi/4)

    // decomposition of toffoli q[0], q[1], q[2]
    qc.applyControlledGate("x","q0,q1",2);
    /*
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
    */

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
#include <bits/stdc++.h>
#include "mimiq.h"

Experiment quantum_teleportation_prtcl_var1(mimiqHandler* handler) 
{
    
    Qcircuit qc(handler, 3,3); 
    
    qc.applyGate("u", 0, {0.3,0.2,0.1}); // U(0.3,0.2,0.1) on qbit 0 
    qc.applyGate("h",1);
    qc.applyControlledGate("x","q1",2);
    qc.applyControlledGate("x","q0",1);
    qc.applyGate("h",0);  
    //qc.print_state();
    qc.measure(0,0); // measure 0th qbit on 0th classical bit // alice measurement
   // qc.print_state();
   qc.applyControlledGate("h","c0",0);
    qc.measure(1,1); // alice measurement
    qc.applyControlledGate("x","c1",2);
    qc.applyControlledGate("z","c0",2);
    qc.measure(2,2); // bob measurement 
     
    return qc.simulation();
}


// the core principles compactible to QASM but a little spice 


Experiment BB98_QKD(mimiqHandler* handler)
{
    Qcircuit qc(handler,2,6);
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

void myLab(mimiqHandler* handler)
{
    int common_basis = 0 ,leaked_bits = 0, eve_attacks = 0 , checker=0;
    std::string key = "";

    for(int j = 1 ; j < 100 ; j++)
    {
        result res = simulate(handler, BB98_QKD);
        if(res.creg[1] == res.creg[2])
        {
            common_basis++;
            if(res.creg[0] != res.creg[3]) leaked_bits++;
            else 
            key += (res.creg[0] ? "1":"0" );
            if(res.creg[4])
            {
                checker++;
                eve_attacks++;
            }
            // ASSERT(checker == leaked_bits)
                
        }
    }
    std::cout<< "common bases: "<< common_basis << " eveattks: "<< eve_attacks << " differecnes1: " << leaked_bits << std::endl;
    std::cout<<"\nFINAL KEY: "<< key <<" length: "<< key.length()<< std::endl;
    handler->writeInPdf("hi ! ");
    handler->writeInPdf("hey ");
}

Experiment madman(mimiqHandler* h)
{
    Qcircuit qc(h,1,1);
    for(int i = 0 ; i < 10 ; i++)
    qc.applyGate("u",0,{0.3,0.2,0.1});
    qc.applyGate("h",0);
    return qc.simulation();
}

void tmpLab(mimiqHandler* handler)
{
    result res = simulate(handler, madman,10);
    res.print_counts();
}

int main()
{   
    mimiqHandler* myhandler= new mimiqHandler("C:/Users/rushi/Downloads/mimiq/"); // path to store report, default is current directory
    
    tmpLab(myhandler);
    
    myhandler->generateReport();
    delete myhandler;
    return 0 ;
}

  






/*
Experiment trial3(mimiqHandler* handler)
{
    
    Qcircuit qc(handler,1,5);
    
    // alice
    if(std::rand()%2)
    {
        qc.applyGate("x",0);
        qc.setCbit1(0);
    }
    if(std::rand()%2)
    {
        qc.applyGate("h",0);
        qc.setCbit1(1);
    }

    /// eve attack plan
    if(std::rand()%4 == 1)
    {
        qc.measure(0,4);
        qc.setCbit1(4);
        if(std::rand()%2)
        {
            qc.applyGate("x",0);
        }
        if(std::rand()%2)
        {
            qc.applyGate("h",0);
        }

    }

    //bob 
    if(std::rand()%2)
    {
        qc.setCbit1(2);
      //  std::cout<<"basis TRUE";
    }
    //else std::cout<<"basis FALSE";
    qc.measure(0,3,0,qc.accessCreg(2));
    
    if( qc.accessCreg(1) == qc.accessCreg(2) )
    {
        if( qc.accessCreg(0) != qc.accessCreg(3) )
        {
            if(qc.accessCreg(4)== 0)
            std::cout<<" LOGIC MISTAKE ------------------------------"<< qc.accessCreg(0)<< " " << qc.accessCreg(1) 
            << " "<< qc.accessCreg(2) << " " << qc.accessCreg(3) <<" \n";
        }
    }
    
    

    return qc.simulation();
    
}

Experiment trial3_v2(mimiqHandler* handler)
{
    
    Qcircuit qc(handler,2,6);
    
    // alice
    qc.applyGate("h",1);
    qc.measure(1,5);
    qc.applyControlledGate("x","c5",1);

    if(qc.accessCreg(5))
    {
        qc.applyGate("x",0);
        qc.setCbit1(0);
    }

    qc.applyGate("h",1);
    qc.measure(1,5);
    qc.applyControlledGate("x","c5",1);

    if(qc.accessCreg(5))
    {
        qc.applyGate("h",0);
        qc.setCbit1(1);
    }

    /// eve attack plan
    qc.applyGate("h",1);
    qc.measure(1,5);
    qc.applyControlledGate("x","c5",1);

    if(qc.accessCreg(5))
    {
        qc.applyGate("h",1);
        qc.measure(1,5);
        qc.applyControlledGate("x","c5",1);

        if(qc.accessCreg(5))
        {
            qc.measure(0,4);
            qc.setCbit1(4);

            qc.applyGate("h",1);
            qc.measure(1,5);
            qc.applyControlledGate("x","c5",1);

            if(qc.accessCreg(5))
            {
                qc.applyGate("x",0);
            }

            qc.applyGate("h",1);
            qc.measure(1,5);
            qc.applyControlledGate("x","c5",1);

            if(qc.accessCreg(5))
            {
                qc.applyGate("h",0);
            }
        }
    }

    qc.applyGate("h",1);
    qc.measure(1,5);
    qc.applyControlledGate("x","c5",1);

    //bob 
    if(qc.accessCreg(5))
    {
        qc.setCbit1(2);
      //  std::cout<<"basis TRUE";
    }
    //else std::cout<<"basis FALSE";
    qc.measure(0,3,0,qc.accessCreg(2));
    
    if( qc.accessCreg(1) == qc.accessCreg(2) )
    {
        if( qc.accessCreg(0) != qc.accessCreg(3) )
        {
            if(qc.accessCreg(4)== 0)
            std::cout<<" LOGIC MISTAKE ------------------------------"<< qc.accessCreg(0)<< " " << qc.accessCreg(1) 
            << " "<< qc.accessCreg(2) << " " << qc.accessCreg(3) <<" \n";
        }
    }
    
    

    return qc.simulation();
    
}

*/
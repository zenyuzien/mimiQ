#include <bits/stdc++.h>

typedef void Lab; // semantic
#include "structures/structures.h"
#include "gates/gates.h"
#include "Experiments/myLab.h"

Lab myLab(mimiqHandler* handler)
{
    result res; // a result is needed to hold a simulationn result
    /*
        Quantum teleportation - variation 1 - qiskit example 
    */
    res = simulate(handler,quantum_teleportation_prtcl_var1, 1024 );
    res.print_counts();
    handler->clean(); // clean the handler for every new Experiment
    /*
        Quantum teleportation - variation 2 - IBM quantum learning example
    */
    res = simulate(handler,quantum_teleportation_prtcl_var2, 1024 );
    res.print_counts();
    handler->clean();
    /*
        Quantum superdense coding - IBM quantum learning 1st example 
    */
    res = simulate(handler,quantum_superdense, 1024 );
    res.print_counts();
    handler->clean();
    /*
        Quantum superdense coding - IBM quantum learning 2nd example 
    */
    res = simulate(handler,quantum_superdense_random, 1024 );
    res.print_counts();
    handler->clean();

    /*
        Quantum key distribution - BB98 protocol - own implementation
    */
    int common_basis = 0 ,leaked_bits = 0, eve_attacks = 0 , checker=0;
    std::string key = "";
    int keySize = 50 ;

    while(keySize)
    {
        res = simulate(handler, BB98_QKD);
        //handler->clean();
        if(res.creg[1] == res.creg[2])
        {
            common_basis++;

            if(res.creg[0] != res.creg[3]) 
                leaked_bits++;

            else 
            {
                key += (res.creg[0] ? "1":"0" );
                keySize--;
            }

            if(res.creg[4])
            {
                checker++;
                eve_attacks++;
            }
            // ASSERT(checker == leaked_bits)
        }
    }
    std::cout<< "common bases: "<< common_basis << " eveattacks: "<< eve_attacks << " detected eve attcks: " << leaked_bits << std::endl;
    std::cout<<"\nFINAL KEY: "<< key <<" length: "<< key.length()<< std::endl;
    handler->writeInPdf("Resultant key: " + key );
    handler->clean();
    
    /*
        Grovers search (quantum-inspire example)
    */
    res = simulate(handler,Grover2, 1000 );
    res.print_counts();
    
    handler->clean(); // clean the handler for every new Experiment
    
}

int main()
{  
    mimiqHandler* myhandler= new mimiqHandler("/home/zenyuzien/Downloads/"); // path to store report, default is current directory
    
    myLab(myhandler);
    
    myhandler->generateReport();
    delete myhandler;
    return 0 ;
}

// tasks for future:
/*
1. fix state.print()
2. toffoli, drawing of it
3. circuit into multiple
4. more examples
5. parallelism
6. cqasm code gen
7. openqasm code gen
8. add more features from net
*/

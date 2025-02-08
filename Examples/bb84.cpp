#include "../src/structures/structures.h"
#include "../src/gates/gates.h"

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


int main()
{
    
    mimiqHandler* handler= new mimiqHandler("/home/zenyuzien/Downloads/"); // path to store report, default is current directory
    result res; 

    int common_basis = 0 ,leaked_bits = 0, eve_attacks = 0 , checker=0;
    std::string key = "";
    int keySize = 50 ;

    while(keySize)
    {
        res = simulate(handler, BB84_QKD);
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
        }
    }
    std::cout<< "common bases: "<< common_basis << " eveattacks: "<< eve_attacks << " detected eve attcks: " << leaked_bits << std::endl;
    std::cout<<"\nFINAL KEY: "<< key <<" length: "<< key.length()<< std::endl;
    handler->writeInPdf("Resultant key: " + key );
    handler->clean();

    handler->generateReport();
    delete handler;
    return 0;

}
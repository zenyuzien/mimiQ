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

Experiment BB98_demonstration(mimiqHandler* handler)
{
    //TODO - replace rand() with extra qbit,cbit and use H gate for measuring
    // also extra classical bits 
    
    Qcircuit qc(handler, 1,2);
    // 1 qbit, 2 classical bits, 1 for eaves to measure, and 1 for bob

    std::cout<<"_";

    // ALICE SECTION
    std::vector<std::pair<bool,bool>> alice_table ,bob_table; 
    bool alice_basis = std::rand()%2; // alice has the choice to apply X gate before sending the qbit
    bool alice_bit = std::rand()%2;
    if(alice_bit == 1)
        qc.applyGate("x",0);
    if(alice_basis == 1)
        qc.applyGate("h",0);
    alice_table.push_back({alice_basis,alice_bit});

    // EAVES SECTION
    bool eves_basis_to_receive = std::rand()%2 ; 
    bool eaves_basis_to_retransmit = std::rand()%2;
    qc.measure(0,1,0, eves_basis_to_receive);
    if(eaves_basis_to_retransmit == 1)
        qc.applyGate("h",0);

    // BOB SECTION
    bool bob_basis = std::rand()%2;
    qc.measure(0,2,0,bob_basis);

    std::cout<<"*";
    return qc.simulation();
}

Experiment trial(mimiqHandler * handler)
{
    Qcircuit qc(handler, 2,5);

    qc.applyGate("h",0);
    qc.measure(0,0);
    qc.applyControlledGate("x","c0",1);
    // now c0 has either 0 or 1
    
    return qc.simulation();
}

void myLab(mimiqHandler* handler)
{
    for(int j = 1 ; j < 2 ; j++)
    {
        result res = simulate(handler, quantum_teleportation_prtcl_var1,1024);
        res.print_counts();
        //std::cout<<"result from a shot: ";
        //for(auto i: res.creg) std::cout<< i;
        //std::cout << std::endl;
    }

}


int main()
{   
    mimiqHandler* myhandler= new mimiqHandler("C:/Users/rushi/Downloads/mimiq/"); // path to store report, default is current directory
    
    myLab(myhandler);
    
    myhandler->generateReport();
    delete myhandler;
    return 0 ;
}

  
    //struct Qcircuit qc(1,1);
    //qc.applyGate("h",0);
    //qc.print_state();
    //auto res  = qc.simulate();
    //auto p = res.state.measureAlong(0,1);
    //std::cout<< p.first.amplitude() << " " << p.second.amplitude() << std::endl;
    //qc.state.print(); // non-order version, empty 
    //qc.state.printprobs(); // need to put into order 
    
    //res.print_counts();

/*

std::vector<std::vector<Coeff>> g1,g2,g3; 
    g3.resize(2, std::vector<Coeff>(2));
    
    g1 = GATES::unitaryGate(PI/2, 0 , PI);
    GATES::printMatrix(g1);

    g2 = GATES::unitaryGate(0,0,PI/2);
    GATES::printMatrix(g2);

 GATES::matrixmultiply2x2(g1,g2,g3);
    GATES::printMatrix(g3);

    state_vector sv(1);
    sv.coeffs[0] = {1,0}; 
    sv.coeffs[1] = {0,0};
sv.print();
    std::vector<Coeff> result (2, { 0, 0 });
    for (int ii = 0; ii < 2; ++ii)
    {
        for (int j = 0; j < 2; ++j)
        result[ii] = result[ii] + (g3[ii][j]* sv.coeffs[j]);
    }
    for(auto i = 0 ; i < 2 ; i++)
    std::cout << "-> " << result[i].real << " " << result[i].complex << std::endl;
    sv.coeffs = result;   
    sv.print();
*/
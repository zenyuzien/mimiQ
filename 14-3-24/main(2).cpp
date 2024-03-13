#include <bits/stdc++.h>
#include "mimiq.h"

/*
int main2()
{
    std::cout << "mimiq.h" << std::endl;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    
    Circuit qc ;
    qc.printQ();
    qc.mimiQ() ;
    qc.sv[qc.nQ -1].print();
    
    return 0 ;
}*/

int main()
{                           
    Qcircuit qc(3,3); 
    /*
    auto m = GATES::unitaryGate(0.3,0.2,0.1) ;
    printMatrix(m);
    state_vector sv ;
    ketZero(sv);
    
    int n = 1 << 1;
    std::vector<Coeff> result(n);
    for (uint32_t i = 0; i < n; ++i) 
        for (uint32_t j = 0; j < n; ++j) 
            result[i] = result[i] + ( m[i][j] * sv.coeffs[j]);
    sv.coeffs = result ; 
    sv.print();*/

    qc.applyGate("u", 0, {0.3,0.2,0.1}); // U(0.3,0.2,0.1) on qbit 0 

    qc.applyGate("h",1);
    qc.applyControlledGate("x",1,2);

    qc.applyControlledGate("x",0,1);
    qc.applyGate("h",0);  
        
    // alice measurement 
    qc.measure(0,0); // measure 0th qbit on 0th classical bit 
    qc.measure(1,1);
        
    qc.applyControlledGate("x",1,2);
    qc.applyControlledGate("z",0,2);
    
    // bob measurement
    qc.measure(2,2,0);
    
    //qc.printVector();
    
    auto res = qc.simulate(1000); // shots = 1024
    
    res.print_state();
    res.get_counts();

    return 0 ;
}


/*
replace the control bit segments of matrix with identity matrix .

for 4x4 matrix when control bit: 1st qbit 

    replace 
    [0][0] [0][1] 
    [1][0] [1][1] 
    with Identity matrix == diagnol matrix of all 1's 

     00      01     10    11 
00 [0][0] [0][1] [0][2] [0][3] 
01 [1][0] [1][1] [1][2] [1][3] 
10 [2][0] [2][1] [2][2] [2][3] 
11 [3][0] [3][1] [3][2] [3][3] 

for 8x8 matrix when control bit is 2nd qbit 

    replace the submatrices with indentity matrices 
    
    [2][2] [2][3]
    [3][2] [3][3]
    with Identity 
    
    [2][6] [2][7]
    [3][6] [3][7]
    with Identity
    
    [6][2] [6][3]
    [7][2] [7][3]
    with Identity
    
    [6][6] [6][7]
    [7][6] [7][7]
    with Identity
    
     000     001    010    011    100    101   110    111
000 [0][0] [0][1] [0][2] [0][3] [0][4] [0][5] [0][6] [0][7] 
001 [1][0] [1][1] [1][2] [1][3] [1][4] [1][5] [1][6] [1][7] 
010 [2][0] [2][1] [2][2] [2][3] [2][4] [2][5] [2][6] [2][7] 
011 [3][0] [3][1] [3][2] [3][3] [3][4] [3][5] [3][6] [3][7] 
100 [4][0] [4][1] [4][2] [4][3] [4][4] [4][5] [4][6] [4][7] 
101 [5][0] [5][1] [5][2] [5][3] [5][4] [5][5] [5][6] [5][7] 
110 [6][0] [6][1] [6][2] [6][3] [6][4] [6][5] [6][6] [6][7] 
111 [7][0] [7][1] [7][2] [7][3] [7][4] [7][5] [7][6] [7][7] 

nah dude, the above logic fails. we are yet to discover the actual logic
*/

/*
int main()
{
    std::cout << "mimiq.h" << std::endl;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    
    struct Qcircuit qc(3,3); 

    
    
    qc.clear(1);
    //qc.state.print();                   //          alpha 0    0      0    beta  0     0  0 
    qc.applyGate("h",1);
    //qc.state.print();                   // 1/root2* alpha 0    0     alpha beta beta   0  0
    qc.applyControlledGate("x",1,2);
    //qc.state.print();                   // 1/root2* alpha 0    0     alpha beta  0     0  beta
        
    qc.applyControlledGate("x",0,1);    // 1/root2* alpha 0    0    alpha   0   beta   beta 0
    //qc.state.print();
    qc.applyGate("h",0);                // 1/2*     alpha beta beta alpha alpha -beta -beta alpha
    //qc.state.print();
        
    // alice measurement 
    qc.measure(0,0); // measure 0th qbit on 0th classical bit 
    qc.measure(1,1);
        
    qc.applyControlledGate("x",1,2);    // 1/2*     alpha beta alpha beta alpha -beta alpha -beta
    //qc.state.print();
    qc.applyControlledGate("z",0,2);    // 1/2*     alpha beta alpha beta alpha beta alpha beta
       
    // bob measurement
    qc.measure(2,2);
        
    qc.simulate(1000);
        
        
    qc.state.print();
    //for(auto i = m.begin() ; i != m.end(); i++)
    //std::cout << i->first << ": "<< i->second << std::endl ;

    return 0 ;
}*/



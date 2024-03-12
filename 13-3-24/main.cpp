
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
    std::cout << "mimiq.h" << std::endl;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    
    struct Qcircuit qc(3); 

    std::map<int,int> m;
    
    for(int shots = 0 ; shots < 1024 ; shots++)
    {
        qc.clear(true);
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
        int C1= qc.prob_qbit(0);
        int C2= qc.prob_qbit(1);
        
        qc.applyControlledGate("x",1,2);    // 1/2*     alpha beta alpha beta alpha -beta alpha -beta
        //qc.state.print();
        qc.applyControlledGate("z",0,2);    // 1/2*     alpha beta alpha beta alpha beta alpha beta
       
        // bob measurement
        int C3 = qc.prob_qbit(2);
        
        int res = ( C1 << 2 ) + (C2 << 1 ) + C3 ;
        m[res]++;
    }
    
    qc.state.print();
    for(auto i = m.begin() ; i != m.end(); i++)
    std::cout << i->first << ": "<< i->second << std::endl ;

    return 0 ;
}


/* 

currently useless functions:

int measure(const std::vector<double>& probabilities) {
    
    // Use the probabilities as weights for a discrete distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());

    return distribution(gen);
}

double prob(state_vector &sv, int qbit)
{
    if (qbit < 0 || qbit >= sv.n_qbits)
    {
        std::cerr << "Invalid qubit index." << std::endl;
        return -1.0; // indicating an error
    }

    // Calculate the probability of measuring |1⟩ in the specified qubit
    double probability = 0.0;

    // Iterate over the state vector elements and sum up the probabilities
    for (size_t i = 0; i < sv.coeffs.size(); ++i)
    {
        // Check if the qubit is in the |1⟩ state
        if ((i >> qbit) & 1)
        {
            // Accumulate the probability
            probability += sv.coeffs[i].amplitude();
        }
    }

    return probability;
}

void distribute( std::vector<double>& probabilities, const state_vector& sv)
{
    double normalization_factor = 0.0;
    probabilities.clear();
    for (const auto& c : sv.coeffs) {
        probabilities.push_back( c.amplitude() );
        normalization_factor +=  c.amplitude()  ;
    }

    // Normalize the probabilities
    for (double& probability : probabilities) {
        probability /= normalization_factor;
    }
}

*/

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




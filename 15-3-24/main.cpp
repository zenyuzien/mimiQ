#include <bits/stdc++.h>
#include "mimiq.h"

int main()
{   
    Qcircuit qc(1,1);
    qc.applyGate("iu1",0,{0.3});
    qc.measure(0,0,0); // bob measurement //qc.applyGate("iu", 2, {0.3,0.2,0.1}); // inverse U(0.3,0.2,0.1) on qbit 2 //if uncommented this line, 000 , 010, 100, 110 equal distributions, else, as per U gate 
     
    qc.printVector();
    auto res = qc.simulate(1024);
    
    res.print_state();
    res.get_counts();

    return 0 ;
}

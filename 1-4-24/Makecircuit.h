#ifndef _makecirc
#define _makecirc
std::vector < double > generateUnitVector(int n) {
    std::vector < double > result;
    result.push_back(5);
    // Calculate the unit vector elements
    double sumOfSquares = 0.0;
    for(int i = 1; i <= n; ++i) {
        double value = static_cast < double > (rand()) / RAND_MAX; // Random value between 0 and 1
        result.push_back(value);
        sumOfSquares += value * value;
    }
    // Normalize the vector to make the sum of squares equal to 1
    double normalizationFactor = 1.0 / std::sqrt(sumOfSquares);
    for(int i = 1; i <= n; ++i) {
        result[i] *= normalizationFactor;
    }
    return result;
}
struct Circuit {
    std::vector < std::vector < std::vector < Coeff > > > gord; // gate order.... do I really need this ?
    std::vector < std::vector < std::vector < uint32_t > > > qlines;
    std::vector < std::vector < double > > table;
    std::vector < state_vector > sv;
    //std::vector< uint32_t > table ;
    int nQ, nC; // number of q/classical bits //! replace with uint8_t to get mindfuck results (white exp snippet)
    uint64_t creg; // 64 maximum classical bit lines 
    int gorder; // pointer to gate order 
    Circuit() {
        std::ifstream ip("circ.txt");
        std::string si, sj;
        std::getline(ip, si);
        std::istringstream iss0(si);
        int nS;
        iss0 >> nQ;
        iss0 >> nC;
        iss0 >> nS;
        std::cout << nQ << " " << nC << " " << nS << std::endl;
        int i = 0, j = 0;
        uint32_t utility;
        for(int z = 1; z <= nQ; z++) {
            std::getline(ip, si);
            j = 0;
            std::cout << si << std::endl;
            std::istringstream iss1(si);
            std::vector < std::vector < uint32_t > > otmp;
            while(iss1 >> sj) {
                std::vector < uint32_t > tmp;
                std::istringstream iss2(sj);
                iss2 >> utility;
                iss2.ignore();
                tmp.push_back(utility);
                if(utility == 1 || utility == 3) {
                    iss2 >> utility;
                    tmp.push_back(utility);
                }
                else if(utility == 2) {
                    for(int x = 0; x < 4; x++) {
                        iss2 >> utility;
                        iss2.ignore();
                        tmp.push_back(utility);
                    }
                }
                otmp.push_back(tmp);
            }
            qlines.push_back(otmp);
        }
        while(std::getline(ip, si)) {
            std::istringstream iss3(si);
            std::string token;
            int first;
            iss3 >> first;
            if(first == 5) // random val; 
            {
                iss3.ignore();
                iss3 >> first; // 2nd val -- no.of qbits for random state 
                first = 2 * (1 << first); // for each state -> 2 needed ( a + ib )
                table.push_back(generateUnitVector(first));
            }
            else if(first == 1) // alpha beta 
            {
                iss3.ignore();
                getline(iss3, token, '_');
                double value1 = std::stod(token);
                getline(iss3, token, '_');
                double value2 = std::stod(token);
                getline(iss3, token, '_');
                double value3 = std::stod(token);
                getline(iss3, token, '_');
                double value4 = std::stod(token);
                table.push_back({ 1, value1, value2, value3, value4 });
            }
        }
        std::cout << "table: \n";
        for(int r = 0; r < table.size(); r++) {
            for(int r2 = 0; r2 < table[r].size(); r2++) std::cout << table[r][r2] << " ";
            std::cout << std::endl;
        }
    }
    void printQ() {
        std::cout << "\nnumber of qbits: " << nQ << " number of classical bits: " << nC << " commands: " << qlines[0].size() << std::endl; //<< "no. of table entries: " << nS ;
        for(int i = 0; i < nQ; i++) {
            std::cout << "qbit: " << i << ": ";
            for(int j = 0; j < qlines[i].size(); ++j) {
                for(int k = 0; k < qlines[i][j].size(); ++k) std::cout << qlines[i][j][k] << " ";
                std::cout << " ; ";
            }
            std::cout << std::endl;
        }
    }
    state_vector mimiQ();
};
/*

state_vector Circuit::mimiQ()
{
    sv.resize(nQ);
    for(int i = 0 ; i < nQ; i++)
    ketZero(sv[i]);
    
    //sv[nQ-1].print();
    
    if(nQ > 1 )
    for(int i = 0 ; i < nQ-1; i++) 
    sv[i+1] = sv[i]*sv[i+1];
    
    sv[nQ-1].print();
    
    int svptr = nQ - 1; 
    
    int stage = 0 ;
    while(stage< qlines[0].size())
    {
        for(int i = 0 ; i < nQ ; i++)
        {
            if(qlines[i][stage][0] == 1) // apply gate 
            {
                std::vector<std::vector<Coeff>> g ; 
                if(qlines[i][stage][0] == 1)
                {
                    kroneckerBUFF(GATES::hadamard, g, nQ, i );
                    applyGate(g, sv[nQ-1]);
                }
                else if(qlines[i][stage][0] == 2)
                {
                    kroneckerBUFF(GATES::pauli_X, g, nQ, i );
                    applyGate(g, sv[nQ-1]);
                }
                else if(qlines[i][stage][0] == 3)
                {
                    kroneckerBUFF(GATES::pauli_Y, g, nQ, i );
                    applyGate(g, sv[nQ-1]);
                }
                else if(qlines[i][stage][0] == 4)
                {
                    kroneckerBUFF(GATES::pauli_Z, g, nQ, i );
                    applyGate(g, sv[nQ-1]);
                }
            }
            else if(qlines[i][stage][0] == 2) // conditional
            {
                if(qlines[i][stage][1] == 1) //qbit ( nor classical)
                {
                    int qbit = qlines[i][stage][2] ; // qbit specified ( control bit)
                    
                    if(qlines[i][stage][3] == 2) // apply gate 
                    {
                        if(qlines[i][stage][4] == 1)
                        GATES::controlled_hadamard(sv[nQ-1],qbit,i);
                        else if(qlines[i][stage][4] == 2)
                        GATES::controlled_pauli_X(sv[nQ-1],qbit,i);
                        else if(qlines[i][stage][4] == 4)
                        GATES::controlled_pauli_Z(sv[nQ-1],qbit,i);
                    }
                    
                }
                else //classical bit 
                {
                    ;
                }
            }
            else if(qlines[i][stage][0] == 3) // table
            {
                // checking table entry 
                if(table[qlines[i][stage][1]][0] == 5)
                {
                    // assign random 
                    for(int r = 1 ; r < table[qlines[i][stage][1]].size(); )
                    {
                        sv[nQ-1].coeffs[(r-1)/2].real = table[qlines[i][stage][1]][r] ;
                        std::cout << nQ -1 << " " << (r-1)/2 << " real: " << sv[nQ-1].coeffs[(r-1)/2].real << std::endl ;
                        sv[nQ-1].coeffs[(r-1)/2].complex = table[qlines[i][stage][1]][r+1] ;
                        std::cout << nQ -1 << " " << (r-1)/2 << " complex: " << sv[nQ-1].coeffs[(r-1)/2].complex << std::endl ;
                        r+=2 ;
                    }
                }
            }
            
        }
        //std::cout<<"\nafter stage "<< stage ;
        //sv[nQ-1].print();
        ++stage;
    }
        
    return sv[nQ -1] ;
}*/
#endif
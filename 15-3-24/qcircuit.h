#ifndef _qcirc
#define _qcirc

struct Qcircuit 
{
    int n_qbits , n_cbits, order_ptr; 
    state_vector state, ZERO ; 
    uint64_t c_reg ; 
    
    std::vector< std::vector< int > > ORDER ;
    std::vector<std::array<double, 3>> euler_container ;
    
    void printVector() 
    {
        std::cout << "ORDER: \n";
        for (const auto& row : ORDER) {
            for (int element : row) {
                std::cout << element << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl ; 
    }
    
    Qcircuit( int nQ , int nC = 0 ) : n_qbits( nQ ), n_cbits(nC), c_reg(0) , order_ptr(0) {
        
        std::vector<struct state_vector> sv ; 
        sv.resize(nQ);
        
        for(int i = 0 ; i < nQ; i++)
            ketZero(sv[i]);
        
        if(nQ > 1 )
            for(int i = 0 ; i < nQ-1; i++) 
                sv[i+1] = sv[i]*sv[i+1];
        
        ZERO = state = sv[nQ - 1] ;
        
        ORDER.push_back({4,true}); // false means hard reset 
    } 
    
    Qcircuit( std::string path) : c_reg(0) , order_ptr(0) {
        
        std::vector<struct state_vector> sv ; 
        
        ORDER = downloadCircuit("C:/Users/rushi/Downloads/mimiQ++/ORDERS/" + path );
        sv.resize(n_qbits);
        
        for(int i = 0 ; i < n_qbits; i++)
            ketZero(sv[i]);
        
        if(n_qbits > 1 )
            for(int i = 0 ; i < n_qbits -1; i++) 
                sv[i+1] = sv[i]*sv[i+1];
        
        ZERO = state = sv[n_qbits - 1] ;
        
        //ORDER.push_back({4,true}); // false means hard reset 
    }
    
    
    void applyGate(const std::string& gate , int t_qbit)
    {
        int GATE ;
        
        if(gate == "H" || gate == "h")
            GATE = 1 ;
        else if(gate == "x" || gate == "X")
            GATE = 2 ;
        else if(gate == "z" || gate == "Z")
            GATE = 4 ;
        else if(gate == "y" || gate == "Y")
            GATE = 3 ;
        
        ORDER.push_back({1,t_qbit,GATE});
        return ;
    }
    
    void applyGate(const std::string& gate , int t_qbit, std::vector<double> arr) // overloading applyGate 
    {
        if(gate == "U3" || gate == "u3" || gate == "u" || gate == "U")
        {
            int size = euler_container.size() ;
            euler_container.push_back({arr[0],arr[1],arr[2]});
            ORDER.push_back({1,t_qbit,5,size}); // send index
        }
        else if(gate == "IU3" || gate == "iu3" || gate == "IU" || gate == "iu") // todo: generalise u1 u2 u3 inverse 
        {
            int size = euler_container.size() ;
            euler_container.push_back({arr[0],arr[1],arr[2]});
            ORDER.push_back({1,t_qbit,6,size}); // send index
        }
        else if(gate == "U1" || gate == "u1")
        {
            int size = euler_container.size() ;
            euler_container.push_back({arr[0]});
            ORDER.push_back({1,t_qbit,7,size}); // send index
        }
        else if(gate == "U2" || gate == "u2")
        {
            int size = euler_container.size() ;
            euler_container.push_back({arr[0],arr[1]});
            ORDER.push_back({1,t_qbit,8,size}); // send index
        }
    }
    
    struct result simulate(int shots = 1, std::string path = "")
    {
        if(path != "")
        {
            std::cout<<"-";
            uploadCircuit(ORDER,path);
        }
        std::map<int,int> m ;
        while(shots--)
        {
            for(int i = 0; i < ORDER.size() ; i++)
            {
                if( ORDER[i][0] == 1 ) // applygate 
                {
                    std::vector<std::vector<Coeff>> matrix, buffedmatrix ; 
                    
                    switch (ORDER[i][2]) 
                    {
                        case 1:
                            matrix = GATES::hadamard;
                            break;
                        case 2:
                            matrix = GATES::pauli_X;
                            break;
                        case 3:
                            matrix = GATES::pauli_Y;
                            break;
                        case 4:
                            matrix = GATES::pauli_Z;
                            break;
                        case 5:
                            matrix = GATES::unitaryGate(euler_container[ORDER[i][3]][0], euler_container[ORDER[i][3]][1], euler_container[ORDER[i][3]][2]);
                            break;
                        case 6:
                            matrix = inverse2x2(GATES::unitaryGate(euler_container[ORDER[i][3]][0], euler_container[ORDER[i][3]][1], euler_container[ORDER[i][3]][2]));
                            break;
                        case 7: 
                            matrix = GATES::unitaryGate(0, 0, euler_container[ORDER[i][3]][0]);
                            break;
                        case 8: 
                            matrix = GATES::unitaryGate( acos(0.0), euler_container[ORDER[i][3]][0], euler_container[ORDER[i][3]][1]);
                            break;
                    }

                    kroneckerBUFF(matrix, buffedmatrix , n_qbits, ORDER[i][1] );
            
                    int n = 1 << n_qbits;
                    std::vector<Coeff> result(n);
                    for (uint32_t i = 0; i < n; ++i) 
                        for (uint32_t j = 0; j < n; ++j) 
                            result[i] = result[i] + ( buffedmatrix[i][j] * state.coeffs[j]);
            
                    state.coeffs = result;
                }
                
                else if( ORDER[i][0] == 2) // controlledgate 
                {
                    if( ORDER[i][1] == 1 )
                        GATES::controlled_hadamard(state, ORDER[i][2], ORDER[i][3]);
                    else if( ORDER[i][1] == 2 )
                        GATES::controlled_pauli_X(state, ORDER[i][2], ORDER[i][3]);
                    else if( ORDER[i][1] == 4 ) // 4 - pauli z , 3 - pauli y 
                        GATES::controlled_pauli_Z(state, ORDER[i][2], ORDER[i][3]);
                }
                
                else if( ORDER[i][0] == 3 ) // measure 
                {
                    //std::cout << "measuring " << ORDER[i][1]  <<  " \n" ;
                    int cnd = 1 << n_qbits ;
                    auto tmask  =  1 << (n_qbits-1-ORDER[i][1]);
                    double prob = 0.0 ; 
                    
                    for(int i = 0 ; i < cnd ; i++)
                        if( i & tmask )
                            prob += state.coeffs[i].amplitude();
                    
                    double randomValue = static_cast<double>(rand()) / RAND_MAX;
                    bool res = (randomValue <= prob );
                    
                    if(ORDER[i][3] == 1)
                    {
                        std::cout << "Probability of measuring |1> on qubit " << ORDER[i][1] << ": " << prob << std::endl;
                        std::cout << "Probability of measuring |0> on qubit " << ORDER[i][1] << ": " << 1-prob << std::endl;
                        std::cout << "rand/prob: " << randomValue << "/ "<< prob << " res: " << res << " " << (1 << (n_qbits - 1 - ORDER[i][2])) <<std::endl ;
                       // std::cout << 
                    }
                    
                    if(res)
                    c_reg = c_reg | (1 << (n_qbits - 1 - ORDER[i][2])) ;
                    //std::cout<<"update = "<< c_reg << std::endl;
                    
                    // let us say we measured q0 as 1 , we want to do qbits - 1 - q0 = 2 = 100
                    
                    // std::cout<<tbit<<"-> "<<res<<" creg: "<<c_reg<<std::endl ;
                    // 1 << n represents nth cbit 
                }
                
                else if( ORDER[i][0] == 4) // clear 
                {
                    if(ORDER[i][1] == 0)  //hard reset
                    {
                        n_qbits = n_cbits = order_ptr = 0 ; 
                        state.n_qbits =0 ;
                        state.coeffs.clear();
                        for (auto& inner_vector : ORDER) 
                            inner_vector.clear();
                        ORDER.clear();
                    }
                    else // back to ZERO 
                    {
                        order_ptr = 0 ;
                        state = ZERO ;
                    } 
                    c_reg = 0 ;
                    euler_container.clear();
                }
            }
            
            m[c_reg]++;
        }
        return {n_qbits,m,state} ;
    }
    
    void applyControlledGate( const std::string& gate, int c_qbit, int t_qbit )
    {
        if(gate == "H" || gate == "h")
        ORDER.push_back({2,1,c_qbit,t_qbit});
        else if(gate == "x" || gate == "X")
        ORDER.push_back({2,2,c_qbit,t_qbit});
        else if(gate == "z" || gate == "Z")
        ORDER.push_back({2,4,c_qbit,t_qbit});
    }
    
    void measure( int tbit, int cbit, int toPrint = 0 ) //cbit = classical // TODO: int->false 
    {
        ORDER.push_back({3,tbit, cbit,toPrint});
    }
    
    void uploadCircuit(const std::vector<std::vector<int>>& ORDER, const std::string& path) 
    {
        std::ofstream outfile(path + ".txt");
        outfile << n_qbits <<" "<< n_cbits << " "<< ORDER.size() << std::endl;
        if (outfile.is_open()) 
        {
            for (const auto& inner : ORDER) {
                for (int val : inner) {
                    outfile << val << " ";
                }
                outfile << "\n";
            }
            
            for (const auto& array : euler_container) {
                for (const auto& value : array) {
                    outfile << value << " "; // Write each value followed by a space
                }
                outfile << "\n"; // Write newline after each array
            }
            
            outfile.close();
        } 
        else 
            std::cerr << "Error: Unable to open file " << path << ".txt\n";
            
    }

    std::vector<std::vector<int>> downloadCircuit(const std::string& path) 
    {
        int ordersize ; 
        std::vector<std::vector<int>> result;
        std::ifstream infile(path);
        if (infile.is_open()) 
        {
            std::string line;
            std::getline(infile,line);
            std::istringstream iss0(line);
            iss0 >> n_qbits;
            iss0 >> n_cbits;
            iss0 >> ordersize ; 
            while (ordersize--) 
            {
                std::getline(infile, line);
                std::istringstream iss(line);
                int val;
                std::vector<int> inner;
                
                while (iss >> val) 
                    inner.push_back(val);
                
                result.push_back(inner);
            }
            double theta,lam,psi ;
            while(std::getline(infile, line))
            {
                std::istringstream iss1(line);
                iss1 >> theta ;
                iss1 >> lam ; 
                iss1 >> psi ;
                euler_container.push_back({theta,lam,psi});
            }
            
            infile.close();
        } 
        else
            std::cerr << "Error: Unable to open file " << path << ".txt\n";
            
        return result;
    }

};


#endif 
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
        
        ORDER = downloadCircuit("C:/Users/rushi/Downloads/15-3-24/ORDERS/" + path + ".txt");
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

    void print_state()
    {
        ORDER.push_back({9});
    }
    
    void print_creg()
    {
        ORDER.push_back({11});
    }

    struct result simulate(int shots = 1, std::string path = "", bool describe = false)
    {
        if(describe)std::cout<<"simulation starting with shots: "<<shots<<std::endl ;

        if(path != "")
        {
            if(describe)
            std::cout<<"uploading circuit to "<< path << std::endl ;
            uploadCircuit(ORDER,path);
        }
        std::map<int,int> m ;
        int tmp1,tmp2 ; 
        while(shots--)
        {
            for(int i = 0; i < ORDER.size() ; i++)
            {
                if( ORDER[i][0] == 9 ) // print state vector 
                {
                    if(describe) std::cout<< i << ") printing state vector" <<std::endl;
                    state.print();
                }
                else if( ORDER[i][0] == 11) // print creg 
                {
                    if(describe) std::cout<<i<<") printing creg" <<std::endl;
                    std::cout<< "creg: "<<  c_reg<< std::endl;
                }
                else if( ORDER[i][0] == 1 ) // applygate 
                {
                     if(describe) std::cout<<i<<") applygate" <<std::endl;
                    tmp1 = ORDER[i][1]; // qbit
                    tmp2 = ORDER[i][2]; // gate
                    gater:;
                    std::vector<std::vector<Coeff>> matrix, buffedmatrix ; 
                     if(describe) std::cout<<"to qbit and gate "<< tmp1 << ", "<<tmp2 <<std::endl;
                    switch (tmp2) 
                    {
                        case 1:
                            matrix = GATES::unitaryGate(PI/2, 0, PI);
                            break;
                        case 2:
                            matrix = GATES::unitaryGate(PI, 0, PI);
                            break;
                        case 3:
                            matrix = GATES::unitaryGate(PI, PI/2, PI/2);
                            break;
                        case 4:
                            matrix = GATES::unitaryGate(0, 0, PI);
                            break;
                        case 5:
                            matrix = GATES::unitaryGate(euler_container[ORDER[i][3]][0], euler_container[ORDER[i][3]][1], euler_container[ORDER[i][3]][2]);
                            break;
                        case 6:
                            matrix = GATES::inverse2x2(GATES::unitaryGate(euler_container[ORDER[i][3]][0], euler_container[ORDER[i][3]][1], euler_container[ORDER[i][3]][2]));
                            break;
                        case 7: 
                            matrix = GATES::unitaryGate(0, 0, euler_container[ORDER[i][3]][0]);
                            break;
                        case 8: 
                            matrix = GATES::unitaryGate( PI, euler_container[ORDER[i][3]][0], euler_container[ORDER[i][3]][1]);
                            break;
                    }
                    
                    GATES::kroneckerBUFF(matrix, buffedmatrix , n_qbits, tmp1 );
                    /*if(i==14)
                    {
                        GATES::printMatrix(buffedmatrix); 
                        state.print();
                    }*/
                                    

                    int n = 1 << n_qbits;
                    std::vector<Coeff> result(n,{0,0});
                    for (int ii = 0; ii < n; ++ii) 
                    {
                        for (int j = 0; j < n; ++j) 
                            result[ii] = result[ii] + ( buffedmatrix[ii][j] * state.coeffs[j]);
                        //std::cout<< "intermediete "<<ii<<" "<<result[ii].real<<" "<<result[ii].complex<<std::endl;
                    }
                        

                    state.coeffs = result;
                    /*if(i==14)
                    {
                        state.print();
                    }*/
                }
                
                else if( ORDER[i][0] == 2) // controlledgate 
                {
                     if(describe) std::cout<<i<<") apply controlled gate " <<std::endl;
                    if(ORDER[i][2] == 1) // type: q
                    {
                        if(describe) std::cout<< "type:q with gate "<< ORDER[i][1] <<"control and targets: " << ORDER[i][3] <<" "<< ORDER[i][4]<<std::endl;
                        if( ORDER[i][1] == 1 ) // gate
                            GATES::controlled_hadamard(state, ORDER[i][3], ORDER[i][4]); // qbit, tqbit
                        else if( ORDER[i][1] == 2 )
                            GATES::controlled_pauli_X(state, ORDER[i][3], ORDER[i][4]);
                        else if( ORDER[i][1] == 4 ) // 4 - pauli z , 3 - pauli y 
                            GATES::controlled_pauli_Z(state, ORDER[i][3], ORDER[i][4]);
                    }
                    else if(ORDER[i][2]== 0 )
                    {
                        if(describe) std::cout<<"control by cbit.. creg: "<< c_reg <<" mask: "<<(1 << (n_cbits-1 -ORDER[i][3]))<<std::endl;
                        if (c_reg & (1 << (n_cbits-1 -ORDER[i][3])) )
                        {
                            //std::cout<<"going to gater "<<std::endl ;
                            tmp2 = ORDER[i][1];
                            tmp1 = ORDER[i][4]; 
                            goto gater ; 
                        }
                        else if(describe) std::cout<<" cbit control false \n";
                        
                    }

                }
                
                else if( ORDER[i][0] == 3 ) // measure 
                {
                     if(describe) std::cout<<i<<") measuring " << ORDER[i][1]  <<  " \n" <<std::endl;
                    //std::cout << "measuring " << ORDER[i][1]  <<  " \n" ;
                    int cnd = 1 << n_qbits ; // for 2 qbits , cnd = 4 
                    auto tmask  =  1 << (n_qbits-1-ORDER[i][1]); // qbit: 0 means tmask 2-1-0 = 1 shit to left.
                    double prob = 0.0 ; 
                    
                    for(int i = 0 ; i < cnd ; i++)
                        if( i & tmask )
                            prob += state.coeffs[i].amplitude();
                    
                    double randomValue = static_cast<double>(rand()) / RAND_MAX;
                    bool res = (randomValue <= prob );
                    
                    if(ORDER[i][3] == 1 || describe)
                    {
                        std::cout << "Probability of measuring |1> on qubit " << ORDER[i][1] << ": " << prob << std::endl;
                        std::cout << "Probability of measuring |0> on qubit " << ORDER[i][1] << ": " << 1-prob << std::endl;
                        std::cout << "rand/prob: " << randomValue << "/ "<< prob << " res: " << res << " " << (1 << (n_cbits - 1 - ORDER[i][2])) <<std::endl ;
                       // std::cout << 
                    }
                    
                    if(res)
                    c_reg = c_reg | (1 << (n_cbits - 1 - ORDER[i][2])) ;
                     if(describe) std::cout << "updated creg: "<<c_reg << std::endl;
                    
                    // let us say we measured q0 as 1 , we want to do qbits - 1 - q0 = 2 = 100
                    
                    // std::cout<<tbit<<"-> "<<res<<" creg: "<<c_reg<<std::endl ;
                    // 1 << n represents nth cbit 
                }
                
                else if( ORDER[i][0] == 4) // clear 
                {
                    if(describe) std::cout << i << ") clearing with hardreset: "<< !ORDER[i][1] << std::endl;
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
                else if(ORDER[i][0]== 10) //set c bit 
                {
                    if(describe) std::cout<<i<<") set c bit given bit position: "<<ORDER[i][2] << " which is mask "<< (1 << (n_cbits - 1 - ORDER[i][2])) << std::endl;
                    if(ORDER[i][1] == 1)
                        c_reg = c_reg | (1 << (n_cbits - 1 - ORDER[i][2])); 

                    else if(ORDER[i][1]==2)
                        c_reg = c_reg & (~(1 << (n_cbits - 1 - ORDER[i][2]))); 
                
                    if(describe)
                    std::cout<<"currently c_reg: "<<c_reg<<std::endl;
                }
            }
            
            m[c_reg]++;
        }
        return {n_cbits,m,state} ;
    }
    
    void applyControlledGate( const std::string& gate, const std::string c_qbit_str , int t_qbit )
    {
        // 2    appplycontrolledgate
        // 2    gate 
        // 1    type: q
        // 1    qbit
        // 2    t qbit
        std::istringstream iss(c_qbit_str);
        char cq ; 
        int n ; 
        iss >> cq ; 
        iss >> n ;

        int third  = (cq == 'q') ? 1 : 0 ; 

        if(gate == "H" || gate == "h")
            ORDER.push_back({2,1,third,n,t_qbit});
        else if(gate == "x" || gate == "X")
            ORDER.push_back({2,2,third,n,t_qbit});
        else if(gate == "z" || gate == "Z")
            ORDER.push_back({2,4,third,n,t_qbit});
    }
    
    void measure( int tbit, int cbit, int toPrint = 0 ) //cbit = classical // TODO: int->false 
    {
        ORDER.push_back({3,tbit, cbit,toPrint});
    }

    void setCbit1(int cbit)
    {
        ORDER.push_back({10,1,cbit});
        
    }
    void setCbit0(int cbit)
    {
        ORDER.push_back({10,2,cbit});
        
    }
    
    void uploadCircuit(const std::vector<std::vector<int>>& ORDER, const std::string& path) 
    {
        std::string full = "C:/Users/rushi/Downloads/15-3-24/ORDERS/" + path + ".txt" ; 
        std::ofstream outfile( full);
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
            
            std::cout<<"circuit saved! \n";
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

// @author: ZENyuzien Zapolean

#ifndef mimiq
#define mimiq 

#define root2 std::sqrt(2)

struct GlobalInitializer {
    GlobalInitializer() {
        std::cout<<"mimiQ++ \n";
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
    }
} initializer;

struct Coeff 
{
    double real;
    double complex;

    Coeff() : real(0), complex(0) {}
    Coeff(double r, double c) : real(r), complex(c) {}

    Coeff operator*(const Coeff& other) const {
        Coeff result;
        result.real = real * other.real - complex * other.complex;
        result.complex = real * other.complex + complex * other.real;
        return result;
    }

    Coeff operator+(const Coeff& other) const {
        Coeff result;
        result.real = real + other.real;
        result.complex = complex + other.complex;
        return result;
    }

    Coeff operator-(const Coeff& other) const {
        Coeff result;
        result.real = real - other.real;
        result.complex = complex - other.complex;
        return result;
    }

    double amplitude() const { //amp squared actually
        return (real * real) + (complex * complex);
    }
};

struct state_vector
{
    int n_qbits ;
    std::vector< Coeff > coeffs;
    
    state_vector(): n_qbits(0) {}
    state_vector(int n)
    {
        n_qbits = n ; 
        coeffs.resize(1<<n,{0,0});
    }
    state_vector(const state_vector& sv)
    {
        n_qbits = sv.n_qbits;
        coeffs = sv.coeffs;
    }
    
    state_vector operator*(const state_vector& other) // returns a tensor product of 2 state vectors 
    {
        state_vector sv3; //(*this).n_qbits == n_qbits 
        sv3.n_qbits = n_qbits + other.n_qbits; 
        int cnd1 = 1 << n_qbits , cnd2 = 1 << other.n_qbits ;
        for(int i = 0 ; i < cnd1 ; i++)
        for(int j = 0 ; j < cnd2 ; j++)
        sv3.coeffs.push_back( coeffs[i] * other.coeffs[j] );
        return sv3 ;
    }
    
    void print()
    {
        std::cout<<std::endl<<"number of qbits = "<< n_qbits<< std::endl ;
        int cnd = 1 << n_qbits ;
        for(int i = 0 ; i < cnd ; i++)
        {
            std::string binaryString;
            for (int j = n_qbits - 1; j >= 0; --j) {
                int bit = (i >> j) & 1;
                binaryString += (bit == 0) ? '0' : '1';
            }
            std::cout<<binaryString << " => probability: "<< coeffs[i].amplitude()*100 <<"% ( " << coeffs[i].real << " |0> + "<<coeffs[i].complex << " i |1> )\n";
        }
        std::cout<<std::endl;
    
    }
};

void ketZero(struct state_vector& sv)
{
    sv.n_qbits = 1 ;
    sv.coeffs.push_back( {1,0} );
    sv.coeffs.push_back( {0,0} );
}

void ketONE(struct state_vector& sv)
{
    sv.n_qbits = 1 ;
    sv.coeffs.push_back( {0,0} );
    sv.coeffs.push_back( {1,0} );
}

namespace GATES {// NEED TO CONST !!!!!!!!!!!!!! maybe convert to array as well ?

    std::vector<std::vector<Coeff>> hadamard = {
        {{1 / root2, 0}, {1 / root2, 0}},
        {{1 / root2, 0}, {-1 / root2, 0}}
    };

    std::vector<std::vector<Coeff>> pauli_X = {
        {{0, 0}, {1, 0}},
        {{1, 0}, {0, 0}}
    };

    std::vector<std::vector<Coeff>> pauli_Z = {
        {{1, 0}, {0, 0}},
        {{0, 0}, {-1, 0}}
    };

    std::vector<std::vector<Coeff>> pauli_Y = {
        {{0, 0}, {0, -1}},
        {{0, 1}, {0, 0}}
    };

    std::vector<std::vector<Coeff>> identity = {
        {{1, 0}, {0, 0}},
        {{0, 0}, {1, 0}}
    };
    
    std::vector<std::vector<Coeff>> unitaryGate(double theta, double psi, double lam) 
    {
        return {
            {
                {cos(theta/2), 0},
                {cos(lam)*sin(theta/2), sin(lam)*sin(theta/2)}
            },
            {
                {cos(lam + psi)*sin(theta/2), sin(lam+ psi)*sin(theta/2)},
                {cos(theta/2), 0}
            }
        };
    }

    
    void controlled_pauli_X( struct state_vector& sv, int cbit, int tbit )
    {
        int cnd = 1 << sv.n_qbits , cmask = 1 << (sv.n_qbits-1-cbit) , tmask = 1 << (sv.n_qbits-1-tbit) , other ; 
        //std::cout<< cnd << " " << cmask << " " << tmask << std::endl; 
        std::vector<bool> visit( cnd , false);
        double tmp1, tmp2 ;
        
        for(int i = 0 ; i< cnd ; i++)
        {
            if( (i & cmask) && (visit[ i ]== false) )
            {
                if(i & tmask ) // if tbit 1 -> other = xxxx & 11011 => 0 at tbit 
                other = i & ~(tmask); 
                else 
                other  = i | tmask ; // else, other = xxxx | 00100 => 1 at tbit 

                tmp1 = sv.coeffs[i].real ;
                tmp2 = sv.coeffs[i].complex ;

                sv.coeffs[i].real = sv.coeffs[other].real;
                sv.coeffs[i].complex = sv.coeffs[other].complex;

                sv.coeffs[other].real = tmp1 ;
                sv.coeffs[other].complex = tmp2; 
                
                visit[i] = visit[other]= true ;
            }
        }
    }
    
    void controlled_hadamard(struct state_vector& sv, int cbit, int tbit)
    {
        int cnd = 1 << sv.n_qbits ,
        cmask = 1 << (sv.n_qbits-1-cbit) , 
        tmask = 1 << (sv.n_qbits-1-tbit) ,
        other ;
        struct state_vector tmp(sv.n_qbits);
        
        for(int i = 0 ; i< cnd ; i++)
        {
            if( i & cmask )
            {
                //std::cout << "in "<< i << std::endl ;
                struct Coeff coeff1(1/root2, 0); // 1/root2 
                
                if( i & tmask ) // if target bit is 1 here, x|1>x --> 1/root2 *[ x|0>x - x|1>x ]
                {
                    other = i & ~(tmask); 
                  //  std::cout << "other "<< other << "\n" << i << " - ";
                    //std::cout <<coeff1.real << " * " << sv.coeffs[i].real << " = "<< ((coeff1)*sv.coeffs[i]).real << std::endl ;
                    tmp.coeffs[i] = tmp.coeffs[i] - (coeff1)*sv.coeffs[i];
                    //std::cout << "updated tmp["<<i<<"] " << tmp.coeffs[i].real << std::endl ;
                }
                else // if target bit is 0 here, x|0>x --> 1/root2 *[ x|0>x + x|1>x ]
                {
                    other  = i | tmask ;
                    //std::cout << "other "<< other << "\n" << i << " + ";
                    //std::cout <<coeff1.real << " * " << sv.coeffs[i].real << " = "<< ((coeff1)*sv.coeffs[i]).real  << std::endl ;
                    tmp.coeffs[i] = tmp.coeffs[i] + (coeff1)*sv.coeffs[i];
                    //std::cout << "updated tmp["<<i<<"] " << tmp.coeffs[i].real << std::endl ;
                }
               // std::cout << "other "<< other << "+ " <<coeff1.real << " * " << sv.coeffs[i].real << " = "<< ((coeff1)*sv.coeffs[i]).real  << std::endl ;
                tmp.coeffs[other] = tmp.coeffs[other] + (coeff1)*sv.coeffs[i];
                //std::cout << "updated tmp["<<other<<"] " << tmp.coeffs[i].real << std::endl ;
            }
            else 
            tmp.coeffs[i] = sv.coeffs[i] ; 
        }
        sv = tmp ; 
        return ;
    }
    
    void controlled_pauli_Z( struct state_vector& sv, int cbit, int tbit )
    {
        int cnd = 1 << sv.n_qbits , cmask = 1 << (sv.n_qbits-1-cbit) , tmask = 1 << (sv.n_qbits-1-tbit) ; 
        struct Coeff tmp(-1,0);
        
        for(int i = 0 ; i< cnd ; i++)
            if((i & cmask)&&(i & tmask) )// if cbit and tbit 1 flip sign.
                sv.coeffs[i] = sv.coeffs[i] * tmp ; 
    }
    
}

void printMatrix(const std::vector<std::vector<Coeff>>& matrix)
{
    if(matrix.size() ==0 )
    {
        std::cout << "empty " << std::endl ;
        return ;
    }
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            std::cout << "{" << element.real << ", " << element.complex << "} ";
        }
        std::cout << std::endl;
    }std::cout << std::endl;
}

void kroneckerSUB(const std::vector<std::vector<Coeff>>& m1, const std::vector<std::vector<Coeff>>& m2, std::vector<std::vector<Coeff>>& m3)
{
    //printMatrix(m1);
    //printMatrix(m2);
    
    if(m1.size()== 0)
    {
        //std::cout<<"m1 empty "<< std::endl;
        m3 = m2 ;
        //std::cout<<"same?";
        return ;
    }
    
    m3.resize( m1.size() * m2.size() );
    for(int i = 0 ; i < m3.size(); i++)
    m3[i].resize(m1[0].size() * m2[0].size());
    
    for(int i = 0 ; i < m1.size() ; i++)  
        for(int j = 0 ; j < m1[i].size(); j++) 
            for(int k = 0 ; k < m2.size() ; k++)
                for(int l = 0 ; l < m2[k].size(); l++)
                    m3[i * m2.size() + k][j * m2[k].size() + l] = m1[i][j] * m2[k][l] ; // this logic is provided by gpt
}

void kroneckerBUFF(const std::vector<std::vector<Coeff>>& gate , std::vector<std::vector<Coeff>>& res, int n_qbits, int tbit )
{
    if(n_qbits == 1)
    {
        res = gate; 
        return ;
    }

    std::vector<std::vector<Coeff>> tmp ;

    for(int i = 0 ; 1 ; )
    {
        if(i == tbit)
        kroneckerSUB( tmp, gate ,res);
        else
        kroneckerSUB( tmp, GATES::identity, res);

        if(++i < n_qbits)
        tmp = res;
        else 
        break ;
    }
    
    return ;
}

std::vector<double> generateUnitVector(int n) 
{
    std::vector<double> result;
    result.push_back(5);
    
    // Calculate the unit vector elements
    double sumOfSquares = 0.0;
    for (int i = 1; i <= n; ++i) {
        double value = static_cast<double>(rand()) / RAND_MAX;  // Random value between 0 and 1
        result.push_back(value);
        sumOfSquares += value * value;
    }

    // Normalize the vector to make the sum of squares equal to 1
    double normalizationFactor = 1.0 / std::sqrt(sumOfSquares);
    for (int i = 1; i <= n; ++i) {
        result[i] *= normalizationFactor;
    }

    return result;
}

struct Circuit {

    std::vector< std::vector< std::vector<Coeff > > > gord ; // gate order.... do I really need this ?
    std::vector< std::vector< std::vector<uint32_t> > > qlines ; 
    std::vector< std::vector<double> > table;
    std::vector<state_vector> sv ;
    //std::vector< uint32_t > table ;
    int nQ,nC; // number of q/classical bits //! replace with uint8_t to get mindfuck results (white exp snippet)
    uint64_t creg ; // 64 maximum classical bit lines 
    int gorder; // pointer to gate order 
    
    Circuit()
    {
        std::ifstream ip("circ.txt");
        std::string si,sj ; 
        std::getline(ip, si);
        std::istringstream iss0(si);
        int nS ; 
        iss0 >> nQ;
        iss0 >> nC; 
        iss0 >> nS; 
        std::cout << nQ << " " << nC << " " << nS << std::endl; 
        int i = 0 , j =0;
        uint32_t utility ; 
        for(int z = 1 ; z<= nQ; z++)
        {
            std::getline(ip, si) ;
            j=0;
            std::cout<< si <<std::endl; 
            std::istringstream iss1(si);
            
            std::vector< std::vector< uint32_t > > otmp ; 
            while( iss1 >> sj )
            {
                std::vector<uint32_t> tmp ;
                std::istringstream iss2(sj);
                iss2 >> utility;
                iss2.ignore();
                
                tmp.push_back(utility);
                
                if(utility == 1 || utility == 3)
                {
                    iss2 >> utility;
                    tmp.push_back(utility);
                }
                else if(utility == 2)
                {
                    for(int x = 0 ; x < 4 ; x++)
                    {
                        iss2 >> utility;
                        iss2.ignore();
                        tmp.push_back(utility);
                    }
                    
                }
                otmp.push_back(tmp);
            }
            qlines.push_back(otmp);
        }
        
        while( std::getline(ip, si) )
        {
            std::istringstream iss3(si);
            std::string token;
            int first; 
            iss3 >> first ;
            if(first == 5) // random val; 
            {
                iss3.ignore();
                iss3 >> first ; // 2nd val -- no.of qbits for random state 
                first = 2*( 1 << first ) ; // for each state -> 2 needed ( a + ib )
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
                table.push_back({1,value1,value2,value3,value4});
            }
            
        }
        std::cout<<"table: \n";
        for(int r = 0 ; r < table.size(); r++)
        {
            for(int r2 = 0 ; r2 < table[r].size(); r2++)
            std::cout << table[r][r2] <<" ";
            std::cout<<std::endl;
        }
    }
    
    void printQ()
    {
        std::cout<<"\nnumber of qbits: " << nQ<< " number of classical bits: " << nC << " commands: "<< qlines[0].size()<< std::endl ; //<< "no. of table entries: " << nS ;
        for(int i =0 ; i< nQ; i++)
        {
            std::cout << "qbit: "<<i <<": "; 
            for(int j = 0 ; j < qlines[i].size(); ++j)
            {
                for(int k = 0 ; k < qlines[i][j].size(); ++k)
                std::cout <<qlines[i][j][k]<<" ";
                std::cout << " ; ";
            }
            std::cout<<std::endl ;
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

struct result 
{ 
    int n_qbits; 
    std::map<int ,int> m ; // TODO rename, restructure 
    struct state_vector state;
    
    void print_state() {
        state.print();
    }
    
    void get_counts()
    {
        std::cout<<"classical register readings for the simulation: "<<std::endl ; 
        for(auto i = m.begin() ; i != m.end(); i++)
        {
            std::string binaryString;
            for (int j = n_qbits - 1; j >= 0; --j) {
                int bit = (i->first >> j) & 1;
                binaryString += (bit == 0) ? '0' : '1';
            }
            std::cout << binaryString << ": "<< i->second << std::endl ;
        }
        std::cout<<std::endl;
    }
};

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
    
    void applyGate(const std::string& gate , int t_qbit,std::array<double, 3> arr) // overloading applyGate 
    {
        if(gate == "U" || gate == "u")
        {
            int size = euler_container.size() ;
            euler_container.push_back({arr[0],arr[1],arr[2]});
            ORDER.push_back({1,t_qbit,5,size}); // send index
        }
    }
    
    struct result simulate(int shots = 1)
    {
        std::map<int,int> m ;
        while(shots--)
        {
            for(int i = 0; i < ORDER.size() ; i++)
            {
                if( ORDER[i][0] == 1 ) // applygate 
                {
                    std::vector<std::vector<Coeff>> matrix, buffedmatrix ; 
                    if(ORDER[i][2] == 1 )
                        matrix = GATES::hadamard ; 
                    else if(ORDER[i][2] == 2 )
                        matrix = GATES::pauli_X ;
                    else if(ORDER[i][2] == 3 )
                        matrix = GATES::pauli_Y ;
                    else if(ORDER[i][2] == 4 )
                        matrix = GATES::pauli_Z ;
                    else if(ORDER[i][2] == 5 )
                        matrix = GATES::unitaryGate( euler_container[ORDER[i][3]][0], euler_container[ORDER[i][3]][1], euler_container[ORDER[i][3]][2] );
                    
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
    
};


#endif
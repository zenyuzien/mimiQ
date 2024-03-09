#include <bits/stdc++.h>
#define root2 std::sqrt(2)

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
            std::cout<<"coefficient of " << binaryString << " is " << coeffs[i].real << " + "<<coeffs[i].complex << " i"<< std::endl ;
            std::cout<< binaryString << ": " << coeffs[i].amplitude()*100 <<"%"<< std::endl ;
        }
    
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


void applyGate(const std::vector<std::vector<Coeff>>& matrix, struct state_vector& sv)
{
    int n = 1 << sv.n_qbits;
    std::vector<Coeff> result(n);
    for (uint32_t i = 0; i < n; ++i) 
        for (uint32_t j = 0; j < n; ++j) 
            result[i] = result[i] + ( matrix[i][j] * sv.coeffs[j]);
        
    // Update the result in the input vector
    sv.coeffs = result;
    return ;
}

void kroneckerSUB(std::vector<std::vector<Coeff>>& m1, std::vector<std::vector<Coeff>>& m2, std::vector<std::vector<Coeff>>& m3)
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

void kroneckerBUFF( std::vector<std::vector<Coeff>>& gate , std::vector<std::vector<Coeff>>& res, int n_qbits, int tbit )
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

    std::vector< std::vector< std::vector<Coeff > > > gord ; // gate order 
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
                table.push_back({1,value1,value2});
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

state_vector Circuit::mimiQ()
{
    sv.resize(nQ);
    for(int i = 0 ; i < nQ; i++)
    ketZero(sv[i]);
    //sv[2].print();
    
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
                        sv[nQ-1].coeffs[(r-1)/2].complex = table[qlines[i][stage][1]][r+1] ;
                        r+=2 ;
                    }
                }
            }
            
        }
        std::cout<<"\nafter stage "<< stage ;
        sv[nQ-1].print();
        ++stage;
    }
        
    return sv[nQ -1] ;
}

// Change the return type to void or int, depending on your requirements
int measure(state_vector &sv, int qbit) {
    

    if (qbit < 0 || qbit >= sv.n_qbits) {
        std::cerr << "Invalid qubit index." << std::endl;
        return -1;
    }

    // Calculate the probability of measuring |1⟩ in the specified qubit
    double probability = 0.0;

    for (size_t i = 0; i < sv.coeffs.size(); ++i) {
        if ((i >> qbit) & 1) {
            probability += sv.coeffs[i].amplitude();
        }
    }

    probability *= 100;
    int randomValue = std::rand() % 100;

    //std::cout << randomValue << " " << probability << std::endl;

    return (randomValue < probability) ? 1 : 0;
}

void prob_qbit(const state_vector& sv, int tbit)
{
    int cnd = 1 << sv.n_qbits ;
    tbit=  1 << (sv.n_qbits-1-tbit);
    double prob = 0.0 ; 
    for(int i = 0 ; i < cnd ; i++)
    {
        if( i & tbit )
        {
            prob += sv.coeffs[i].amplitude();
        }
    }
    std::cout << "Probability of measuring |1> on qubit " << tbit << ": " << prob << std::endl;
    std::cout << "Probability of measuring |0> on qubit " << tbit << ": " << 1-prob << std::endl;
    return ;
}


int main()
{
    std::cout << "mimiq.h" << std::endl;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    
    Circuit qc ;
    qc.printQ();
    qc.mimiQ() ;
    qc.sv[qc.nQ -1].print();
    
    return 0 ;
}

/*
    struct state_vector q1,q2 ;
    
    ketZero(q1);
    ketZero(q2);
    applyGate(GATES::hadamard, q1);
    
    state_vector q12; 
    tensorProd( q1, q2, q12 );
    
    //printState(q12);
    GATES::controlled_pauli_X(q12,0,1);
    
    state_vector q0 ; 
    q0.n_qbits = 1;
    //a = 0.6413, b = 0.7542, c = 0.0745, d = 0.3930
    q0.coeffs.push_back({0.6413,0.7542});
    q0.coeffs.push_back({0.0745,0.3930});
    std::cout<<"intial q0 state: \n";
    printState(q0);
    state_vector q012; 
    tensorProd(q0,q12,q012);

    GATES::controlled_pauli_X(q012, 0,1);
    //printState(q012);
    
    std::vector<std::vector<Coeff>> h,z,x ; 
    kroneckerBUFF(GATES::hadamard, h, 3, 0 );
    kroneckerBUFF(GATES::pauli_Z, z, 3, 2 );
    kroneckerBUFF(GATES::pauli_X, x, 3, 2 );
    
    //printMatrix(h);
    
    applyGate(h,q012);
    std::cout<<"final q012 state: \n";
    printState(q012);
    
    //apply correction
    
    auto save = q012 ; 
    int shots = 1000 , randv; 
    
    std::map<int,int> map ;
    
    for( int i = 0 ; i < shots ; i++ )
    {
        q012 = save; 
        int fin = 0 ;
        
        auto prob_q0 = prob(q012, 0);
        
        int randv = std::rand() % 100;
        int q0BIT = (randv < prob_q0) ? 1 : 0;
        
        auto prob_q1 = prob(q012, 1);
        randv = std::rand() % 100;
        int q1BIT = (randv < prob_q1) ? 1 : 0;
        
        if(q0BIT)
        applyGate( z, q012 );
        if(q1BIT)
        applyGate( x, q012 );
        
        auto prob_q2 = prob(q012, 2);
        randv = std::rand() % 100;
        int q2BIT = (randv < prob_q2) ? 1 : 0;
        
        fin  |= q0BIT ;
        fin  = fin << 1 ;
        fin  |= q1BIT ;
        fin = fin << 1 ;
        fin |= q2BIT ;
        
        map[fin]++;
    }
    
    for(auto i = map.begin() ; i != map.end(); i++)
    std::cout << i->first << " : " << i->second << std::endl ;
    */

 //std::vector<double> probabilities ;
    //distribute(probabilities, q012);
    
    //int shots = 100 , measurement ; 
    /*
    for(int i = 0 ; i < shots ; i++)
    {

        measurement = measure( probabilities ) ;
        map[measurement]++;
    }
    
    measurement = 0; 
    std::cout<<"measurement = " << measurement << std::endl;
    
    
    if( measurement & 4)
    applyGate(z,q012); 
    if( measurement & 2)
    applyGate(x,q012);
    
    
    printState(q012);
    */
   // print_map_binary(map);

    
/*
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

void print_map_binary(const std::map<int, int>& mymap) {
  if (mymap.empty()) {
    std::cerr << "The map is empty." << std::endl;
    return;
  }
  std::cout<<std::endl;

  // Determine the number of bits required to represent the largest key
  int max_key = mymap.rbegin()->first;
  int bits_needed = 0;
  while (max_key > 0) {
    max_key >>= 1;
    bits_needed++;
  }

  // Iterate through the map and print key-value pairs
  for (const auto& pair : mymap) {
    // Specify the number of bits for the bitset (using a fixed size)
    std::cout << std::bitset<32>(pair.first).to_string().substr(32 - bits_needed)
              << " " << pair.second << std::endl;
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


/*
int main()
{
    std::cout << "mimiq.h" << std::endl;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    
    
    /*struct state_vector q0,q1,q2 ,q01, q012 ; 
    ketZero(q0);
    ketONE(q1);
    ketZero(q2);
    q01 = q0*q1;
    q012 = q01* q2; 
    q012.print();* /
    
    //Circuit qc ;
    //qc.printQ();
    //qc.mimiQ() ;
    
    
    std::unordered_map<int,int> m;
    
    for(int shots = 0 ; shots < 1 ; shots++)
    {
        state_vector q0(1),q1,q2,q01,q012;
    
        double a = 0.984, b = 0, c = 0.173, d = 0 ;
        q0.coeffs[0].real = a ;
        q0.coeffs[0].complex = b ; 
        q0.coeffs[1].real = c ;
        q0.coeffs[1].complex = d ; 
        //q0.print();
        
       
        ketZero(q1);
        ketZero(q2);
        
        q01 = q0*q1 ; 
        q012 = q01*q2 ; 
        
        //prob_qbit(q012,0);
        
       // q012.print();
        
        std::vector<std::vector<Coeff>> h ; // alternative to controlled bits ,z,x; 
        kroneckerBUFF(GATES::hadamard, h, 3, 1 ); // on 2nd qbit 
        applyGate(h,q012);
        GATES::controlled_pauli_X(q012,1,2);
        GATES::controlled_pauli_X(q012,0,1);
        
        kroneckerBUFF(GATES::hadamard, h, 3, 0 ); // on 1st qbit 
        
        //kroneckerBUFF(GATES::pauli_X, z, 3,2  ); // on 3rd qbit  alternative to controlled bits 
        //kroneckerBUFF(GATES::pauli_Z, x, 3,2 ); // on 3rd qbit  alternative to controlled bits 
        
        applyGate(h,q012);
        int C1 = measure(q012,0);
        int C2 = measure(q012,1);
         
        /*  alternative to controlled gates
        if(C2)
        applyGate(x,q012);
        if(C1)
        applyGate(z,q012);
        * /
        prob_qbit(q012,0);
        prob_qbit(q012,1);
        
        GATES::controlled_pauli_X(q012,1,2);
        GATES::controlled_pauli_Z(q012,0,2);
        
        int C3 = measure(q012,2);
        q012.print();
        
        prob_qbit(q012,2);
        
        int res = ( C1 << 2 ) + (C2 << 1 ) + C3 ;
        m[res]++;
    }
    
    for(auto i = m.begin() ; i != m.end(); i++)
    std::cout << i->first << ": "<< i->second << std::endl ;

    return 0 ;
}
*/




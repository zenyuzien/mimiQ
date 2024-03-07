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

namespace GATES {

 //the reason for vectors is to generalise 


 // NEED TO CONST !!!!!!!!!!!!!!

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
    
    void CNOT( struct state_vector& sv, int cbit, int tbit )
    {
        int cnd = 1 << sv.n_qbits , cmask = 1 << (sv.n_qbits-1-cbit) , tmask = 1 << (sv.n_qbits-1-tbit) , other ; 
        //std::cout<< cnd << " " << cmask << " " << tmask << std::endl; 
        std::vector<bool> visit( cnd , false);
        double tmp1, tmp2 ;
        
        // [ 00 01 10 11 ]
        
        for(int i = 0 ; i< cnd ; i++)
        {
           // std::cout<< i << std::endl;
            if( (i & cmask) && (visit[ i ]== false) )
            {
                if(i & tmask ) 
                other = i & ~(tmask);
                else 
                other  = i | tmask ; 
                //std::cout << "swap "<< i <<" with "<< other << std::endl ;

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
    //auto tmp = matrix[0][0]*sv.coeffs[0] + matrix[0][1]*sv.coeffs[1] ;
    //sv.coeffs[1] = matrix[1][0]*sv.coeffs[0] + matrix[1][1]*sv.coeffs[1] ;
    //sv.coeffs[0] = tmp ;
    
    int n = 1 << sv.n_qbits;
    std::vector<Coeff> result(n);
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < n; ++j) {
            result[i] = result[i] + ( matrix[i][j] * sv.coeffs[j]);
        }
    }

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
    {
        for(int j = 0 ; j < m1[i].size(); j++) 
        {
            for(int k = 0 ; k < m2.size() ; k++)
            {
                for(int l = 0 ; l < m2[k].size(); l++)
                {
                    m3[i * m2.size() + k][j * m2[k].size() + l] = m1[i][j] * m2[k][l] ; // this logic is provided by gpt
                }
            }
            
        }
    }
}

void kroneckerBUFF( std::vector<std::vector<Coeff>>& gate , std::vector<std::vector<Coeff>>& res, int n_qbits, int tbit )
{
   // int dim = 1 << n_qbits ;
   // res.resize(dim, std::vector<Coeff>(dim));

    if(n_qbits == 1)
    {
        res = gate; 
        return ;
    }

    std::vector<std::vector<Coeff>> tmp ;
    //  printMatrix(tmp);

    for(int i = 0 ; 1 ; )
    {
        //std::cout << i << std::endl;
        if(i == tbit)
        kroneckerSUB( tmp, gate ,res);
        else
        kroneckerSUB( tmp, GATES::identity, res);

        // printMatrix(res);

        if(++i < n_qbits)
        tmp = res;
        else 
        break ;
    }
    
    // std::cout<<"done";
    return ;
}

struct Circuit {

    std::vector< std::vector< std::vector<Coeff > > > gord ; // gate order 
    std::vector< std::vector< std::vector<uint32_t> > > qlines ; 
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
            ;// later.
        }
    }
    
    void printQ()
    {
        std::cout<<"\nnumber of qbits: " << nQ<< " number of classical bits: " << nC << std::endl ; //<< "no. of table entries: " << nS ;
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
    std::vector<state_vector> sv(nQ) ;
    for(int i = 0 ; i < nQ; i++)
    ketZero(sv[i]);
    //sv[2].print();
    
    if(nQ > 1 )
    for(int i = 0 ; i < nQ-1; i++) 
    sv[i+1] = sv[i]*sv[i+1];
        
    
        
    return sv[nQ -1] ;
}


int main()
{
    std::cout << "mimiq.h" << std::endl;
    
    /*
    struct state_vector q0,q1,q2 ,q01, q012 ; 
    ketZero(q0);
    ketZero(q1);
    ketZero(q2);
    tensorProd(q0,q1,q01);
    tensorProd(q012,q2,q012);
    */
    
    Circuit qc ;
    qc.printQ();
    qc.mimiQ() ;
    
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
    GATES::CNOT(q12,0,1);
    
    state_vector q0 ; 
    q0.n_qbits = 1;
    //a = 0.6413, b = 0.7542, c = 0.0745, d = 0.3930
    q0.coeffs.push_back({0.6413,0.7542});
    q0.coeffs.push_back({0.0745,0.3930});
    std::cout<<"intial q0 state: \n";
    printState(q0);
    state_vector q012; 
    tensorProd(q0,q12,q012);

    GATES::CNOT(q012, 0,1);
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
111 [7][0] [7][1] [7][2] [7][3] [7][4] [7][5] [7][6] [7][7] */




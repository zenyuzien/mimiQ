#ifndef _util 
#define _util 

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

std::vector<std::vector<Coeff>> inverse2x2(const std::vector<std::vector<Coeff>>& matrix)
{
    std::vector<std::vector<Coeff>> res; 
    res.resize(2, std::vector<Coeff>(2));
    Coeff num,det = ( (matrix[0][0]*matrix[1][1])-(matrix[0][1]*matrix[1][0]) ) ;
    num.real = det.real ;
    num.complex = -det.complex ;
    auto denom = (det.real*det.real) - (det.complex*det.complex);
    num.real /= denom ;
    num.complex /= denom ;
    
    //1 / a + i b  
    //a - i b /  a^2 - b^2 
    
    res[0][0] = num*matrix[1][1] ;
    
    res[0][1].real = - matrix[0][1].real ;
    res[0][1].complex = - matrix[0][1].complex ;
    res[0][1] = num*res[0][1] ; 
    
    res[1][0].real = - matrix[1][0].real ;
    res[1][0].complex = - matrix[1][0].complex ;
    res[1][0] = num*res[1][0] ; 
    
    res[1][1] = num* matrix[0][0] ;
    
    return res ;
}


#endif 
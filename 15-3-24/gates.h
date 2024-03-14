#ifndef _gates
#define _gates

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
                {-cos(psi)*sin(theta/2), -sin(psi)*sin(theta/2)}
            },
            {
                {cos(psi)*sin(theta/2), sin(psi)*sin(theta/2)},
                {cos(theta/2)*cos(lam + psi), cos(theta/2)*sin(lam + psi)}
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

#endif
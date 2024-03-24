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

#endif 
#ifndef _gates
#define _gates

#include "vector"
#include "../structures/structures.h"
// namespace for gates and their related functions
// every gate is a vector of vector of struct Coeff -- so a 3d array basically, the 3rd dimension for complex number
namespace GATES 
{
    std::vector < std::vector < Coeff >> unitaryGate(double theta, double psi, double lam);
    void controlled_pauli_X(struct state_vector & sv, int cbit, int tbit);
    void toffoli(struct state_vector & sv, int cbit1, int cbit2, int tbit);
    void controlled_hadamard(struct state_vector & sv, int cbit, int tbit);
    void matrixmultiply2x2(const std::vector<std::vector<Coeff>>& g1, const std::vector<std::vector<Coeff>>& g2, std::vector<std::vector<Coeff>>& g3) ;
    void controlled_pauli_Z(struct state_vector & sv, int cbit, int tbit) ;
    void printMatrix(const std::vector < std::vector < Coeff >> & matrix) ;
    void kroneckerSUB(const std::vector < std::vector < Coeff >> & m1,
    const std::vector < std::vector < Coeff >> & m2, std::vector < std::vector < Coeff >> & m3);
    void kroneckerBUFF(const std::vector < std::vector < Coeff >> & gate, std::vector < std::vector < Coeff >> & res, int n_qbits, int tbit) ;
    std::vector < std::vector < Coeff >> inverse2x2(const std::vector < std::vector < Coeff >> & matrix);
}


#endif


#include "gates.h"
#include "../structures/structures.h"
#include <array>
#include <sstream>
#include <vector>

    // all basisc gates like unitary gate, hadamard gate, pauli's x,y,z are all derieved from this unitary gate which takes euler angles as input
    // the mathematical formula for this matrix is based on qiskit's implementation
    std::vector < std::vector < Coeff >> GATES::unitaryGate(double theta, double psi, double lam) {
        return {
            {
                { cos(theta / 2), 0 }, {-cos(lam) * sin(theta / 2), -sin(lam) * sin(theta / 2) }
            },
            {
                { cos(psi) * sin(theta / 2), sin(psi) * sin(theta / 2) },
                { cos(theta / 2) * cos(lam + psi), cos(theta / 2) * sin(lam + psi) }
            }
        };
    }

    // while applying a nomral independant gate can be generalized, controlled gates are not. Is that even possible to do os ?
    // controlled pauli X gate - when control qbit is 1 on {|0>, |1>} basis, then target qbit is flipped.
    void GATES::controlled_pauli_X(struct state_vector & sv, int cbit, int tbit) {

        // the cnd is just the number of states possible from n_qbits, which is 2^n_qbits
        // the mask is created from the specified positions of control bit and target bit
        // logical operation with the mask will result in 1 or 0
        int cnd = 1 << sv.n_qbits, cmask = 1 << (sv.n_qbits - 1 - cbit), tmask = 1 << (sv.n_qbits - 1 - tbit), other;
        std::vector < bool > visit(cnd, false);
        // the reason why I have a visit array is because I want to swap only once for a given set 1xxt 
        // where t(target bit) is 0 or 1 and the most signioficant 1 is assumed to be the control bit
        // so if I encounter either of 1xx0 or 1xx1, I will swap the Coeffs and then I will set flag of both states as true to prevent re-swapping which will be pointless
        double tmp1, tmp2;
        for(int i = 0; i < cnd; i++) {
            if((i & cmask) && (visit[i] == false)) { // when unvisited and mask condition satisfies -- meaning the position is right
                if(i & tmask) // if tbit 1 -> other = xxxx & 11011 => 0 at tbit 
                    other = i & ~(tmask); // retrieving the basis state of the complement of the current 1xxt
                else other = i | tmask; // else, other = xxxx | 00100 => 1 at tbit 

                // swapping 
                tmp1 = sv.coeffs[i].real;
                tmp2 = sv.coeffs[i].complex;
                sv.coeffs[i].real = sv.coeffs[other].real;
                sv.coeffs[i].complex = sv.coeffs[other].complex;
                sv.coeffs[other].real = tmp1;
                sv.coeffs[other].complex = tmp2;
                visit[i] = visit[other] = true;
            }
        }
    }
    void GATES::toffoli(struct state_vector & sv, int cbit1, int cbit2, int tbit) {

        // the cnd is just the number of states possible from n_qbits, which is 2^n_qbits
        // the mask is created from the specified positions of control bit and target bit
        // logical operation with the mask will result in 1 or 0
        int cnd = 1 << sv.n_qbits, 
        cmask1 = 1 << (sv.n_qbits - 1 - cbit1),
        cmask2 = 1 << (sv.n_qbits - 1 - cbit2),
        tmask = 1 << (sv.n_qbits - 1 - tbit), other;
        std::vector < bool > visit(cnd, false);
        // the reason why I have a visit array is because I want to swap only once for a given set 1xxt 
        // where t(target bit) is 0 or 1 and the most signioficant 1 is assumed to be the control bit
        // so if I encounter either of 1xx0 or 1xx1, I will swap the Coeffs and then I will set flag of both states as true to prevent re-swapping which will be pointless
        double tmp1, tmp2;
        for(int i = 0; i < cnd; i++) {
            if(((i & cmask1)&&(i & cmask2)) && (visit[i] == false)) { // when unvisited and mask condition satisfies -- meaning the position is right
                if(i & tmask) // if tbit 1 -> other = xxxx & 11011 => 0 at tbit 
                    other = i & ~(tmask); // retrieving the basis state of the complement of the current 1xxt
                else other = i | tmask; // else, other = xxxx | 00100 => 1 at tbit 

                // swapping 
                tmp1 = sv.coeffs[i].real;
                tmp2 = sv.coeffs[i].complex;
                sv.coeffs[i].real = sv.coeffs[other].real;
                sv.coeffs[i].complex = sv.coeffs[other].complex;
                sv.coeffs[other].real = tmp1;
                sv.coeffs[other].complex = tmp2;
                visit[i] = visit[other] = true;
            }
        }
    }

    // controlled hadamard-  qbit in control and target also qbit, along with state vector of course
    void GATES::controlled_hadamard(struct state_vector & sv, int cbit, int tbit) {
        // masking as explaine in controlled pauli x
        int cnd = 1 << sv.n_qbits,
            cmask = 1 << (sv.n_qbits - 1 - cbit),
            tmask = 1 << (sv.n_qbits - 1 - tbit),
            other;
        // I am not directly updating the passed state vector because I need to use all the states passed before updating any, so I just create a new vector and assign it in the end
        struct state_vector tmp(sv.n_qbits);
        for(int i = 0; i < cnd; i++) {
            if(i & cmask) {
                //std::cout << "in "<< i << std::endl ;
                // in high level, we apply hadamard gate to a qbit. A single qbit, hence whatever is the Coeff value, eg. C|1> or C|0>
                // it will be 1/root2 * C|0> +/- 1/root2* C|1>
                //so in the new vector, we add the 2 parts in their respective places
                struct Coeff coeff1(1 / root2, 0); // 1/root2 
                if(i & tmask) // if target bit is 1 here, x|1>x --> 1/root2 *[ x|0>x - x|1>x ]
                {
                    other = i & ~(tmask);
                    tmp.coeffs[i] = tmp.coeffs[i] - (coeff1) * sv.coeffs[i];
                }
                else // if target bit is 0 here, x|0>x --> 1/root2 *[ x|0>x + x|1>x ]
                {
                    other = i | tmask;
                    tmp.coeffs[i] = tmp.coeffs[i] + (coeff1) * sv.coeffs[i];
                }
                tmp.coeffs[other] = tmp.coeffs[other] + (coeff1) * sv.coeffs[i];
            }
            // if not involved in control bit, just assign whats existing as there wont be any change required
            else tmp.coeffs[i] = sv.coeffs[i];
        }
        // assign to make changes
        sv = tmp;
        return;
    }

    // simple 2x2 matrix multiplication
    void GATES::matrixmultiply2x2(const std::vector<std::vector<Coeff>>& g1, const std::vector<std::vector<Coeff>>& g2, std::vector<std::vector<Coeff>>& g3) 
    {
        g3[0][0] = (g1[0][0]*g2[0][0]) + (g1[0][1]*g2[1][0]) ;
        g3[0][1] = (g1[0][0]*g2[0][1]) + (g1[0][1]*g2[1][1]) ;
        g3[1][0] = (g1[1][0]*g2[0][0]) + (g1[1][1]*g2[1][0]) ;
        g3[1][1] = (g1[1][0]*g2[0][1]) + (g1[1][1]*g2[1][1]) ;
    }

    // when control bit is 1, then target bit coeff should be multiplied with (-1 + 0i)
    void GATES::controlled_pauli_Z(struct state_vector & sv, int cbit, int tbit) {
        int cnd = 1 << sv.n_qbits, cmask = 1 << (sv.n_qbits - 1 - cbit), tmask = 1 << (sv.n_qbits - 1 - tbit);
        struct Coeff tmp(-1, 0); 
        for(int i = 0; i < cnd; i++)
            if((i & cmask) && (i & tmask)) // if cbit and tbit 1 flip sign.
                sv.coeffs[i] = sv.coeffs[i] * tmp;
    }

    // to print matrix
    void GATES::printMatrix(const std::vector < std::vector < Coeff >> & matrix) {
        if(matrix.size() == 0) {
            std::cout << "empty " << std::endl;
            return;
        }
        for(const auto & row: matrix) {
            for(const auto & element: row) {
                std::cout << "{" << element.real << ", " << element.complex << "} ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    // kroneckerSUB is the utility of KnroneckerBUFF, it will create a cross product (or kronecker product ) of 2 gates and assign it to m3
    void GATES::kroneckerSUB(const std::vector<std::vector<Coeff>> &m1,
                            const std::vector<std::vector<Coeff>> &m2,
                            std::vector<std::vector<Coeff>> &m3) {
        // Handle edge case
        if (m1.empty()) {
            m3 = m2;
            return;
        }

        // Precompute dimensions
        size_t m1_rows = m1.size();
        size_t m1_cols = m1[0].size();
        size_t m2_rows = m2.size();
        size_t m2_cols = m2[0].size();

        size_t m3_rows = m1_rows * m2_rows;
        size_t m3_cols = m1_cols * m2_cols;

        // Allocate memory for m3
        m3.resize(m3_rows);
        for (auto &row : m3) {
            row.resize(m3_cols);
        }

        // Compute Kronecker product
        for (size_t i = 0; i < m1_rows; ++i) {
            for (size_t j = 0; j < m1_cols; ++j) {
                for (size_t k = 0; k < m2_rows; ++k) {
                    for (size_t l = 0; l < m2_cols; ++l) {
                        m3[i * m2_rows + k][j * m2_cols + l] = m1[i][j] * m2[k][l];
                    }
                }
            }
        }
    }


    // when a gate has to be applied to a qbit.. but we have a state vector of multiple states of multiple qbits, 
    // we need to build a matrix which when multiplied with the state vector, the areas of state_vector which correspond to the atrget qbit only has to be changed.
    // the kroneckerbuff will create the matrix.
    // for a statevector which has 5 qbits like xxxxx, if we want to apply a unitary gate U to 2nd qbit, 
    // then our matrix would like I x U x I x I x I which is multiplying identity gates using cross product in appropriate order
    void GATES::kroneckerBUFF(const std::vector < std::vector < Coeff >> & gate, std::vector < std::vector < Coeff >> & res, int n_qbits, int tbit) {
        if(n_qbits == 1) {
        // if only 1 qbit, then no need of buffing so just directly return gate 
            res = gate;
            return;
        }
        std::vector < std::vector < Coeff >> tmp; // tmp will be buffed from empty, and will be assigned to res on completion
        for(int i = 0; 1;) { 
            // loop for n_qbit times, cross product with Identity every iteration except for the position of target qbit
            if(i == tbit) kroneckerSUB(tmp, gate, res);
            else kroneckerSUB(tmp, GATES::unitaryGate(0, 0, 0), res);
            if(++i < n_qbits) tmp = res;
            else break;
        }
        return;
    }

    // given a 2x2 mateix of complex numbers (or type Coeff), this function will return the inverse of it
    std::vector < std::vector < Coeff >> GATES::inverse2x2(const std::vector < std::vector < Coeff >> & matrix) {
        std::vector < std::vector < Coeff >> res;
        res.resize(2, std::vector < Coeff > (2));
        Coeff num, det = ((matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]));
        num.real = det.real;
        num.complex = -det.complex;
        auto denom = (det.real * det.real) + (det.complex * det.complex);
        num.real /= denom;
        num.complex /= denom;
        //std::cout<<"1/det: "<<num.real<<std::endl;
        //1 / a + i b  
        //a - i b /  a^2 - (ib)^2 
        res[0][0] = matrix[1][1];
        res[0][1].real = -matrix[0][1].real;
        res[0][1].complex = -matrix[0][1].complex;
        res[1][0].real = -matrix[1][0].real;
        res[1][0].complex = -matrix[1][0].complex;
        res[1][1].real = matrix[0][0].real;
        res[1][1].complex = (matrix[0][0].complex);
        res[0][0] = num * res[0][0];
        res[0][1] = num * res[0][1];
        res[1][0] = num * res[1][0];
        res[1][1] = num * res[1][1];
        return res;
    }

    
    void Qcircuit::applyUtility(const std::string &gate, int& G, int& esize ,std::vector<double> arr = std::vector<double>(),std::vector<std::vector<Coeff> >& matrix = empty_matrix)
    {
        if (gate == "H" || gate == "h")
        {
            G = 1;
            matrix = GATES::unitaryGate (PI / 2, 0, PI);
        }
            
        else if (gate == "x" || gate == "X")
        {
            G = 2;
            matrix = GATES::unitaryGate (PI, 0, PI);
        }
        else if (gate == "y" || gate == "Y")
        {
            G = 3;
            matrix = GATES::unitaryGate (PI, PI / 2, PI / 2);
        }
         else if (gate == "z" || gate == "Z")
        {
            G = 4;
            matrix = GATES::unitaryGate (0, 0, PI);
        }   
        else if (gate == "s" || gate == "S")
        {
            G = 5;
            matrix = GATES::unitaryGate (0, 0, PI/2);
        }      
        else if (gate == "sd" || gate == "Sd")
        {
            G = 6;
            matrix = GATES::unitaryGate (0, 0, -PI/2);
        }   
        else if (gate == "T" || gate == "t")
        {
            G = 7;
            matrix = GATES::unitaryGate (0, 0, PI/4);
        }
        else if (gate == "Td" || gate == "td")
        {
            G = 8;
            matrix = GATES::unitaryGate (0, 0, -PI/4);
        }   


        else if (gate == "rx" || gate == "Rx")
        {
            G = 105;
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], 0, 0 });
                
                matrix = GATES::unitaryGate (
                                        euler_container[esize][0], -PI/2, PI/2);
        }     
        else if (gate == "ry" || gate == "Ry")
        {
            G = 106;
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], 0, 0 });
                
                matrix = GATES::unitaryGate (
                                        euler_container[esize][0], 0, 0);
        }   
        else if (gate == "U3" || gate == "u3" || gate == "u" || gate == "U")
            {
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], arr[2] });
                G= 101;
                matrix = GATES::unitaryGate (
                                        euler_container[esize][0],
                                        euler_container[esize][1],
                                        euler_container[esize][2]);
            }
        else if (gate == "IU3" || gate == "iu3" || gate == "IU"
                 || gate == "iu") // todo: generalise u1 u2 u3 inverse
            {
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], arr[2] });
                G = 102;
                matrix = GATES::inverse2x2 (
                                        GATES::unitaryGate (
                                            euler_container[esize][0],
                                            euler_container[esize][1],
                                            euler_container[esize][2]));
            }
        else if (gate == "U1" || gate == "u1")
            {
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], 0, 0 });
                G=103;
                matrix = GATES::unitaryGate (
                                        0, 0, euler_container[esize][0]);
            }
        else if (gate == "U2" || gate == "u2")
            {
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], 0 });
                G =104;
                matrix = GATES::unitaryGate (
                                        PI/2, euler_container[esize][0],
                                        euler_container[esize][1]);
            }
    }

    void
    Qcircuit::applyGate (const std::string &gate, int t_qbit,
               std::vector<double> arr, bool ORDERinc )
    {
        std::vector<std::vector<Coeff> > matrix,buffedmatrix;
        int esize = -1,G;
   
        applyUtility(gate,G,esize,arr,matrix);
            
            //std::cout<<"p: "<<1 <<" "<<t_qbit << " " << G << " " << esize << std::endl;
            if(ORDERinc)
            ORDER.push_back ({ 1, t_qbit, G, esize }); // send index
                            
                            GATES::kroneckerBUFF (matrix, buffedmatrix, n_qbits,
                                                  t_qbit);
                            int n = 1 << n_qbits;
                            std::vector<Coeff> result (n, { 0, 0 });
                            for (int ii = 0; ii < n; ++ii)
                                {
                                    for (int j = 0; j < n; ++j)
                                        result[ii] = result[ii]
                                                     + (buffedmatrix[ii][j]
                                                        * state.coeffs[j]);
                                }
                            state.coeffs = result;
    }
      void
    Qcircuit::applyControlledGate (const std::string &gate, const std::string c_qbit_str,
                         int t_qbit,  std::vector<double> arr)
    {
        std::istringstream iss (c_qbit_str);
        char cq1,cq2,d='0';
        int control_bit1,control_bit2,c_or_q1,c_or_q2,G;
        iss >> cq1 >> control_bit1 >> d >> cq2 >> control_bit2 ;

        if(d=='0')
        {
            c_or_q1 = (cq1 == 'q') ? 1 : 0 ;

            if(c_or_q1 ==0)
            {
                                        if(c_reg & (1 << (n_cbits - 1 - control_bit1)))
                                                applyGate(gate,t_qbit,arr,false);
            }

            int esize;
            applyUtility(gate,G,esize,arr);

            ORDER.push_back ({ 2, t_qbit, G, c_or_q1, control_bit1, esize});
            
            if(c_or_q1 == 1)
            { // CONTROL ON ZERO WITH Q needs to be done. TODO
                                        if (G == 1) // gate
                                            GATES::controlled_hadamard (
                                                state, control_bit1,
                                                t_qbit); // qbit, tqbit

                                        else if (G== 2)
                                            GATES::controlled_pauli_X (
                                                state, control_bit1, t_qbit);

                                        else if (G
                                                == 4) // 4 - pauli z , 3 - pauli y
                                            GATES::controlled_pauli_Z (
                                                state, control_bit1, t_qbit);
            }
        }
        else if(d==',')
        {
            c_or_q1 = (cq1 == 'q') ? 1 : 0;
            c_or_q2 = (cq2 == 'q') ? 1 : 0;
            if((c_or_q1 == c_or_q2) &&(c_or_q1==0))
            {
                if((c_reg & (1 << (n_cbits - 1 - control_bit1)))&&(c_reg & (1 << (n_cbits - 1 - control_bit2))))
                    applyGate(gate,t_qbit,arr,false);
            }
            int esize;
            applyUtility(gate,G,esize,arr);
            // circuit
            ORDER.push_back ({ 4, t_qbit, G, c_or_q1, control_bit1, control_bit2, esize});
            
            if((c_or_q1 == c_or_q2) &&(c_or_q1==1))
            { // CONTROL ON ZERO WITH Q needs to be done. TODO
                                        if (G== 2)
                                            GATES::toffoli (
                                                state, control_bit1,control_bit2, t_qbit);
            }
        }


        
    }
    
#ifndef _gates
#define _gates
namespace GATES {
    std::vector < std::vector < Coeff >> unitaryGate(double theta, double psi, double lam) {
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
    void controlled_pauli_X(struct state_vector & sv, int cbit, int tbit) {
        int cnd = 1 << sv.n_qbits, cmask = 1 << (sv.n_qbits - 1 - cbit), tmask = 1 << (sv.n_qbits - 1 - tbit), other;
        std::vector < bool > visit(cnd, false);
        double tmp1, tmp2;
        for(int i = 0; i < cnd; i++) {
            if((i & cmask) && (visit[i] == false)) {
                if(i & tmask) // if tbit 1 -> other = xxxx & 11011 => 0 at tbit 
                    other = i & ~(tmask);
                else other = i | tmask; // else, other = xxxx | 00100 => 1 at tbit 
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
    void controlled_hadamard(struct state_vector & sv, int cbit, int tbit) {
        int cnd = 1 << sv.n_qbits,
            cmask = 1 << (sv.n_qbits - 1 - cbit),
            tmask = 1 << (sv.n_qbits - 1 - tbit),
            other;
        struct state_vector tmp(sv.n_qbits);
        for(int i = 0; i < cnd; i++) {
            if(i & cmask) {
                //std::cout << "in "<< i << std::endl ;
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
            else tmp.coeffs[i] = sv.coeffs[i];
        }
        sv = tmp;
        return;
    }
    void controlled_pauli_Z(struct state_vector & sv, int cbit, int tbit) {
        int cnd = 1 << sv.n_qbits, cmask = 1 << (sv.n_qbits - 1 - cbit), tmask = 1 << (sv.n_qbits - 1 - tbit);
        struct Coeff tmp(-1, 0);
        for(int i = 0; i < cnd; i++)
            if((i & cmask) && (i & tmask)) // if cbit and tbit 1 flip sign.
                sv.coeffs[i] = sv.coeffs[i] * tmp;
    }
    void printMatrix(const std::vector < std::vector < Coeff >> & matrix) {
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
    void kroneckerSUB(const std::vector < std::vector < Coeff >> & m1,
        const std::vector < std::vector < Coeff >> & m2, std::vector < std::vector < Coeff >> & m3) {
        if(m1.size() == 0) {
            m3 = m2;
            return;
        }
        m3.resize(m1.size() * m2.size());
        for(long long unsigned int i = 0; i < m3.size(); i++) 
            m3[i].resize(m1[0].size() * m2[0].size());
        for(long long unsigned int i = 0; i < m1.size(); i++)
            for(long long unsigned int j = 0; j < m1[i].size(); j++)
                for(std::vector<std::vector<Coeff> >::size_type k = 0; k < m2.size(); k++)
                    for(long long unsigned int l = 0; l < m2[k].size(); l++) 
                        m3[i * m2.size() + k][j * m2[k].size() + l] = m1[i][j] * m2[k][l];
    }
    void kroneckerBUFF(const std::vector < std::vector < Coeff >> & gate, std::vector < std::vector < Coeff >> & res, int n_qbits, int tbit) {
        if(n_qbits == 1) {
            res = gate;
            return;
        }
        std::vector < std::vector < Coeff >> tmp;
        for(int i = 0; 1;) {
            if(i == tbit) kroneckerSUB(tmp, gate, res);
            else kroneckerSUB(tmp, GATES::unitaryGate(0, 0, 0), res);
            if(++i < n_qbits) tmp = res;
            else break;
        }
        return;
    }
    std::vector < std::vector < Coeff >> inverse2x2(const std::vector < std::vector < Coeff >> & matrix) {
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
}
#endif
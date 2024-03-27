#ifndef _structures
#define _structures
#define root2 std::sqrt(2)
#define PI acos(0.0) * 2

struct Coeff {
    double real;
    double complex;
    Coeff(): real(0), complex(0) {}
    Coeff(double r, double c): real(r), complex(c) {}

    Coeff operator * (const Coeff & other) const {
        Coeff result;
        result.real = real * other.real - complex * other.complex;
        result.complex = real * other.complex + complex * other.real;
        return result;
    }
    Coeff operator + (const Coeff & other) const {
        Coeff result;
        result.real = real + other.real;
        result.complex = complex + other.complex;
        return result;
    }
    Coeff operator - (const Coeff & other) const {
        Coeff result;
        result.real = real - other.real;
        result.complex = complex - other.complex;
        return result;
    }
    Coeff operator / (const Coeff & other) const {
        Coeff result;
        result.real = ((real * other.real) + (complex * other.complex)) / ((other.real * other.real) + (other.complex * other.complex));
        result.complex = ((other.complex * real) - (real * other.complex)) / ((other.real * other.real) + (other.complex * other.complex));
        return result;
    }
    double amplitude() const { //amp squared actually
        return (real * real) + (complex * complex);
    }
};
struct state_vector {
    int n_qbits; // 4 
    std::vector < Coeff > coeffs; // 16 * (1<< nq) // for 2 qbits -> 544 cbits, for 3qbits -> 1056 cbits
    state_vector(): n_qbits(0) {}
    state_vector(int n) {
        n_qbits = n;
        coeffs.resize(1 << n, { 0, 0 });
    }
    state_vector& operator=(const state_vector& other) {
        if (this != &other) { // Check for self-assignment
            n_qbits = other.n_qbits;
            coeffs = other.coeffs; 
        }
        return *this;
    }
    state_vector(const struct state_vector & sv) {
        n_qbits = sv.n_qbits;
        coeffs = sv.coeffs;
    }
    struct state_vector operator * (const struct state_vector & other) // returns a tensor product of 2 state vectors 
    {
        struct state_vector sv3; //(*this).n_qbits == n_qbits 
        sv3.n_qbits = n_qbits + other.n_qbits;
        int cnd1 = 1 << n_qbits, cnd2 = 1 << other.n_qbits;
        for(int i = 0; i < cnd1; i++)
            for(int j = 0; j < cnd2; j++) sv3.coeffs.push_back(coeffs[i] * other.coeffs[j]);
        return sv3;
    }
    void print() {
        wr<< "\\[\n\\begin{array}{@{}llll@{}}\n";
        std::cout << std::endl << "number of qbits = " << n_qbits << std::endl;
    //    std::cout << std::endl << "number of qbits = " << n_qbits << "\\\\";
        int cnd = 1 << n_qbits;
        for(int i = 0; i < cnd; i++) {
            std::string binaryString;
            for(int j = n_qbits - 1; j >= 0; --j) {
                int bit = (i >> j) & 1;
                binaryString += (bit == 0) ? '0' : '1';
            }
           wr << "\\text{" << binaryString << ":} & " << coeffs[i].amplitude() * 100 << "\\% & " << coeffs[i].real << " |0\\rangle &  " ;
           wr << coeffs[i].complex << " |1\\rangle \\\\" << std::endl;
           std::cout << binaryString << " =>" << coeffs[i].amplitude() * 100 << "% ( " << coeffs[i].real << " |0> + " << coeffs[i].complex << " i |1> )\n";
        }
        wr << "\\end{array}\n\\]\n";
        std::cout << std::endl;
/*
\[
\begin{array}{@{}llll@{}}
\text{000:} & 50\% & 0.707107|0\rangle & - 8.65927e-17 |1\rangle \\
\text{001:} & 4.95606e-31\% & -5.55112e-17|0\rangle & + 4.32964e-17 |1\rangle \\
\text{010:} & 1.87458e-31\% & 2.65105e-33|0\rangle & + 4.32964e-17 |1\rangle \\
\text{011:} & 1.65236e-63\% & -2.65105e-33|0\rangle & + 3.08149e-33 |1\rangle \\
\text{100:} & 3.85186e-32\% & 1.96262e-17|0\rangle & + 1.54074e-32 |1\rangle \\
\text{101:} & 2.07487e-30\% & 1.37383e-16|0\rangle & + 4.32964e-17 |1\rangle \\
\text{110:} & 1.42005e-30\% & 1.11022e-16|0\rangle & + 4.32964e-17 |1\rangle \\
\text{111:} & 50\% & 0.707107|0\rangle &  - 8.65927e-17 |1\rangle \\
\end{array}
\]

*/
    }
    void printprobs()
    {
        wr << "\\begin{align*}\n";
        int cnd = 1 << n_qbits ;
        for(int bit = 0; bit <  n_qbits ; bit++)
        {
            auto tmask = 1 << (n_qbits - 1 - bit);
            double prob = 0.0;
            for(int i = 0; i < cnd; i++)
                if(i & tmask) prob += coeffs[i].amplitude();
            std::cout << "for qbit " << bit << ": " <<"Prob of |1> : " << prob << ",  ";
            std::cout << "Prob of |0> : " << 1 - prob << std::endl;
            wr << "\\text{For qbit "<< bit <<"} \\quad & \\text{Probability of} |1\\rangle: " << prob << ", \\text{Probability of} |0\\rangle: "
            << 1 - prob << " \\\\\n";
        }
        std::cout<<std::endl;
        wr<<"\\end{align*}\n";
    }
/*
\begin{align*}
\text{For qubit 0:} \quad & \text{Probability of } |1\rangle: 0.5, \text{ Probability of } |0\rangle: 0.5 \\
\text{For qubit 1:} \quad & \text{Probability of } |1\rangle: 0.5, \text{ Probability of } |0\rangle: 0.5 \\
\end{align*}
*/
};
struct result {
    int n_cbits;
    std::map < int, int > m; // TODO rename, restructure 
    struct state_vector state;
    void print_state() {
        state.print();
    }
    void get_counts() {
        wr << "\\begin{align*}\n\\text{Classical register readings for the simulation:} \\\\\n\\begin{array}{@{}ll@{}}\n";
        std::cout << "classical register readings for the simulation: " << std::endl;
        for(auto i = m.begin(); i != m.end(); i++) {
            std::string binaryString;
            for(int j = n_cbits - 1; j >= 0; --j) {
                int bit = (i -> first >> j) & 1;
                binaryString += (bit == 0) ? '0' : '1';
            }
            std::cout << binaryString << ": " << i -> second << std::endl;  
            wr << binaryString << ": & "<< i -> second << "\\\\\n";
        }
        std::cout << std::endl;
        wr<< "\\end{array}\n\\end{align*}\n";
    }/*
    \begin{align*}
\text{Classical register readings for the simulation:} \\
\begin{array}{@{}ll@{}}
000: & 251 \\
001: & 4 \\
010: & 261 \\
011: & 4 \\
100: & 239 \\
101: & 5 \\
110: & 253 \\
111: & 7 \\
\end{array}
\end{align*}
*/
    void get_probs()
    {
        state.printprobs();
    }
};
std::vector < std::vector < Coeff >> ma1, ma2, ma3;
#endif
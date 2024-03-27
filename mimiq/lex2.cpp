#include <iostream>
#include <fstream>
#include <cstdlib> // for system()
#include <vector>

struct Coeff {
    double real;
    double complex;
    Coeff(): real(0), complex(0) {}
    Coeff(double r, double c): real(r), complex(c) {}

    Coeff operator =(const Coeff & other) {
        real=other.real;
        complex= other.complex;

    }

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
        std::cout << std::endl << "number of qbits = " << n_qbits << std::endl;
        int cnd = 1 << n_qbits;
        for(int i = 0; i < cnd; i++) {
            std::string binaryString;
            for(int j = n_qbits - 1; j >= 0; --j) {
                int bit = (i >> j) & 1;
                binaryString += (bit == 0) ? '0' : '1';
            }
            std::cout << binaryString << " =>" << coeffs[i].amplitude() * 100 << "% ( " << coeffs[i].real << " |0> + " << coeffs[i].complex << " i |1> )\n";
        }
        std::cout << std::endl;
    }
    void printprobs()
    {
        int cnd = 1 << n_qbits ;
        for(int bit = 0; bit <  n_qbits ; bit++)
        {
            auto tmask = 1 << (n_qbits - 1 - bit);
            double prob = 0.0;
            for(int i = 0; i < cnd; i++)
                if(i & tmask) prob += coeffs[i].amplitude();
            std::cout << "for qbit " << bit << ": " <<"Prob of |1> : " << prob << ",  ";
            std::cout << "Prob of |0> : " << 1 - prob << std::endl;
        }
        std::cout<<std::endl;
    }
};
void generateLatexProbs(const state_vector& sv, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return;
    }

    file << "\\documentclass{article}\n";
    file << "\\begin{document}\n\n";

    int cnd = 1 << sv.n_qbits;
    for (int bit = 0; bit < sv.n_qbits; bit++) {
        auto tmask = 1 << (sv.n_qbits - 1 - bit);
        double prob = 0.0;
        for (int i = 0; i < cnd; i++)
            if (i & tmask) prob += sv.coeffs[i].amplitude();

        // Generate LaTeX code for probabilities
        file << "For qbit " << bit << ": Probability of $|1\\rangle$: " << prob << ", Probability of $|0\\rangle$: " << (1 - prob) << "\n\n";
    }

    file << "\\end{document}\n";
    file.close();
}

void generatePdf(const std::string& texFilename, const std::string& pdfFilename) {
    std::string compileCommand = "pdflatex bot_circuit.tex > NUL 2>&1";
    int exitCode = std::system(compileCommand.c_str());

    if (exitCode == 0) {
        std::string moveCommand = "mv " + texFilename + ".pdf " + pdfFilename;
        std::system(moveCommand.c_str());
        std::cout << "PDF generated: " << pdfFilename << std::endl;
    } else {
        std::cerr << "Error occurred during compilation." << std::endl;
    }
}

int main() {
    // Create a state vector
    state_vector sv(2);
    sv.coeffs[0] = Coeff(0.5, 0.5);
    sv.coeffs[1] = Coeff(0.3, 0.7);
    sv.coeffs[2] = Coeff(0.1, 0.9);
    sv.coeffs[3] = Coeff(0.6, 0.4);

    // Generate LaTeX code for probabilities
    generateLatexProbs(sv, "probs.tex");

    // Compile LaTeX file into PDF
    generatePdf("probs.tex", "probs.pdf");

    return 0;
}

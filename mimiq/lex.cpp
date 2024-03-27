#include <bits/stdc++.h>

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
// Generate the LaTeX code for the quantum circuit

void generateLatexProbs(const state_vector& sv, std::ofstream& file) {

    if (!file.is_open()) {
        std::cerr << "Error: file open chesi saavu." << std::endl;
        return;
    }

    int cnd = 1 << sv.n_qbits;
    for (int bit = 0; bit < sv.n_qbits; bit++) {
        auto tmask = 1 << (sv.n_qbits - 1 - bit);
        double prob = 0.0;
        for (int i = 0; i < cnd; i++)
            if (i & tmask) prob += sv.coeffs[i].amplitude();

        // Generate LaTeX code for probabilities
        file << "For qbit " << bit << ": Probability of $|1\\rangle$: " << prob << ", Probability of $|0\\rangle$: " << (1 - prob) << "\n\n";
    }

}

int main() {


    std::ofstream of("trial.tex") ;

    //file << "\\documentclass{article}\n";
    of << "\\documentclass{article}\n";
    of << "\\usepackage{qcircuit}\n";
    of << "\\begin{document}\n";
    of << "Let ZENyuzien cook\n";

    of << "\\begin{displaymath}\n";
    of << "\\Qcircuit @C=1em @R=1em {\n";
    of << "& \\gate{H} & \\ctrl{1} & \\qw \\\n";
    of << "& \\qw & \\targ & \\qw \\ \n";
    of << "\\end{displaymath}\n";


    state_vector sv(2);
    sv.coeffs[0] = Coeff(0.5, 0.5);
    sv.coeffs[1] = Coeff(0.3, 0.7);
    sv.coeffs[2] = Coeff(0.1, 0.9);
    sv.coeffs[3] = Coeff(0.6, 0.4);

    // Generate LaTeX code for probabilities
    generateLatexProbs(sv, of);

    of << "\\end{document}\n";
    of.close();
    // Compile LaTeX file into PDF
    std::string compileCommand = "pdflatex trial.tex > NUL 2>&1";
    int exitCode = std::system(compileCommand.c_str());

    if (exitCode == 0) {
        std::cout << "PDF generated: trial.pdf" << std::endl;
    } else {
        std::cerr << "Error occurred during compilation." << std::endl;
    }

    return 0;
}


/*#include <iostream>
#include <fstream>

int main() {
    // Create a new ofstream object to write to a .tex file
    std::ofstream outFile("test.tex");

    // Check if the file is successfully opened
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return 1; // Exit with error code
    }

    // Write the LaTeX document content to the file
    outFile << "\\documentclass{article}\n";
    outFile << "\\begin{document}\n";
    outFile << "Hello, LaTeX! This is a test.\n";
    outFile << "\\end{document}\n";

    // Close the file
    outFile.close();

    // Compile the LaTeX document using pdflatex
    system("pdflatex test.tex");

    // Return 0 to indicate successful completion
    return 0;
}
*/
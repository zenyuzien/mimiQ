#include <algorithm>
#include <array>
#include <cmath>
#include <ctime>
#include <functional>
#include <sstream>  // If you use std::istringstream
#include <string>
#include <utility>
#include <vector>

#include "structures.h"
#define R3(x) (std::round((x) * 1000.0) / 1000.0)

mimiqHandler::mimiqHandler(std::string path) {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    dir_path = path;
    wr.open(dir_path + "report.tex");
    wr << "\\documentclass{article}\n";
    wr << "\\usepackage[margin=1in]{geometry}\n";
    wr << "\\usepackage{datetime}\n";
    wr << "\\usepackage{qcircuit}\n";
    wr << "\\usepackage{amsmath}\n";
    wr << "\\usepackage{pgfplots}\n";
    wr << "\\pgfplotsset{compat=1.18}\n";  // Ensure compatibility
    wr << "\\usepackage[utf8]{inputenc}\n";
    wr << "\\begin{document}\n";
    wr << "\\noindent\n\\text{mimiQ++ report - \\today\\ \\currenttime} \n\n";
    circuittDrawn = false;
}
void mimiqHandler::clean() { circuittDrawn = false; }

void mimiqHandler::writeInPdf(std::string msg) {
    if (wr.is_open() && wr.good()) wr << msg << "\n\n";
}

void mimiqHandler::generateReport() {
    if (wr.is_open()) {
        reportGenerated = true;
        wr << "\\end{document}";
        wr.close();
    }

    else
        std::cerr << "Unable to generate report as either already generated earlier (or) writer closed\n";
    std::string latexFilename = dir_path + "report.tex";
    // Compile LaTeX file into PDF
    std::string compileCommand = "pdflatex -output-directory=" + dir_path + " " + latexFilename + " > nul 2>&1";
    int exitCode = std::system(compileCommand.c_str());

    if (exitCode == 0)
        std::cout << "PDF generated: " << latexFilename.substr(0, latexFilename.find_last_of('.')) << ".pdf" << std::endl;
    else
        std::cerr << "Error occurred during latex compilation." << std::endl;

    std::string tmp2 = dir_path + "report.log", tmp3 = dir_path + "report.aux";
    std::remove(tmp2.c_str());
    std::remove(tmp3.c_str());
}

// a qbit as from our understanding can be represented as a sphere will 3 angles which correspond to the rotations from reference axes x,y,z
// in this project, a qbit is represented in the form of a|0> + b|1>
// the reason is for its simplicity and ease to manipulate it mathematically
// a is called the real part, and b is called the complex part (b is imaginary part where a,b together is called complex number, but bare with the naming sense lol)

// prime basis state: |0> , |1>
// for other basis sets: {|+>, |->}, {|i>, |-i>}, spontaneous computation is done as they are not often used

// default constructors

// overloading operators for simplifying arithmetic operations with Coeff type structures
Coeff Coeff::operator*(const Coeff& other) const {
    Coeff result;
    result.real = real * other.real - complex * other.complex;
    result.complex = real * other.complex + complex * other.real;
    return result;
}
Coeff Coeff::operator+(const Coeff& other) const {
    Coeff result;
    result.real = real + other.real;
    result.complex = complex + other.complex;
    return result;
}
Coeff Coeff::operator-(const Coeff& other) const {
    Coeff result;
    result.real = real - other.real;
    result.complex = complex - other.complex;
    return result;
}
Coeff Coeff::operator/(const Coeff& other) const {
    Coeff result;
    result.real = ((real * other.real) + (complex * other.complex)) / ((other.real * other.real) + (other.complex * other.complex));
    result.complex = ((other.complex * real) - (real * other.complex)) / ((other.real * other.real) + (other.complex * other.complex));
    return result;
}

// will return the square of amplitude though the name of the function is amplitude.
double Coeff::amp_sq() const {  // amp squared actually
    return (real * real) + (complex * complex);
}

// our prime basis states as we told are {|0>, |1>}
// so the number of basis states from n qbits is 2^n

// the state_vector is the array of the Coefficients of all the basis states from a system of n_qbits - number of qbits

state_vector::state_vector(int n) {
    n_qbits = n;
    coeffs.resize(1 << n, {0, 0});
}

state_vector::state_vector(const struct state_vector& sv) {
    n_qbits = sv.n_qbits;
    coeffs = sv.coeffs;
}

// appropriate overloading for ease in arithmetic operations
struct state_vector& state_vector::operator=(const struct state_vector& other) {
    // Check for self-assignment
    if (this == &other) {
        return *this;
    }

    // Copy the number of qubits
    n_qbits = other.n_qbits;

    // Deep copy the vector of coefficients
    coeffs = other.coeffs;

    return *this;
}
struct state_vector state_vector::operator*(const struct state_vector& other)  // returns a tensor product of 2 state vectors
{
    struct state_vector sv3;  //(*this).n_qbits == n_qbits
    sv3.n_qbits = n_qbits + other.n_qbits;
    int cnd1 = 1 << n_qbits, cnd2 = 1 << other.n_qbits;
    for (int i = 0; i < cnd1; i++)
        for (int j = 0; j < cnd2; j++) sv3.coeffs.push_back(coeffs[i] * other.coeffs[j]);
    return sv3;
}
void state_vector::print() {
    handler->wr << "\\text{The state vector for the last shot is as follows: }";
    handler->wr << "\\[\n\\begin{array}{@{}llll@{}}\n";
    std::cout << std::endl << "number of qbits = " << n_qbits << std::endl;
    int cnd = 1 << n_qbits;
    for (int i = 0; i < cnd; i++) {
        std::string binaryString;
        for (int j = n_qbits - 1; j >= 0; --j) {
            int bit = (i >> j) & 1;
            binaryString += (bit == 0) ? '0' : '1';
        }
        handler->wr << "\\text{" << binaryString << ":} & " << coeffs[i].amp_sq() * 100 << "\\% & " << coeffs[i].real << " |0\\rangle &  ";
        handler->wr << coeffs[i].complex << " |1\\rangle \\\\" << std::endl;
        std::cout << binaryString << " =>" << coeffs[i].amp_sq() * 100 << "% ( " << coeffs[i].real << " |0> + " << coeffs[i].complex << " i |1> )\n";
    }
    handler->wr << "\\end{array}\n\\]\n";
    std::cout << std::endl;
}

std::pair<Coeff, Coeff> state_vector::measureAlong(int bit) {
    int cnd = 1 << n_qbits;
    Coeff z0, z1;
    auto tmask = 1 << (n_qbits - 1 - bit);
    // std::cout << "mask: "<<tmask<<std::endl;
    for (int i = 0; i < cnd; i++) {
        // std::cout << "i: "<< i << " tmask: "<< tmask << "__ i | tmask: " << (i | tmask ) << "__ i & tmask: "<< (i & tmask) << "__i &~ tmask: "<< (i & (~tmask)) <<std::endl;
        if (i & tmask) {
            z1.real += (coeffs[i].real * coeffs[i].real);
            z1.complex += (coeffs[i].complex * coeffs[i].complex);
        } else {
            z0.real += (coeffs[i].real * coeffs[i].real);
            z0.complex += (coeffs[i].complex * coeffs[i].complex);
        }
    }
    z0.real = std::sqrt(z0.real);
    z0.complex = std::sqrt(z0.complex);
    z1.real = std::sqrt(z1.real);
    z1.complex = std::sqrt(z1.complex);
    return {z0, z1};
}
void state_vector::printprobs() {
    handler->wr << "\\begin{align*}\n";
    int cnd = 1 << n_qbits;
    for (int bit = 0; bit < n_qbits; bit++) {
        auto tmask = 1 << (n_qbits - 1 - bit);
        double prob = 0.0;
        for (int i = 0; i < cnd; i++)
            if (i & tmask) prob += coeffs[i].amp_sq();
        std::cout << "for qbit " << bit << ": " << "Prob of |1> : " << prob << ",  ";
        std::cout << "Prob of |0> : " << 1 - prob << std::endl;
        handler->wr << "\\text{For qbit " << bit << "} \\quad & \\text{Probability of} |1\\rangle: " << prob << ", \\text{Probability of} |0\\rangle: " << 1 - prob << " \\\\\n";
    }
    std::cout << std::endl;
    handler->wr << "\\end{align*}\n";
}

std::string reverse(const std::string& str) {
    std::string reversed = str;                      // Create a copy of the input string
    std::reverse(reversed.begin(), reversed.end());  // Reverse the string
    return reversed;
}
std::string generateLatexBarChart(const std::vector<std::pair<std::string, std::string>>& data) {
    std::ostringstream latex;
    latex << "\\begin{center}\n";
    latex << "    \\begin{tikzpicture}\n";
    latex << "        \\begin{axis}[\n";
    latex << "            ybar,\n";
    latex << "            symbolic x coords={";

    // Add x-coordinates (Quantum States)
    for (size_t i = 0; i < data.size(); ++i) {
        latex << data[i].first;
        if (i != data.size() - 1) latex << ", ";
    }
    latex << "},\n";

    latex << "            xtick=data,\n";
    latex << "            xlabel={Quantum States},\n";
    latex << "            ylabel={Counts},\n";
    latex << "            ymin=0,\n";
    latex << "            bar width=20pt,\n";
    latex << "            width=10cm,\n";
    latex << "            height=7cm,\n";
    latex << "            nodes near coords,\n";
    latex << "            nodes near coords align={vertical},\n";
    latex << "            enlarge x limits=0.3,\n";
    latex << "            title={Quantum State Counts}\n";
    latex << "        ]\n";

    // Add plot coordinates
    latex << "        \\addplot coordinates {";
    for (const auto& [state, count] : data) {
        latex << "(" << state << "," << count << ") ";
    }
    latex << "};\n";

    latex << "        \\end{axis}\n";
    latex << "    \\end{tikzpicture}\n";
    latex << "\\end{center}\n";

    return latex.str();
}

void result::print_counts() {
    std::vector<std::pair<std::string, std::string>> gr;
    handler->wr << "\n\n\n\n\n\n\\text{Classical register readings (left to right: cn,cn-1,..c2,c1,c0) for the simulation:} \n\n\n\n\n\n";
    std::cout << "classical register readings for the simulation: " << std::endl;
    for (auto i = m.begin(); i != m.end(); i++) {
        std::string binaryString;
        for (int j = n_cbits - 1; j >= 0; --j) {
            int bit = (i->first >> j) & 1;
            binaryString += (bit == 0) ? '0' : '1';
        }
        auto rev = reverse(binaryString);
        std::cout << rev << ": " << i->second << std::endl;
        handler->wr << rev << ": " << i->second << "\n\n\n";
        gr.push_back({rev, std::to_string(i->second)});
    }
    handler->wr << generateLatexBarChart(gr);
    std::cout << std::endl;
    handler->wr << "\n\n";
}

std::pair<int, int> result::get_counts_of(int cbit) {
    int cmask = n_cbits - 1 - cbit;
    int ones = 0, zeroes = 0;
    int cnd = 1 << state.n_qbits;

    for (int i = 0; i < cnd; i++) {
        if (i & cmask)  // x1x val = 1
            ones += m[i];
        else
            zeroes += m[i];
    }
    return {zeroes, ones};
}

struct wire {
    std::vector<std::pair<std::string, int>> line;
    int pos;
    bool nature;  // TODO
};

std::string trimZeroes(std::string str) {
    size_t decimalIndex = str.find('.');
    if (decimalIndex != std::string::npos) {
        size_t lastNonZeroIndex = std::string::npos;
        for (auto i = str.size() - 1; i > decimalIndex; --i) {
            if (str[i] != '0') {
                lastNonZeroIndex = i;
                break;
            }
        }
        if (lastNonZeroIndex != std::string::npos) {
            if (lastNonZeroIndex != decimalIndex) {
                return str.substr(0, lastNonZeroIndex + 1);
            } else {
                return str.substr(0, decimalIndex);
            }
        } else {
            return str.substr(0, decimalIndex);
        }
    } else {
        return str;
    }
}

void Qcircuit::printVector() {
    std::cout << "ORDER: \n";
    for (const auto& row : ORDER) {
        for (int element : row) std::cout << element << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Qcircuit::initialize(mimiqHandler* h, int nQ, int nC, std::string n) {
    name = n;
    pdfDone = false;
    n_qbits = nQ;
    n_cbits = nC;
    c_reg = 0;
    handler = h;
    std::vector<struct state_vector> sv;
    sv.resize(nQ);
    for (int i = 0; i < nQ; i++) {
        sv[i].n_qbits = 1;
        sv[i].coeffs.push_back({1, 0});
        sv[i].coeffs.push_back({0, 0});
    }
    if (nQ > 1)
        for (int i = 0; i < nQ - 1; i++) sv[i + 1] = sv[i] * sv[i + 1];
    state = sv[nQ - 1];
    state.handler = handler;
}
Qcircuit::Qcircuit(mimiqHandler* h, int nQ, std::string name = "") { initialize(h, nQ, 0, name); }
Qcircuit::Qcircuit(mimiqHandler* h, int nQ, int nC, std::string name = "") { initialize(h, nQ, nC, name); }

void Qcircuit::print_creg() { std::cout << "creg: " << c_reg << std::endl; }

void Qcircuit::measure(int tbit, int cbit,  // cbit = classical // TODO: int->false
                       int toPrint,
                       int basis)  // 0 basis  = along 0 or 1 , 1 bais = along + or - , 2 basis = along +i or -i
{
    std::pair<Coeff, Coeff> res1;
    if (basis == 1) {
        std::vector<double> placebo;
        applyGate("h", tbit, placebo, false);
        res1 = state.measureAlong(tbit);
        applyGate("h", tbit, placebo, false);
    } else
        res1 = state.measureAlong(tbit);

    if (describe || toPrint) {
        std::cout << "for 0: " << res1.first.real << " " << res1.first.complex << "\n";
        std::cout << "for 1: " << res1.second.real << " " << res1.second.complex << "\n";
    }
    // res = a | 0> + b | 1 > ;
    int cnd = 1 << n_qbits;                  // for 2 qbits , cnd = 4
    auto tmask = 1 << (n_qbits - 1 - tbit);  // qbit: 0 means tmask
                                             // 2-1-0 = 1 shit to left.

    double prob = res1.second.amp_sq();  // b*b
    double randomValue;
    for (int l = 1; l <= 3; l++) randomValue = static_cast<double>(rand()) / RAND_MAX;
    bool res = (randomValue <= prob);
    if (toPrint == 1 || describe) {
        std::cout << "Probability of measuring |1> "
                     "on qubit "
                  << tbit << ": " << prob << std::endl;
        std::cout << "Probability of measuring |0> "
                     "on qubit "
                  << tbit << ": " << 1 - prob << std::endl;
        std::cout << "rand/prob: " << randomValue << "/ " << prob << " res: " << res << " " << (1 << (n_cbits - 1 - cbit)) << std::endl;
    }
    if (res)
        c_reg = c_reg | (1 << (n_cbits - 1 - cbit));
    else
        c_reg = c_reg & (~(1 << (n_cbits - 1 - cbit)));

    // now update state vector to cause collapse
    // TODO: the search is exhaustive, u only need to check limited. Find a better way if possible
    Coeff normaliser = {0, 0};
    if (res)  // if res== 1 , take the zero value and add it to 1 value
    {
        if (describe) std::cout << "measured 1 so.. tmask " << tmask << " \n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
            {
                if (describe) std::cout << "for i =  " << i << "adding to " << (i | (tmask)) << " from " << (i & (~tmask)) << std::endl;

                normaliser.real = normaliser.real + state.coeffs[i | tmask].amp_sq();
                state.coeffs[i & (~tmask)] = {0, 0};
            }
        }
    }

    else  // res == 0
    {
        if (describe) std::cout << "measured 0 so.. tmask is " << tmask << "\n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
            {
                if (describe) std::cout << "for i =  " << i << " adding to " << (i & (~tmask)) << " from " << (i | (tmask)) << std::endl;

                normaliser.real = normaliser.real + state.coeffs[i & (~tmask)].amp_sq();
                state.coeffs[i | tmask] = {0, 0};
            }
        }
    }
    normaliser.real = std::sqrt(normaliser.real);
    if (res)  // if res== 1 , take the zero value and add it to 1 value
    {
        if (describe) std::cout << "measured 1 so.. tmask " << tmask << " \n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
                state.coeffs[i | tmask] = state.coeffs[i | tmask] / normaliser;
        }
    }

    else  // res == 0
    {
        if (describe) std::cout << "measured 0 so.. tmask is " << tmask << "\n";
        for (int i = 0; i < cnd; ++i) {
            if (i & (tmask))  // we got x1x
                state.coeffs[i & (~tmask)] = state.coeffs[i & (~tmask)] / normaliser;
        }
    }
    if (describe) std::cout << "updated creg: " << c_reg << std::endl;

    // std::cout<<"m "<< tbit << " "<<cbit << std::endl;
    ORDER.push_back({3, tbit, cbit, toPrint, basis});
}
void Qcircuit::setCbit1(int cbit) { c_reg = c_reg | (1 << (n_cbits - 1 - cbit)); }
void Qcircuit::setCbit0(int cbit) { c_reg = c_reg & (~(1 << (n_cbits - 1 - cbit))); }

std::pair<std::string, int> Qcircuit::drawCircuitUtility(int key, int eulerIndex) {
    int length = 1;
    // std::cout<<"received with key: "<< key <<" , index: "<< eulerIndex ;
    std::string res;
    auto& euler_array = euler_container[eulerIndex];  // Explicitly reference the array
    switch (key) {
        case 1:
            res = "\\gate{H} & ";
            break;
        case 2:
            res = "\\gate{X} & ";
            break;
        case 3:
            res = "\\gate{Y} & ";
            break;
        case 4:
            res = "\\gate{Z} & ";
            break;
        case 5:
            res = "\\gate{S} & ";
            break;
        case 6:
            res = "\\gate{S^\dag} & ";
            break;
        case 7:
            res = "\\gate{T} & ";
            break;
        case 8:
            res = "\\gate{T^\dag} & ";
            break;
        case 101:
            res = "\\gate{U(" + trimZeroes(std::to_string(R3(euler_array[0]))) + "," + trimZeroes(std::to_string(R3(euler_container[eulerIndex][1]))) + "," + trimZeroes(std::to_string(R3(euler_container[eulerIndex][2]))) + ")} & ";
            length += 3;
            break;
        case 102:
            res = "\\gate{IU3(" + trimZeroes(std::to_string(R3(euler_container[eulerIndex][0]))) + "," + trimZeroes(std::to_string(R3(euler_container[eulerIndex][1]))) + "," + trimZeroes(std::to_string(R3(euler_container[eulerIndex][2]))) + ")} & ";
            length += 3;
            break;
        case 103:
            res = "\\gate{U1(" + trimZeroes(std::to_string(R3(euler_container[eulerIndex][0]))) + ")} & ";
            length++;
            break;
        case 104:
            res = ("\\gate{U2(" + trimZeroes(std::to_string(R3(euler_container[eulerIndex][0]))) + "," + trimZeroes(std::to_string(R3(euler_container[eulerIndex][1]))) + ")} & ");
            length += 2;
            break;

        case 105:
            res = "\\gate{Rx(" + trimZeroes(std::to_string(R3(euler_container[eulerIndex][0]))) + ")} & ";
            length++;
            break;
        case 106:
            res = "\\gate{Ry(" + trimZeroes(std::to_string(R3(euler_container[eulerIndex][0]))) + ")} & ";
            length++;
            break;
    }
    // std::cout<<" returns: "<<res<<std::endl;
    return {res, length};
}
void Qcircuit::drawCircuit(int num = 0) {
    /*
        this function is responsible for drawing the circuit using Qcircuit typesetting langauge (latex)

    */
    // printVector();
    int maxt = 1;
    struct wire qw[n_qbits], cw[n_cbits];
    // init the lines
    for (int i = 0; i < n_qbits; i++) {
        qw[i].line.push_back({("\\lstick{q" + std::to_string(i) + "} & "), 1});
        qw[i].pos = 1;
    }
    for (int i = 0; i < n_cbits; i++) {
        cw[i].line.push_back({("\\lstick{c" + std::to_string(i) + "} & "), 1});
        cw[i].pos = 1;
    }

    // int tmp1, tmp2; // utilities
    for (std::vector<std::vector<int>>::size_type i = 0; i < ORDER.size(); i++) 
    {
        // std::cout<<"current: ";
        // for(auto q : ORDER[i]) std::cout<< q <<" ";
        // std::cout<<std::endl;
        if (ORDER[i][0] == 1) {
            qw[ORDER[i][1]].pos++;
            if (qw[ORDER[i][1]].pos > maxt) {
                maxt = qw[ORDER[i][1]].pos;
            }
            auto tmp = drawCircuitUtility(ORDER[i][2], ORDER[i][3]);
            qw[ORDER[i][1]].line.push_back({tmp.first, tmp.second});

        }
        // TODO OPTIMISE DRAWCIRCUIT
        else if (ORDER[i][0] == 4) {
            // ORDER.push_back ({ 4, t_qbit, G, c_or_q1, control_bit1, control_bit2, esize});
            int tbit = ORDER[i][1];
            int gate = ORDER[i][2];
            int cq = ORDER[i][3];
            int c1 = ORDER[i][4];
            int c2 = ORDER[i][5];
            int eindex = ORDER[i][6];
            /*
                cases:
                1. c c t
                2. c t c
                3. t c c

                top: min of all 3
                bottom: max of all 3
                mid: the middle one

                maxt++
                from top+1 to bottom-1 excluding mid, fill till maxt
                for top,mid,bottom, fill till maxt-1

                add cntl, target for them

            */
            int top = std::min({tbit, c1, c2});
            int bottom = std::max({tbit, c1, c2});
            int mid = tbit + c1 + c2 - top - bottom;

            for (int i = top + 1; i < bottom; i++) {
                if (i == mid) continue;
                while (qw[i].pos <= maxt) {
                    qw[i].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                    qw[i].pos++;
                }
            }
            while (qw[mid].pos < maxt) {
                qw[mid].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                qw[mid].pos++;
            }
            while (qw[top].pos < maxt) {
                qw[top].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                qw[top].pos++;
            }
            while (qw[bottom].pos < maxt) {
                qw[bottom].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                qw[bottom].pos++;
            }

            qw[c1].pos++;
            qw[c1].line.push_back({"\\ctrl{" + std::to_string(tbit - c1) + "} & ", 1});  // Append control gate

            qw[c2].pos++;
            qw[c2].line.push_back({"\\ctrl{" + std::to_string(tbit - c2) + "} & ", 1});  // Append control gate

            qw[tbit].pos++;
            auto tmp = drawCircuitUtility(gate, 0);
            qw[tbit].line.push_back({tmp.first, tmp.second});  // Append the result of drawCircuitUtility
            maxt++;

        } else if (ORDER[i][0] == 2) {
            //  ORDER.push_back ({ 2, t_qbit, GATE, c_or_q, control_bit});
            if (ORDER[i][3] == 1)  // type: q
            {
                // std::cout<<"control: G:  "<<ORDER[i][2] << " from "<<ORDER[i][4]<< " to "<< ORDER[i][1]  <<" \n";
                int low = std::min(ORDER[i][4], ORDER[i][1]);
                int high = (low == ORDER[i][4]) ? ORDER[i][1] : ORDER[i][4];
                // fill the middle ones completely
                // quantum control has cntrl direct single line
                // classical control has doublle line
                for (int x = low + 1; x < high; ++x) {
                    while (qw[x].pos <= maxt) {
                        qw[x].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                        qw[x].pos++;
                    }
                }

                // Fill for the lower bound
                while (qw[low].pos < maxt) {
                    qw[low].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                    qw[low].pos++;
                }

                // Fill for the higher bound
                while (qw[high].pos < maxt) {
                    qw[high].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                    qw[high].pos++;
                }

                // Add control gate for ORDER[i][4]
                qw[ORDER[i][4]].line.push_back({
                    "\\ctrl{" + std::to_string(ORDER[i][1] - ORDER[i][4]) + "} & ",
                    1  // Length increment
                });
                qw[ORDER[i][4]].pos++;

                // Add the gate for ORDER[i][1]
                qw[ORDER[i][1]].pos++;
                auto tmp = drawCircuitUtility(ORDER[i][2], 0);
                qw[ORDER[i][1]].line.push_back({tmp.first, tmp.second});  // Append the result of drawCircuitUtility
                maxt++;

            } else if (ORDER[i][3] == 0) {
                // ORDER.push_back ({ 2, t_qbit, G, c_or_q, control_bit, esize});

                // [control, target] - add wires till the max
                for (int x = ORDER[i][1]; x < n_qbits; x++) {
                    while (qw[x].pos < maxt) {
                        qw[x].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                        qw[x].pos++;
                    }
                }

                for (int x = 0; x < n_cbits; x++) {
                    while (cw[x].pos < maxt) {
                        cw[x].line.push_back({"\\cw & ", 1});  // Append "\\cw & " with a length increment of 1
                        cw[x].pos++;
                    }
                }

                // The wires between get extra "\\qw \\cwx"
                for (int x = ORDER[i][1] + 1; x < n_qbits; x++) {
                    qw[x].line.push_back({"\\qw \\cwx & ", 1});  // Append "\\qw \\cwx & " with a length increment of 1
                    qw[x].pos++;
                }

                // The control and target get specific control/wire commands
                for (int x = 0; x < n_cbits; x++) {
                    if (x != ORDER[i][4]) {
                        if (x < ORDER[i][4]) {
                            cw[x].line.push_back({"\\cw \\cwx & ", 1});  // Append "\\cw \\cwx & " with a length increment of 1
                        } else {
                            cw[x].line.push_back({"\\cw & ", 1});  // Append "\\cw & " with a length increment of 1
                        }
                        cw[x].pos++;
                    } else {
                        cw[x].line.push_back({"\\control \\cw \\cwx & ", 1});  // Append "\\control \\cw \\cwx & " with a length increment of 1
                        cw[x].pos++;
                    }
                }

                qw[ORDER[i][1]].pos++;
                auto tmp = drawCircuitUtility(ORDER[i][2], ORDER[i][5]);
                qw[ORDER[i][1]].line.push_back({tmp.first, tmp.second});  // Append the result of drawCircuitUtility
                maxt++;
            }
        } else if (ORDER[i][0] == 3)  // measure
        {
            // For the main, put wires
            while (qw[ORDER[i][1]].pos < maxt) {
                qw[ORDER[i][1]].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                qw[ORDER[i][1]].pos++;
            }
            qw[ORDER[i][1]].line.push_back({"\\meter & ", 1});  // Append "\\meter & " with a length increment of 1
            qw[ORDER[i][1]].pos++;

            // Q: Putting wires in between
            for (int x = ORDER[i][1] + 1; x < n_qbits; x++) {
                while (qw[x].pos < maxt) {
                    qw[x].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                    qw[x].pos++;
                }
            }

            // Q: Edgy wires for the between
            for (int x = ORDER[i][1] + 1; x < n_qbits; x++) {
                qw[x].line.push_back({"\\qw \\cwx & ", 1});  // Append "\\qw \\cwx & " with a length increment of 1
                qw[x].pos++;
            }

            // C: Putting cwires
            for (int x = 0; x <= ORDER[i][2]; ++x) {
                while (cw[x].pos < maxt) {
                    cw[x].line.push_back({"\\cw & ", 1});  // Append "\\cw & " with a length increment of 1
                    cw[x].pos++;
                }
            }

            // C: Putting edgy wires
            for (int x = 0; x <= ORDER[i][2]; ++x) {
                cw[x].line.push_back({"\\cw \\cwx & ", 1});  // Append "\\cw \\cwx & " with a length increment of 1
                cw[x].pos++;
            }

            maxt++;
        }
    }
        for (int i = 0; i < n_qbits; ++i) {
            while (qw[i].pos <= maxt) {
                qw[i].line.push_back({"\\qw & ", 1});  // Append "\\qw & " with a length increment of 1
                qw[i].pos++;
            }
            qw[i].line.push_back({"\\\\ ", 0});  // Add a newline in LaTeX format
        }

        for (int i = 0; i < n_cbits; ++i) {
            while (cw[i].pos <= maxt) {
                cw[i].line.push_back({"\\cw & ", 1});  // Append "\\cw & " with a length increment of 1
                cw[i].pos++;
            }
            cw[i].line.push_back({"\\\\ ", 0});  // Add a newline in LaTeX format

        }
    
        // std::cout << "\\[\n\\Qcircuit @C=1em @R=.7em {\n";
        //
        handler->wr << "\\clearpage\n\\begin{figure}[htbp]\n\\[\n\\Qcircuit @C=1em @R=.7em {\n";
        //int max_length = 0;











        std::vector<int> inds;
        int csum= 0 ;
        inds.push_back(0);
        int starter=0;

            for(int z = starter ; z < qw[0].line.size(); z++)
            {
                int max_block = -1; 
                for (int i = 0; i < n_qbits; i++) 
                {
                    if(max_block < qw[i].line[z].second)
                        max_block = qw[i].line[z].second;
                }
                csum += max_block ;
                if(csum>15)
                {
                    inds.push_back(z);
                    csum=0;
                }
            }
            
        
        std::cout<<"finally \n";
        // 0 5 10 17
        uint32_t min_diff = 0xffffffff; 
        for(int i = 0 ; i < inds.size()-1; i++)
        {
            min_diff = ( min_diff > (inds[i+1]-inds[i]) )?((inds[i+1]-inds[i])):(min_diff);
        }
        int tmp = inds.size();
        inds.clear();
        for (int i = 0; i <= qw[0].line.size(); i += min_diff) {
            inds.push_back(i);
            std::cout<<i<<" ";
        }
        std::cout<<std::endl;




        /*std::vector<std::pair< std::string, int>> panel;
        int csum= 0 ;
        uint32_t maxi =0; 
        inds.push_back(0);
        int starter=0;

        again:;
        if(maxi == 0xffffffff)goto end;
        maxi=0xffffffff;
        std::cout<<"starter: "<<starter<<"\n";
        for (int i = 0; i < n_qbits; i++) {
            std::cout << "\nin "<< i << std::endl;
            
            csum=0;
            for(int z = starter ; z < qw[0].line.size(); z++)
            {
                csum += qw[i].line[z].second ;
                std::cout<<csum<<" ";
                if( csum>15 )
                {
                    std::cout<<" index: "<<z<<std::endl;
                    if(maxi>z)
                    maxi = z;
                    std::cout<<" maxi: "<<maxi<<std::endl;
                    break;
                }
            
            }
        }
        std::cout<<"\nend of iteration, \n";
        starter = maxi+1;
        if(starter)
        inds.push_back(starter);
        goto again;
        end:;
        std::cout<<"finally \n";
        for(auto& i: inds)
        std::cout<<i<<" ";
        std::cout<<std::endl;*/

        for(int i = 0 ; i < inds.size()-1; i++) // i about iterations
        {
            for(int j = 0 ; j < n_qbits ; j++)  // j about q/c number 
            {
                std::string combinedLine= qw[j].line[0].first ;
                for(int k = inds[i]; k < inds[i+1]; k++)    // k about line indicies
                    if(k)
                    combinedLine += qw[j].line[k].first;

                handler->wr << combinedLine<< " \\qw & \\\\"  << std::endl ;
                
            }
            for(int j = 0 ; j < n_cbits ; j++)
            {
                std::string combinedLine= cw[j].line[0].first ;
                for(int k = inds[i]; k < inds[i+1]; k++)
                    if(k)
                    combinedLine += cw[j].line[k].first;

                handler->wr << combinedLine << " \\cw & \\\\" << std::endl ;
            }
            handler->wr<<"\\\\ \n\\\\ \n";

        }
            for(int i = 0 ; i < n_qbits ; i++)
            {
                std::string combinedLine= qw[i].line[0].first ;
                for(int j = inds[inds.size()-1]; j < qw[0].line.size(); j++)
                    if(j)
                    combinedLine += qw[i].line[j].first;

                handler->wr << combinedLine  << std::endl ;
            }
            for(int i = 0 ; i < n_cbits ; i++)
            {
                std::string combinedLine= cw[i].line[0].first ;
                for(int j = inds[inds.size()-1]; j < qw[0].line.size(); j++)
                    if(j)
                    combinedLine += cw[i].line[j].first;

                handler->wr << combinedLine  << std::endl ;
            }
            handler->wr<<"\\\\ \n\\\\ \n";
                

        /*


            std::string combinedLine= qw[i].line[0].first ;
            int ut=0, help=0;
            for (size_t j = 1; j < qw[i].line.size(); ++j) {

                help+= csum[j];
                if(help>15)
                {
                    help=0;
                    panel.push_back({ combinedLine ,ut++});
                    combinedLine= qw[i].line[0].first ;
                }
                combinedLine += qw[i].line[j].first;  // Access each string using the index
            }
            if(combinedLine.size())
            panel.push_back({ combinedLine ,ut++});
            if(max_ut < ut)max_ut = ut;
            //handler->wr << combinedLine << std::endl;  // Write the combined line to the handler
        */

        // std::cout << "}\n\\end{displaymath}\n";
        if (handler->wr.good()) {
            handler->wr << "}\n\\]\n";
            if (name != "") handler->wr << "\\caption{" << name << "}\n";
            handler->wr << "\\end{figure}\n";
            handler->wr.flush();
        } else
            std::cerr << "writing issue";

         //std::cout<<"CIRCUIT maxt: "<<maxt << " len: "<< max_length << std::endl;
         //handler->writeInPdf("CIRCUIT maxt: "+ std::to_string(maxt)+ " len: " + std::to_string(max_length));
    
}

bool Qcircuit::accessCreg(int cbit) { return (c_reg & (1 << (n_cbits - 1 - cbit))); }

/*void clear(bool isHardreset = false)
                      {
                        if (isHardreset) // hard reset
                            {
                                n_qbits = n_cbits = order_ptr = 0;
                                state.n_qbits = 0;
                                state.coeffs.clear ();
                                for (auto &inner_vector : ORDER)
                                    inner_vector.clear ();
                                ORDER.clear ();
                            }
                        else // back to ZERO
                            {
                                order_ptr = 0;
                                state = ZERO;
                            }
                        c_reg = 0;
                        euler_container.clear ();
                    }*/

// std::vector<struct experiment> simulation()
Experiment Qcircuit::simulation() {
    if (handler->circuittDrawn == false) {
        drawCircuit();
        handler->circuittDrawn = true;
    }
    Experiment exp1;
    exp1.c_reg_value = c_reg;
    exp1.final_state = state;
    exp1.n_cbits = static_cast<uint64_t>(n_cbits);
    if (!handler->wr.is_open()) std::cerr << "already closed before returning shot";
    // handler->wr = NULL;
    return exp1;
}

struct result simulate(mimiqHandler* handler, std::function<Qcircuit::experiment(mimiqHandler*)> func, int shots) {
    if (!func) {
        std::cerr << "NULL Experiment error \n";
        exit(0);
    }
    struct result res;
    Experiment shot_result;

    while (shots--) {
        shot_result = func(handler);
        if (!handler->wr.is_open()) std::cerr << "end of shot closed";
        res.m[shot_result.c_reg_value]++;
    }

    // last shot results
    res.n_cbits = shot_result.n_cbits;

    res.creg.resize(res.n_cbits);
    for (uint64_t i = 0; i < res.n_cbits; ++i) res.creg[shot_result.n_cbits - 1 - i] = (shot_result.c_reg_value >> i) & 1;
    res.state = shot_result.final_state;  // TODO final state doesnt make sense
    // res.state.print();
    res.handler = handler;
    return res;
}

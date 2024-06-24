#ifndef _structures
#define _structures

#define root2 std::sqrt(2)
#define PI acos(0.0) * 2

struct mimiqHandler 
{
    std::ofstream wr; // to write Latex 
    std::string dir_path ; // where pdf should be stored
    bool circuittDrawn ;

    mimiqHandler(std::string path = "")
    {
        if(GLOBAL_SEEDER == false)
        {
            std::srand(static_cast < unsigned int > (std::time(nullptr)));
            GLOBAL_SEEDER = true;
        }

        dir_path = path ; 
        wr.open(dir_path + "report.tex");
        wr<<"\\documentclass{article}\n\\usepackage[margin=1in]{geometry}\n\\usepackage{datetime}\n\\usepackage{qcircuit}\n\\usepackage{amsmath}\n\\begin{document}\n";
        wr<<"\\noindent\n\\text{mimiQ++ report - \\today\\ \\currenttime} \\\\\n\\\\\n";
        circuittDrawn = false;
    }

    void writeInPdf(std::string msg)
    {
        if(wr.is_open() && wr.good())
        wr << msg << "\\\\";
    } 
    void generateReport();
    

};
// a qbit as from our understanding can be represented as a sphere will 3 angles which correspond to the rotations from reference axes x,y,z
// in this project, a qbit is represented in the form of a|0> + b|1> 
// the reason is for its simplicity and ease to manipulate it mathematically
// a is called the real part, and b is called the complex part (b is imaginary part where a,b together is called complex number, but bare with the naming sense lol)

struct Coeff 
{
    // prime basis state: |0> , |1>
    // for other basis sets: {|+>, |->}, {|i>, |-i>}, spontaneous computation is done as they are not often used
    double real;
    double complex;

    // default constructors
    Coeff(): real(0), complex(0) {}
    Coeff(double r, double c): real(r), complex(c) {}

    // overloading operators for simplifying arithmetic operations with Coeff type structures
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

    // will return the square of amplitude though the name of the function is amplitude. 
    double amp_sq() const { //amp squared actually
        return (real * real) + (complex * complex);
    }
};

static std::vector<std::vector<Coeff> > empty_matrix ; // for applyUtility
// our prime basis states as we told are {|0>, |1>}
// so the number of basis states from n qbits is 2^n

// the state_vector is the array of the Coefficients of all the basis states from a system of n_qbits - number of qbits
struct state_vector {
    int n_qbits; 
    std::vector < Coeff > coeffs; // 16 * (1<< nq) // for 2 qbits -> 544 cbits, for 3qbits -> 1056 cbits
    mimiqHandler* handler; 
    // default constuctor
    state_vector(): n_qbits(0) {}
    state_vector(int n) {
        n_qbits = n;
        coeffs.resize(1 << n, { 0, 0 });
    }
    
    state_vector(const struct state_vector & sv) {
        n_qbits = sv.n_qbits;
        coeffs = sv.coeffs;
    }

    // appropriate overloading for ease in arithmetic operations
    struct state_vector& operator=(const struct state_vector& other) 
    {
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
    struct state_vector operator * (const struct state_vector & other) // returns a tensor product of 2 state vectors 
    {
        struct state_vector sv3; //(*this).n_qbits == n_qbits 
        sv3.n_qbits = n_qbits + other.n_qbits;
        int cnd1 = 1 << n_qbits, cnd2 = 1 << other.n_qbits;
        for(int i = 0; i < cnd1; i++)
            for(int j = 0; j < cnd2; j++) sv3.coeffs.push_back(coeffs[i] * other.coeffs[j]);
        return sv3;
    }
    void print() 
    {
        handler->wr<< "\\[\n\\begin{array}{@{}llll@{}}\n";
        std::cout << std::endl << "number of qbits = " << n_qbits << std::endl;
        int cnd = 1 << n_qbits;
        for(int i = 0; i < cnd; i++) {
            std::string binaryString;
            for(int j = n_qbits - 1; j >= 0; --j) {
                int bit = (i >> j) & 1;
                binaryString += (bit == 0) ? '0' : '1';
            }
           handler->wr << "\\text{" << binaryString << ":} & " << coeffs[i].amp_sq() * 100 << "\\% & " << coeffs[i].real << " |0\\rangle &  " ;
           handler->wr << coeffs[i].complex << " |1\\rangle \\\\" << std::endl;
           std::cout << binaryString << " =>" << coeffs[i].amp_sq() * 100 << "% ( " << coeffs[i].real << " |0> + " << coeffs[i].complex << " i |1> )\n";
        }
        handler->wr << "\\end{array}\n\\]\n";
        std::cout << std::endl;
    }
    
    std::pair<Coeff,Coeff> measureAlong(int bit)
    {
        int cnd = 1 << n_qbits ;
        Coeff z0,z1; 
        auto tmask = 1 << (n_qbits - 1 - bit);
        // std::cout << "mask: "<<tmask<<std::endl; 
        for(int i = 0; i < cnd; i++)
        {
            // std::cout << "i: "<< i << " tmask: "<< tmask << "__ i | tmask: " << (i | tmask ) << "__ i & tmask: "<< (i & tmask) << "__i &~ tmask: "<< (i & (~tmask)) <<std::endl;
            if(i & tmask)
            {
                z1.real += ( coeffs[i ].real * coeffs[i ].real );
                z1.complex += ( coeffs[i ].complex * coeffs[i ].complex );
            }
            else 
            {
                z0.real += ( coeffs[i ].real * coeffs[i ].real );
                z0.complex += ( coeffs[i ].complex * coeffs[i ].complex );
            }
        }
        z0.real = std::sqrt(z0.real);
        z0.complex = std::sqrt(z0.complex);
        z1.real = std::sqrt(z1.real);
        z1.complex= std::sqrt(z1.complex);
        return {z0,z1};
    }
    void printprobs()
    {
        handler->wr << "\\begin{align*}\n";
        int cnd = 1 << n_qbits ;
        for(int bit = 0; bit <  n_qbits ; bit++)
        {
            auto tmask = 1 << (n_qbits - 1 - bit);
            double prob = 0.0;
            for(int i = 0; i < cnd; i++)
                if(i & tmask) prob += coeffs[i].amp_sq();
            std::cout << "for qbit " << bit << ": " <<"Prob of |1> : " << prob << ",  ";
            std::cout << "Prob of |0> : " << 1 - prob << std::endl;
            handler->wr << "\\text{For qbit "<< bit <<"} \\quad & \\text{Probability of} |1\\rangle: " << prob << ", \\text{Probability of} |0\\rangle: "
            << 1 - prob << " \\\\\n";
        }
        std::cout<<std::endl;
        handler->wr<<"\\end{align*}\n";
    }
};
struct result 
{
    uint64_t n_cbits;
    std::map < int, int > m; // TODO rename, restructure 
    std::vector<bool> creg;
    struct state_vector state;
    mimiqHandler* handler; 
    
    void print_counts() {
        handler->wr<< "\\\\\n\\\\\n\\text{Classical register readings (left to right, c0 least significant) for the simulation:} \\\\\n\\\\\n";
        std::cout << "classical register readings for the simulation: " << std::endl;
        for(auto i = m.begin(); i != m.end(); i++) {
            std::string binaryString;
            for(int j = n_cbits - 1; j >= 0; --j) {
                int bit = (i -> first >> j) & 1;
                binaryString += (bit == 0) ? '0' : '1';
            }
            std::cout << binaryString << ": " << i -> second << std::endl;  
            handler->wr << binaryString << ": "<< i -> second << "\\\\\n";
        }
        std::cout << std::endl;
        handler->wr<< "\\\\";
    }
    
    std::pair<int,int> get_counts_of(int cbit)
    {
        int cmask = n_cbits - 1 - cbit ; 
        int ones =0 , zeroes = 0 ;
        int cnd = 1 << state.n_qbits;

        for(int i = 0 ; i < cnd ; i++ )
        {
            if(i & cmask ) // x1x val = 1 
            ones += m[i];   
            else 
            zeroes += m[i];
                
        }
        return {zeroes,ones};
    }
    
};

struct Qcircuit
{
    bool pdfDone, describe = false; // TODO deal with describe
    int n_qbits, n_cbits;
    uint64_t c_reg;
    std::vector<std::vector<int> > ORDER; // for drawing circuits
    std::vector<std::array<double, 3> > euler_container; // for drawing circuits
    struct state_vector state;
    struct experiment
    {
        uint64_t c_reg_value,n_cbits ;
        struct state_vector final_state;
    };
    mimiqHandler* handler;    
    std::string name; 

    void printVector();
    void initialize(mimiqHandler* h, int nQ, int nC, std::string name );
    Qcircuit ( mimiqHandler* h, int nQ, int nC, std::string );
    Qcircuit ( mimiqHandler* h, int nQ, std::string );
    ~Qcircuit()
    {
       // delete handler; TODO analyse
    }
    void applyUtility(const std::string &gate, int& G, int& esize,std::vector<double> arr, std::vector<std::vector<Coeff> >& matrix);
    void applyGate (const std::string &gate, int t_qbit, std::vector<double> arr, bool ORDERinc ); // overloading applyGate
    void print_state();
    void print_creg ();
    void applyControlledGate (const std::string &gate, const std::string c_qbit_str,int t_qbit, std::vector<double> arr );
    void measure (int tbit, int cbit, // cbit = classical // TODO: int->false
             int toPrint ,
             int basis ); // 0 basis  = along 0 or 1 , 1 bais = along + or - , 2 basis = along +i or -i 
    void setCbit1 (int cbit);
    void setCbit0 (int cbit);
    std::string drawCircuitUtility(int key, int eulerIndex, int& length );
    void drawCircuit(int num);
    bool accessCreg(int cbit);
    experiment simulation();
};
typedef Qcircuit::experiment Experiment;

void mimiqHandler::generateReport()
{
    if (wr.is_open()) 
    {
        wr<<"\\end{document}";
        wr.close();
    }
    else std::cerr<<"Unable to generate report: ALready closed \n";
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

#endif
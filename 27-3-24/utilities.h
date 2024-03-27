#ifndef _util
#define _util
bool getPdf = false; 
bool getPdfUtil = false;
void ketZero(struct state_vector & sv) {
    sv.n_qbits = 1;
    sv.coeffs.push_back({ 1, 0 });
    sv.coeffs.push_back({ 0, 0 });
}
void ketONE(struct state_vector & sv) {
    sv.n_qbits = 1;
    sv.coeffs.push_back({ 0, 0 });
    sv.coeffs.push_back({ 1, 0 });
}
void cleanupFunction()
{
    if (wr.is_open()) 
    {
        wr<<"\\end{document}";
        wr.close();
    }
    if(getPdf)
    {
        std::string latexFilename = ABSPATH + "qcirc.tex";
        // Compile LaTeX file into PDF
        //std::string compileCommand = "pdflatex " + latexFilename + " > nul 2>&1";
        std::string compileCommand = "pdflatex -output-directory=" + ABSPATH + " " + latexFilename + " > nul 2>&1";
        int exitCode = std::system(compileCommand.c_str());

        if (exitCode == 0) {
            std::cout << "PDF generated: " << latexFilename.substr(0, latexFilename.find_last_of('.')) << ".pdf" << std::endl;
        } else {
            std::cerr << "Error occurred during latex compilation." << std::endl;
        }
        
        if(!getPdfUtil)
        {
            std::string tmp2 = ABSPATH + "qcirc.log", tmp3 = ABSPATH + "qcirc.aux";
            std::remove(tmp2.c_str());
            std::remove(tmp3.c_str());
        }
        std::remove(latexFilename.c_str());
    }

}
struct GlobalInitializer {
    GlobalInitializer() {
        std::cout << "mimiQ++ \n";
        std::srand(static_cast < unsigned int > (std::time(nullptr)));
        atexit(cleanupFunction);

        wr<<"\\documentclass{article}\n\\usepackage{qcircuit}\n\\usepackage{amsmath}\n\\begin{document}\n";
        /*
        \documentclass{article}
        \usepackage{qcircuit}
        \usepackage{amsmath}
        \begin{document}
        */
    }
}
initializer;

void get_pdf(int num=0)
{
    if(num ==1 )
    {
        //keep temp files 
        getPdfUtil = true; 
    }
    getPdf = true ; 
}

#endif
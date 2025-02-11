#include <bits/stdc++.h>
std::string ABSPATH = "/home/zenyuzien/Downloads/" ;
int main()
{
        std::string latexFilename = ABSPATH + "report.tex";
        // Compile LaTeX file into PDF
        //std::string compileCommand = "pdflatex " + latexFilename + " > nul 2>&1";
        std::string compileCommand = "pdflatex -output-directory=" + ABSPATH + " " + latexFilename + " > nul 2>&1";
        int exitCode = std::system(compileCommand.c_str());

        if (exitCode == 0) {
            std::cout << "PDF generated: " << latexFilename.substr(0, latexFilename.find_last_of('.')) << ".pdf" << std::endl;
        } else {
            std::cerr << "Error occurred during latex compilation." << std::endl;
        }
        return 0;
}
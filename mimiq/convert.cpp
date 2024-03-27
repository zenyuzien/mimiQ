#include <iostream>
#include <cstdlib> // for system()

int main() {
    // Provide the filename of the LaTeX file to convert
    std::string latexFilename = "trial.tex";

    // Compile LaTeX file into PDF
    std::string compileCommand = "pdflatex " + latexFilename + " > nul 2>&1";

    int exitCode = std::system(compileCommand.c_str());

    if (exitCode == 0) {
        std::cout << "PDF generated: " << latexFilename.substr(0, latexFilename.find_last_of('.')) << ".pdf" << std::endl;
    } else {
        std::cerr << "Error occurred during compilation." << std::endl;
    }

    return 0;
}

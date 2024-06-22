#ifndef _util
#define _util
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

struct GlobalInitializer {
    GlobalInitializer() {
        std::cout << "mimiQ++ \n";
    }
}
initializer;

std::string trimZeroes(std::string str) {
    size_t decimalIndex = str.find('.');
    if (decimalIndex != std::string::npos) {
        // Find the last non-zero character index after the decimal point
        size_t lastNonZeroIndex = std::string::npos;
        for (auto i = str.size() - 1; i > decimalIndex; --i) {
            if (str[i] != '0') {
                lastNonZeroIndex = i;
                break;
            }
        }

        if (lastNonZeroIndex != std::string::npos) {
            // If there are non-zero digits after the decimal point
            if (lastNonZeroIndex != decimalIndex) {
                return str.substr(0, lastNonZeroIndex + 1);
            } else {
                return str.substr(0, decimalIndex);
            }
        } else {
            // If all digits after the decimal point are zero
            return str.substr(0, decimalIndex);
        }
    } else {
        // If there is no decimal point
        return str;
    }
}

#endif
#ifndef _qcirc
#define _qcirc

struct wire
{
    std::string line;
    int pos;
};

struct Qcircuit
{

    int n_qbits, n_cbits, order_ptr;
    uint64_t c_reg;
    struct state_vector state, ZERO;

    std::vector<std::vector<int> > ORDER;
    std::vector<std::array<double, 3> > euler_container;

    void
    printVector ()
    {
        std::cout << "ORDER: \n";
        for (const auto &row : ORDER)
            {
                for (int element : row)
                    {
                        std::cout << element << " ";
                    }
                std::cout << std::endl;
            }
        std::cout << std::endl;
    }
    Qcircuit (int nQ, int nC = 0)
        : n_qbits (nQ), n_cbits (nC), order_ptr (0), c_reg (0)
    {
        std::vector<struct state_vector> sv;
        sv.resize (nQ);
        for (int i = 0; i < nQ; i++)
            ketZero (sv[i]);
        if (nQ > 1)
            for (int i = 0; i < nQ - 1; i++)
                sv[i + 1] = sv[i] * sv[i + 1];
        ZERO = state = sv[nQ - 1];
        ORDER.push_back ({ 4, true }); // false means hard reset
    }
    Qcircuit (std::string path) : order_ptr (0), c_reg (0)
    {
        std::vector<struct state_vector> sv;
        ORDER = downloadCircuit (ABSPATH + "ORDERS/" + path);
        sv.resize (n_qbits);
        for (int i = 0; i < n_qbits; i++)
            ketZero (sv[i]);
        if (n_qbits > 1)
            for (int i = 0; i < n_qbits - 1; i++)
                sv[i + 1] = sv[i] * sv[i + 1];
        ZERO = state = sv[n_qbits - 1];
        // ORDER.push_back({4,true}); // false means hard reset
    }
    void
    applyGate (const std::string &gate, int t_qbit)
    {
        int GATE;
        if (gate == "H" || gate == "h")
            GATE = 1;
        else if (gate == "x" || gate == "X")
            GATE = 2;
        else if (gate == "z" || gate == "Z")
            GATE = 4;
        else if (gate == "y" || gate == "Y")
            GATE = 3;
        ORDER.push_back ({ 1, t_qbit, GATE });
        return;
    }
    void
    applyGate (const std::string &gate, int t_qbit,
               std::vector<double> arr) // overloading applyGate
    {
        if (gate == "U3" || gate == "u3" || gate == "u" || gate == "U")
            {
                int size = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], arr[2] });
                ORDER.push_back ({ 1, t_qbit, 5, size }); // send index
            }
        else if (gate == "IU3" || gate == "iu3" || gate == "IU"
                 || gate == "iu") // todo: generalise u1 u2 u3 inverse
            {
                int size = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], arr[2] });
                ORDER.push_back ({ 1, t_qbit, 6, size }); // send index
            }
        else if (gate == "U1" || gate == "u1")
            {
                int size = euler_container.size ();
                euler_container.push_back ({ arr[0], 0, 0 });
                ORDER.push_back ({ 1, t_qbit, 7, size }); // send index
            }
        else if (gate == "U2" || gate == "u2")
            {
                int size = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], 0 });
                ORDER.push_back ({ 1, t_qbit, 8, size }); // send index
            }
    }
    void
    print_state ()
    {
        ORDER.push_back ({ 9 });
    }
    void
    print_creg ()
    {
        ORDER.push_back ({ 11 });
    }

    void
    applyControlledGate (const std::string &gate, const std::string c_qbit_str,
                         int t_qbit)
    {
        std::istringstream iss (c_qbit_str);
        char cq;
        int n;
        iss >> cq;
        iss >> n;
        int third = (cq == 'q') ? 1 : 0;
        if (gate == "H" || gate == "h")
            ORDER.push_back ({ 2, 1, third, n, t_qbit });
        else if (gate == "x" || gate == "X")
            ORDER.push_back ({ 2, 2, third, n, t_qbit });
        else if (gate == "z" || gate == "Z")
            ORDER.push_back ({ 2, 4, third, n, t_qbit });
    }
    void
    measure (int tbit, int cbit,
             int toPrint = 0) // cbit = classical // TODO: int->false
    {
        ORDER.push_back ({ 3, tbit, cbit, toPrint });
    }
    void
    setCbit1 (int cbit)
    {
        ORDER.push_back ({ 10, 1, cbit });
    }
    void
    setCbit0 (int cbit)
    {
        ORDER.push_back ({ 10, 2, cbit });
    }
    void
    uploadCircuit (const std::vector<std::vector<int> > &ORDER,
                   const std::string &path)
    {
        std::string full = ABSPATH + "/ORDERS/" + path;
        std::ofstream outfile (full);
        outfile << n_qbits << " " << n_cbits << " " << ORDER.size ()
                << std::endl;
        if (outfile.is_open ())
            {
                for (const auto &inner : ORDER)
                    {
                        for (int val : inner)
                            {
                                outfile << val << " ";
                            }
                        outfile << "\n";
                    }
                for (const auto &array : euler_container)
                    {
                        for (const auto &value : array)
                            {
                                outfile << value << " "; // Write each value
                                                         // followed by a space
                            }
                        outfile << "\n"; // Write newline after each array
                    }
                std::cout << "circuit saved! \n";
                outfile.close ();
            }
        else
            std::cerr << "Error: Unable to open file " << full << "\n";
    }
    std::vector<std::vector<int> >
    downloadCircuit (const std::string &path)
    {
        int ordersize;
        std::vector<std::vector<int> > result;
        std::ifstream infile (path);
        if (infile.is_open ())
            {
                std::string line;
                std::getline (infile, line);
                std::istringstream iss0 (line);
                iss0 >> n_qbits;
                iss0 >> n_cbits;
                iss0 >> ordersize;
                while (ordersize--)
                    {
                        std::getline (infile, line);
                        std::istringstream iss (line);
                        int val;
                        std::vector<int> inner;
                        while (iss >> val)
                            inner.push_back (val);
                        result.push_back (inner);
                    }
                double theta, lam, psi;
                while (std::getline (infile, line))
                    {
                        std::istringstream iss1 (line);
                        iss1 >> theta;
                        iss1 >> lam;
                        iss1 >> psi;
                        euler_container.push_back ({ theta, lam, psi });
                    }
                infile.close ();
            }
        else
            std::cerr << "Error: Unable to open file " << path << ".txt\n";
        return result;
    }

    struct result simulate (int shots, std::string path, bool describe);
    void printQ ();

    void get_pdf(int num=0)
    {
        if(num ==1 )
        {
            //keep temp files 
            getPdfUtil = true; 
        }
        getPdf = true ; 
        printQ();
    }
};

void
Qcircuit::printQ ()
{
    struct wire qw[n_qbits], cw[n_cbits];
    for (int i = 0; i < n_qbits; i++)
        {
            qw[i].line += ("\\lstick{q" + std::to_string (i) + "} & ");
            qw[i].pos = 1;
        }
    for (int i = 0; i < n_cbits; i++)
        {
            cw[i].line += ("\\lstick{c" + std::to_string (i) + "} & ");
            cw[i].pos = 1;
        }

    int tmp1, tmp2; // utilities
    for (std::vector<std::vector<int> >::size_type i = 0; i < ORDER.size ();
         i++)
        {
            if (ORDER[i][0] == 1)
                {
                    int tmp1 = ORDER[i][1], tmp2 = ORDER[i][2];
                pdf_gater:;
                    qw[tmp1].pos++;
                    switch (tmp2)
                        {
                        case 1:
                            qw[tmp1].line += ("\\gate{H} & ");
                            break;
                        case 2:
                            qw[tmp1].line += ("\\gate{X} & ");
                            break;
                        case 3:
                            qw[tmp1].line += ("\\gate{Y} & ");
                            break;
                        case 4:
                            qw[tmp1].line += ("\\gate{Z} & ");
                            break;
                        case 5:
                            qw[tmp1].line
                                += ("\\gate{U("
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][0]))
                                    + ","
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][1]))
                                    + ","
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][2]))
                                    + ")} & ");
                            break;
                        case 6:
                            qw[tmp1].line
                                += ("\\gate{IU3("
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][0]))
                                    + ","
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][1]))
                                    + ","
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][2]))
                                    + ")} & ");
                            break;
                        case 7:
                            qw[tmp1].line
                                += ("\\gate{U1("
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][0]))
                                    + ")} ");
                            break;
                        case 8:
                            qw[tmp1].line
                                += ("\\gate{U2("
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][0]))
                                    + ","
                                    + trimZeroes (std::to_string (
                                        euler_container[ORDER[i][3]][1]))
                                    + ")} & ");
                            break;
                        }
                }
            else if (ORDER[i][0] == 2)
                {
                    if (ORDER[i][2] == 1) // type: q
                        {
                            int low = std::min (ORDER[i][3], ORDER[i][4]);
                            int high = (low == ORDER[i][3]) ? ORDER[i][4]
                                                            : ORDER[i][3];
                            int maxt = 1;

                            for (int x = low; x <= high; ++x)
                                if (maxt < qw[x].pos)
                                    maxt = qw[x].pos;

                            for (int x = low + 1; x < high; ++x)
                                while (qw[x].pos <= maxt)
                                    {
                                        qw[x].line += "\\qw & ";
                                        qw[x].pos++;
                                    }

                            while (qw[low].pos < maxt)
                                {
                                    qw[low].line += "\\qw & ";
                                    qw[low].pos++;
                                }

                            while (qw[high].pos < maxt)
                                {
                                    qw[high].line += "\\qw & ";
                                    qw[high].pos++;
                                }

                            qw[ORDER[i][3]].line
                                += ("\\ctrl{"
                                    + std::to_string (ORDER[i][4] - ORDER[i][3])
                                    + "} & ");
                            qw[ORDER[i][3]].pos++;
                            qw[ORDER[i][4]].pos++;

                            if (ORDER[i][1] == 1) // gate
                                qw[ORDER[i][4]].line += ("\\gate{H} & ");

                            else if (ORDER[i][1] == 2)
                                qw[ORDER[i][4]].line += ("\\gate{X} & ");

                            else if (ORDER[i][1]
                                     == 4) // 4 - pauli z , 3 - pauli y
                                qw[ORDER[i][4]].line += ("\\gate{Z} & ");
                        }
                    else if (ORDER[i][2] == 0)
                        {

                            ; // work needs to be done
                        }
                }
            else if (ORDER[i][0] == 3) // measure
                {
                    // finding max
                    int maxt = 1;
                    for (int x = ORDER[i][1]; x < n_qbits; x++)
                        if (maxt < qw[x].pos)
                            maxt = qw[x].pos;

                    // for the main, put wires
                    while (qw[ORDER[i][1]].pos < maxt)
                        {
                            qw[ORDER[i][1]].line += "\\qw & ";
                            qw[ORDER[i][1]].pos++;
                        }
                    qw[ORDER[i][1]].line += "\\meter & ";
                    qw[ORDER[i][1]].pos++;

                    // q: putting wires in between
                    for (int x = ORDER[i][1] + 1; x < n_qbits; x++)
                        while (qw[x].pos < maxt)
                            {
                                qw[x].line += ("\\qw & ");
                                qw[x].pos++;
                            }
                    // q: edgy wires for the between
                    for (int x = ORDER[i][1] + 1; x < n_qbits; x++)
                        {
                            qw[x].line += ("\\qw \\cwx & ");
                            qw[x].pos++;
                        }
                    // c: putting wires
                    for (int x = 0; x <= ORDER[i][2]; ++x)
                        while (cw[x].pos < maxt)
                            {
                                cw[x].line += ("\\cw & ");
                                cw[x].pos++;
                            }
                    // c: putting edgy wires
                    for (int x = 0; x <= ORDER[i][2]; ++x)
                        {
                            cw[x].line += ("\\cw \\cwx & ");
                            cw[x].pos++;
                        }
                }
        }
    int maxt = 1;
    for (int i = 0; i < n_qbits; ++i)
        {
            if (qw[i].pos > maxt)
                maxt = qw[i].pos;
        }
    for (int i = 0; i < n_qbits; ++i)
        {
            while (qw[i].pos <= maxt)
                {
                    qw[i].line += "\\qw & ";
                    qw[i].pos++;
                }
            qw[i].line += "\\\\ ";
        }
    for (int i = 0; i < n_cbits; ++i)
        {
            while (cw[i].pos <= maxt)
                {
                    cw[i].line += "\\cw & ";
                    cw[i].pos++;
                }
            cw[i].line += "\\\\ ";
        }

    std::cout << "\\begin{displaymath}\n\\Qcircuit @C=1em @R=.7em {\n";
    wr << "\\begin{displaymath}\n\\Qcircuit @C=1em @R=.7em {\n";
    for (int i = 0; i < n_qbits; i++)
        {
            std::cout << qw[i].line << std::endl;
            wr << qw[i].line << std::endl;
        }

    for (int i = 0; i < n_cbits; i++)
        {
            std::cout << cw[i].line << std::endl;
            wr << cw[i].line << std::endl;
        }

    std::cout << "}\n\\end{displaymath}\n";
    wr << "}\n\\end{displaymath}\n";
}

struct result
Qcircuit::simulate (int shots = 1, std::string path = "", bool describe = false)
{

    if (describe)
        shots = 1;
    if (describe)
        std::cout << "simulation starting with shots: " << shots << std::endl;
    if (path != "")
        {
            if (describe)
                std::cout << "uploading circuit to " << path << std::endl;
            uploadCircuit (ORDER, path);
        }
    std::map<int, int> m;
    int tmp1, tmp2;
    while (shots--)
        {
            for (std::vector<std::vector<int> >::size_type i = 0;
                 i < ORDER.size (); i++)
                {
                    if (ORDER[i][0] == 9) // print state vector
                        {
                            if (describe)
                                std::cout << i << ") printing state vector"
                                          << std::endl;
                            state.print ();
                        }
                    else if (ORDER[i][0] == 11) // print creg
                        {
                            if (describe)
                                std::cout << i << ") printing creg"
                                          << std::endl;
                            std::cout << "creg: " << c_reg << std::endl;
                        }
                    else if (ORDER[i][0] == 1) // applygate
                        {
                            if (describe)
                                std::cout << i << ") applygate" << std::endl;
                            tmp1 = ORDER[i][1]; // qbit
                            tmp2 = ORDER[i][2]; // gate
                        gater:;

                            std::vector<std::vector<Coeff> > matrix,
                                buffedmatrix;
                            if (describe)
                                std::cout << "to qbit and gate " << tmp1 << ", "
                                          << tmp2 << std::endl;
                            switch (tmp2)
                                {
                                case 1:
                                    matrix = GATES::unitaryGate (PI / 2, 0, PI);
                                    break;
                                case 2:
                                    matrix = GATES::unitaryGate (PI, 0, PI);
                                    break;
                                case 3:
                                    matrix = GATES::unitaryGate (PI, PI / 2,
                                                                 PI / 2);
                                    break;
                                case 4:
                                    matrix = GATES::unitaryGate (0, 0, PI);
                                    break;
                                case 5:
                                    matrix = GATES::unitaryGate (
                                        euler_container[ORDER[i][3]][0],
                                        euler_container[ORDER[i][3]][1],
                                        euler_container[ORDER[i][3]][2]);
                                    break;
                                case 6:
                                    matrix = GATES::inverse2x2 (
                                        GATES::unitaryGate (
                                            euler_container[ORDER[i][3]][0],
                                            euler_container[ORDER[i][3]][1],
                                            euler_container[ORDER[i][3]][2]));
                                    break;
                                case 7:
                                    matrix = GATES::unitaryGate (
                                        0, 0, euler_container[ORDER[i][3]][0]);
                                    break;
                                case 8:
                                    matrix = GATES::unitaryGate (
                                        PI, euler_container[ORDER[i][3]][0],
                                        euler_container[ORDER[i][3]][1]);
                                    break;
                                }
                            GATES::kroneckerBUFF (matrix, buffedmatrix, n_qbits,
                                                  tmp1);
                            int n = 1 << n_qbits;
                            std::vector<Coeff> result (n, { 0, 0 });
                            for (int ii = 0; ii < n; ++ii)
                                {
                                    for (int j = 0; j < n; ++j)
                                        result[ii] = result[ii]
                                                     + (buffedmatrix[ii][j]
                                                        * state.coeffs[j]);
                                }
                            state.coeffs = result;
                        }
                    else if (ORDER[i][0] == 2) // controlledgate
                        {
                            /*
                                ORDER[i][1] = gate
                                ORDER[i][2] = type q/c
                                ORDER[i][3] = control bit
                                ORDER[i][4] = target qbit

                            */
                            if (describe)
                                std::cout << i << ") apply controlled gate "
                                          << std::endl;
                            if (ORDER[i][2] == 1) // type: q
                                {
                                    if (describe)
                                        std::cout << "type:q with gate "
                                                  << ORDER[i][1]
                                                  << "control and targets: "
                                                  << ORDER[i][3] << " "
                                                  << ORDER[i][4] << std::endl;
                                    if (ORDER[i][1] == 1) // gate
                                        GATES::controlled_hadamard (
                                            state, ORDER[i][3],
                                            ORDER[i][4]); // qbit, tqbit

                                    else if (ORDER[i][1] == 2)
                                        GATES::controlled_pauli_X (
                                            state, ORDER[i][3], ORDER[i][4]);

                                    else if (ORDER[i][1]
                                             == 4) // 4 - pauli z , 3 - pauli y
                                        GATES::controlled_pauli_Z (
                                            state, ORDER[i][3], ORDER[i][4]);
                                }
                            else if (ORDER[i][2] == 0)
                                {
                                    if (describe)
                                        std::cout
                                            << "control by cbit.. creg: "
                                            << c_reg << " mask: "
                                            << (1
                                                << (n_cbits - 1 - ORDER[i][3]))
                                            << std::endl;
                                    if (c_reg
                                        & (1 << (n_cbits - 1 - ORDER[i][3])))
                                        {
                                            // std::cout<<"going to gater
                                            // "<<std::endl ;
                                            tmp2 = ORDER[i][1];
                                            tmp1 = ORDER[i][4];
                                            goto gater;
                                        }
                                    else if (describe)
                                        std::cout << " cbit control false \n";
                                }
                        }
                    else if (ORDER[i][0] == 3) // measure
                        {
                            if (describe)
                                std::cout << i << ") measuring " << ORDER[i][1]
                                          << " \n"
                                          << std::endl;
                            int cnd = 1 << n_qbits; // for 2 qbits , cnd = 4
                            auto tmask
                                = 1
                                  << (n_qbits - 1
                                      - ORDER[i][1]); // qbit: 0 means tmask
                                                      // 2-1-0 = 1 shit to left.
                            double prob = 0.0;
                            for (int i = 0; i < cnd; i++)
                                if (i & tmask)
                                    prob += state.coeffs[i].amplitude ();
                            double randomValue
                                = static_cast<double> (rand ()) / RAND_MAX;
                            bool res = (randomValue <= prob);
                            if (ORDER[i][3] == 1 || describe)
                                {
                                    std::cout << "Probability of measuring |1> "
                                                 "on qubit "
                                              << ORDER[i][1] << ": " << prob
                                              << std::endl;
                                    std::cout << "Probability of measuring |0> "
                                                 "on qubit "
                                              << ORDER[i][1] << ": " << 1 - prob
                                              << std::endl;
                                    std::cout
                                        << "rand/prob: " << randomValue << "/ "
                                        << prob << " res: " << res << " "
                                        << (1 << (n_cbits - 1 - ORDER[i][2]))
                                        << std::endl;
                                }
                            if (res)
                                c_reg = c_reg
                                        | (1 << (n_cbits - 1 - ORDER[i][2]));
                            else
                                c_reg = c_reg
                                        & (~(1 << (n_cbits - 1 - ORDER[i][2])));
                            if (describe)
                                std::cout << "updated creg: " << c_reg
                                          << std::endl;
                        }
                    else if (ORDER[i][0] == 4) // clear
                        {
                            if (describe)
                                std::cout << i << ") clearing with hardreset: "
                                          << !ORDER[i][1] << std::endl;
                            if (ORDER[i][1] == 0) // hard reset
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
                        }
                    else if (ORDER[i][0] == 10) // set c bit
                        {
                            if (describe)
                                std::cout << i
                                          << ") set c bit given bit position: "
                                          << ORDER[i][2] << " which is mask "
                                          << (1 << (n_cbits - 1 - ORDER[i][2]))
                                          << std::endl;
                            if (ORDER[i][1] == 1)
                                c_reg = c_reg
                                        | (1 << (n_cbits - 1 - ORDER[i][2]));
                            else if (ORDER[i][1] == 2)
                                c_reg = c_reg
                                        & (~(1 << (n_cbits - 1 - ORDER[i][2])));
                            if (describe)
                                std::cout << "currently c_reg: " << c_reg
                                          << std::endl;
                        }
                }
            m[c_reg]++;
        }

    return { n_cbits, m, state };
}

#endif
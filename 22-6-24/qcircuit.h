#ifndef _qcirc
#define _qcirc

#include "utilities.h"
struct wire
{
    std::string line;
    int pos;
    bool nature; // TODO
};
    void Qcircuit::printVector ()
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
    Qcircuit::Qcircuit (mimiqHandler* h, int nQ, int nC = 0)
        : pdfDone(false), n_qbits (nQ), n_cbits (nC), c_reg (0)
    {
        handler = h ;        
        std::vector<struct state_vector> sv;
        sv.resize (nQ);
        for (int i = 0; i < nQ; i++)
            ketZero (sv[i]);
        if (nQ > 1)
            for (int i = 0; i < nQ - 1; i++)
                sv[i + 1] = sv[i] * sv[i + 1];
        state = sv[nQ - 1];
        state.handler=handler;
       
    }
    
    void Qcircuit::applyUtility(const std::string &gate, int& G= dummy, int& esize = dummy ,std::vector<double> arr = std::vector<double>(),std::vector<std::vector<Coeff> >& matrix = empty_matrix)
    {
        if (gate == "H" || gate == "h")
        {
            G = 1;
            matrix = GATES::unitaryGate (PI / 2, 0, PI);
        }
            
        else if (gate == "x" || gate == "X")
        {
            G = 2;
            matrix = GATES::unitaryGate (PI, 0, PI);
        }
            
        else if (gate == "z" || gate == "Z")
        {
            G = 4;
            matrix = GATES::unitaryGate (0, 0, PI);
        }
            
        else if (gate == "y" || gate == "Y")
        {
            G = 3;
            matrix = GATES::unitaryGate (PI, PI / 2, PI / 2);
        }
            
        else if (gate == "U3" || gate == "u3" || gate == "u" || gate == "U")
            {
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], arr[2] });
                G= 5;
                matrix = GATES::unitaryGate (
                                        euler_container[esize][0],
                                        euler_container[esize][1],
                                        euler_container[esize][2]);
            }
        else if (gate == "IU3" || gate == "iu3" || gate == "IU"
                 || gate == "iu") // todo: generalise u1 u2 u3 inverse
            {
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], arr[2] });
                G = 6;
                matrix = GATES::inverse2x2 (
                                        GATES::unitaryGate (
                                            euler_container[esize][0],
                                            euler_container[esize][1],
                                            euler_container[esize][2]));
            }
        else if (gate == "U1" || gate == "u1")
            {
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], 0, 0 });
                G=7;
                matrix = GATES::unitaryGate (
                                        0, 0, euler_container[esize][0]);
            }
        else if (gate == "U2" || gate == "u2")
            {
                esize = euler_container.size ();
                euler_container.push_back ({ arr[0], arr[1], 0 });
                G =8;
                matrix = GATES::unitaryGate (
                                        PI/2, euler_container[esize][0],
                                        euler_container[esize][1]);
            }
    }

    void
    Qcircuit::applyGate (const std::string &gate, int t_qbit,
               std::vector<double> arr = std::vector<double>(), bool ORDERinc = true)
    {
        std::vector<std::vector<Coeff> > matrix,buffedmatrix;
        int esize = -1,G;
   
        applyUtility(gate,G,esize,arr,matrix);
            
            //std::cout<<"p: "<<1 <<" "<<t_qbit << " " << G << " " << esize << std::endl;
            if(ORDERinc)
            ORDER.push_back ({ 1, t_qbit, G, esize }); // send index
                            
                            GATES::kroneckerBUFF (matrix, buffedmatrix, n_qbits,
                                                  t_qbit);
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
    void
    Qcircuit::print_state ()
    {
        state.print();
        //ORDER.push_back ({ 9 });
    }
    void
    Qcircuit::print_creg ()
    {
        std::cout << "creg: " << c_reg << std::endl;
        //ORDER.push_back ({ 11 });
    }

    void
    Qcircuit::applyControlledGate (const std::string &gate, const std::string c_qbit_str,
                         int t_qbit,  std::vector<double> arr = std::vector<double>())
    {
        // control = true => CONTROL ON ZERO
        // control = false => CONTROL ON ONE
        std::istringstream iss (c_qbit_str);
        char cq;
        int control_bit,c_or_q,G;
        iss >> cq;
        iss >> control_bit;
        c_or_q = (cq == 'q') ? 1 : 0;

        if(c_or_q ==0)
        {
                                    if(c_reg & (1 << (n_cbits - 1 - control_bit)))
                                            applyGate(gate,t_qbit,arr,false);
        }

        int esize;
        applyUtility(gate,G,esize,arr);

        ORDER.push_back ({ 2, t_qbit, G, c_or_q, control_bit, esize});
        
        if(c_or_q == 1)
        { // CONTROL ON ZERO WITH Q needs to be done. TODO
                                    if (G == 1) // gate
                                        GATES::controlled_hadamard (
                                            state, control_bit,
                                            t_qbit); // qbit, tqbit

                                    else if (G== 2)
                                        GATES::controlled_pauli_X (
                                            state, control_bit, t_qbit);

                                    else if (G
                                             == 4) // 4 - pauli z , 3 - pauli y
                                        GATES::controlled_pauli_Z (
                                            state, control_bit, t_qbit);
        }
    }
    void
    Qcircuit::measure (int tbit, int cbit, // cbit = classical // TODO: int->false
             int toPrint = 0,
             int basis = 0) // 0 basis  = along 0 or 1 , 1 bais = along + or - , 2 basis = along +i or -i 
    {

                             //    if (describe || toPrint)
                              //  std::cout << i << ") measuring " << tbit
                                   //       << " with basis state " <<basis <<"\n"
                                     //     << std::endl;

                            auto res1 = state.measureAlong(tbit,basis);
                            if (describe || toPrint)
                            {
                                std::cout << "for 0: " << res1.first.real << " " << res1.first.complex << "\n";
                                std::cout << "for 1: " << res1.second.real << " " << res1.second.complex << "\n";
                            }
                            // res = a | 0> + b | 1 > ; 
                            int cnd = 1 << n_qbits; // for 2 qbits , cnd = 4
                            auto tmask
                                = 1
                                  << (n_qbits - 1
                                      - tbit); // qbit: 0 means tmask
                                                      // 2-1-0 = 1 shit to left.
                            // res = a | 0> + b | 1 > ; 
                            
                            double prob = res1.second.amp_sq(); // b*b
                            double randomValue;
                            for(int l = 1 ; l <= 3; l++) randomValue
                                = static_cast<double> (rand ()) / RAND_MAX;
                            bool res = (randomValue <= prob);
                            if (toPrint == 1 || describe)
                                {
                                    std::cout << "Probability of measuring |1> "
                                                 "on qubit "
                                              << tbit << ": " << prob
                                              << std::endl;
                                    std::cout << "Probability of measuring |0> "
                                                 "on qubit "
                                              << tbit << ": " << 1 - prob
                                              << std::endl;
                                    std::cout
                                        << "rand/prob: " << randomValue << "/ "
                                        << prob << " res: " << res << " "
                                        << (1 << (n_cbits - 1 - cbit))
                                        << std::endl;
                                }
                            if (res)
                                c_reg = c_reg
                                        | (1 << (n_cbits - 1 - cbit));
                            else
                                c_reg = c_reg
                                        & (~(1 << (n_cbits - 1 - cbit)));

                            // now update state vector to cause collapse
                            // TODO: the search is exhaustive, u only need to check limited. Find a better way if possible
                            Coeff normaliser = {0,0};
                            if(res)  // if res== 1 , take the zero value and add it to 1 value
                            {
                                if(describe)
                                std::cout << "measured 1 so.. tmask " <<tmask <<" \n";
                                for(int i = 0 ; i < cnd ; ++i)
                                {
                                    if( i & (tmask) ) // we got x1x
                                    {
                                        if(describe)
                                        std::cout <<"for i =  " << i <<  "adding to "<<  (i | (tmask)) << " from " <<  (i & (~tmask)) << std::endl ;
                                    
                                        normaliser.real = normaliser.real + state.coeffs[i | tmask].amp_sq();
                                        state.coeffs[i & (~tmask)]= {0,0};
                                    }
                                }
                            }
                            
                            else // res == 0 
                            {
                                if(describe)
                                std::cout << "measured 0 so.. tmask is " << tmask << "\n";
                                for(int i = 0 ; i < cnd ; ++i)
                                {
                                    
                                    if( i & (tmask) ) // we got x1x
                                    {
                                        if(describe)
                                        std::cout << "for i =  " << i <<  " adding to "<< ( i & (~tmask) ) << " from " <<  (i | (tmask)) << std::endl ;
                                    
                                    normaliser.real = normaliser.real + state.coeffs[i & (~tmask)].amp_sq();
                                    state.coeffs[i | tmask]={0,0};
                                    }
                                }
                            }  
                            normaliser.real = std::sqrt(normaliser.real);
                            if(res)  // if res== 1 , take the zero value and add it to 1 value
                            {
                                if(describe)
                                std::cout << "measured 1 so.. tmask " <<tmask <<" \n";
                                for(int i = 0 ; i < cnd ; ++i)
                                {
                                    if( i & (tmask) ) // we got x1x
                                    state.coeffs[i | tmask] = state.coeffs[i | tmask] / normaliser;
                                }
                            }
                            
                            else // res == 0 
                            {
                                if(describe)
                                std::cout << "measured 0 so.. tmask is " << tmask << "\n";
                                for(int i = 0 ; i < cnd ; ++i)
                                {
                                    
                                    if( i & (tmask) ) // we got x1x
                                    state.coeffs[i & (~tmask)] = state.coeffs[i & (~tmask)] / normaliser ;
                                    
                                }
                            }  
                            if (describe)
                                std::cout << "updated creg: " << c_reg
                                          << std::endl;

        //std::cout<<"m "<< tbit << " "<<cbit << std::endl;
        ORDER.push_back ({ 3, tbit, cbit, toPrint, basis });
    }
    void
    Qcircuit::setCbit1 (int cbit)
    {
        c_reg = c_reg | (1 << (n_cbits - 1 - cbit));
    }
    void
    Qcircuit::setCbit0 (int cbit)
    {
        c_reg = c_reg & (~(1 << (n_cbits - 1 - cbit)));
    }
    
    std::string Qcircuit::drawCircuitUtility(int key, int eulerIndex =  -1 )
    {
        std::cout<<"received with key: "<< key <<" , index: "<< eulerIndex ; 
        std::string res;
                            switch (key)
                            {
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
                                res = "\\gate{U("
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][0]))
                                        + ","
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][1]))
                                        + ","
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][2]))
                                        + ")} & ";
                                break;
                            case 6:
                                res = "\\gate{IU3("
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][0]))
                                        + ","
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][1]))
                                        + ","
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][2]))
                                        + ")} & ";
                                break;
                            case 7:
                                res = "\\gate{U1("
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][0]))
                                        + ")} ";
                                break;
                            case 8:
                                res = ("\\gate{U2("
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][0]))
                                        + ","
                                        + trimZeroes (std::to_string (
                                            euler_container[eulerIndex][1]))
                                        + ")} & ");
                                break;
                            }
        std::cout<<" returns: "<<res<<std::endl;
        return res;
    }
    void Qcircuit::drawCircuit(int num=0)
    {
        printVector();
        handler->wr<< "\\text{Circuit diagram:} \\\\\n\\\\\n";
        int maxt= 1 ;
        struct wire qw[n_qbits], cw[n_cbits];
        // init the lines
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

    // int tmp1, tmp2; // utilities
        for (std::vector<std::vector<int> >::size_type i = 0; i < ORDER.size ();
            i++)
            {
                std::cout<<"current: ";
                for(auto q : ORDER[i]) std::cout<< q <<" ";
                std::cout<<std::endl;
                if (ORDER[i][0] == 1)
                {
                    qw[ORDER[i][1]].pos++;
                    if(qw[ORDER[i][1]].pos > maxt) maxt= qw[ORDER[i][1]].pos; 
                    
                                if(ORDER[i][2] < 5)
                                qw[ORDER[i][1]].line += drawCircuitUtility(ORDER[i][2]);
                                else 
                                qw[ORDER[i][1]].line += drawCircuitUtility(ORDER[i][2],ORDER[i][3]);
                }
                // TODO OPTIMISE DRAWCIRCUIT
                else if (ORDER[i][0] == 2)
                    {
                    //  ORDER.push_back ({ 2, t_qbit, GATE, c_or_q, control_bit});
                        if (ORDER[i][3] == 1) // type: q
                            {
                                //std::cout<<"control: G:  "<<ORDER[i][2] << " from "<<ORDER[i][4]<< " to "<< ORDER[i][1]  <<" \n"; 
                                int low = std::min (ORDER[i][4], ORDER[i][1]);
                                int high = (low == ORDER[i][4]) ? ORDER[i][1]
                                                                : ORDER[i][4];
                                
                                // fill the middle ones completely
                                // quantum control has cntrl direct single line
                                // classical control has doublle line 
                                for (int x = low + 1; x < high; ++x)
                                    while (qw[x].pos <= maxt)
                                        {
                                            qw[x].line += "\\qw & ";
                                            qw[x].pos++;
                                        }
                                // the lower and higher fill 
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

                                qw[ORDER[i][4]].line
                                    += ("\\ctrl{"
                                        + std::to_string (ORDER[i][1] - ORDER[i][4])
                                        + "} & ");
                                qw[ORDER[i][4]].pos++;

                                qw[ORDER[i][1]].pos++;
                                qw[ORDER[i][1]].line += drawCircuitUtility(ORDER[i][2]);
                                maxt++;
                            }
                        else if (ORDER[i][3] == 0)
                            {
                                //ORDER.push_back ({ 2, t_qbit, G, c_or_q, control_bit, esize});

                               // [control, target] - add wires till the max
                                for (int x = ORDER[i][1]; x < n_qbits; x++)
                                    while(qw[x].pos < maxt)
                                    {
                                        qw[x].line += "\\qw & ";
                                        qw[x].pos ++; 
                                    }
                                for (int x = 0; x < n_cbits; x++)
                                    while(cw[x].pos < maxt)
                                    {
                                        cw[x].line += "\\cw & ";
                                        cw[x].pos ++; 
                                    }   
                               // the between get extra wires
                                for(int x = ORDER[i][1] +1; x < n_qbits; x++)
                                {
                                    qw[x].line += "\\qw \\cwx & ";
                                    qw[x].pos ++; 
                                }
                                //the between get extra wires, the control gets control, the target gates gate. 
                                for (int x = 0; x < n_cbits; x++)
                                {
                                    if(x != ORDER[i][4])
                                    {
                                        if(x < ORDER[i][4])
                                        cw[x].line += "\\cw \\cwx & ";
                                        else 
                                        cw[x].line += "\\cw & ";
                                        cw[x].pos ++; 
                                    }
                                    else 
                                    {
                                        cw[x].line += "\\control \\cw \\cwx & ";
                                        cw[x].pos ++; 
                                    }
                                } 
                                qw[ORDER[i][1]].pos++;
                                std::cout<<"before: "<< qw[ORDER[i][1]].line << std::endl;
                                qw[ORDER[i][1]].line += drawCircuitUtility(ORDER[i][2],ORDER[i][5]);
                                std::cout<<"after: "<< qw[ORDER[i][1]].line << std::endl;                             
                                maxt++;
                            }
                    }
                else if (ORDER[i][0] == 3) // measure
                    {
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
                        // c: putting cwires
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
                        maxt++;
                    }
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

        //std::cout << "\\begin{displaymath}\n\\Qcircuit @C=1em @R=.7em {\n";
        handler->wr << "\\begin{displaymath}\n\\Qcircuit @C=1em @R=.7em {\n";
        for (int i = 0; i < n_qbits; i++)
            {
                //std::cout << qw[i].line << std::endl;
                handler->wr << qw[i].line << std::endl;
            }

        for (int i = 0; i < n_cbits; i++)
            {
                //std::cout << cw[i].line << std::endl;
                handler->wr << cw[i].line << std::endl;
            }

        //std::cout << "}\n\\end{displaymath}\n";
        if(handler->wr.good())
        {
            handler->wr << "}\n\\end{displaymath}\n";
            handler->wr.flush();
        }
        else 
            std::cerr<<"writing issue";
        
    }

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

    //std::vector<struct experiment> simulation()
    Experiment Qcircuit::simulation()
    {
        if(handler->circuittDrawn == false)
        {
            drawCircuit();
            handler->circuittDrawn = true;
        }
        Experiment exp1;
        exp1.c_reg_value = c_reg;
        exp1.final_state = state;
        exp1.n_cbits = static_cast<uint64_t>(n_cbits);
                   if(!handler->wr.is_open())
                   std::cerr<<"already closed before returning shot";
                //handler->wr = NULL;
        return exp1;
    }


/*
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
          m[c_reg]++;
        }

    return { n_cbits, m, state };
}
*/

    struct result simulate(mimiqHandler* handler, std::function<Qcircuit::experiment(mimiqHandler* )> func = nullptr , int shots = 1 )
    {
        if(!func)
        {
            std::cerr << "NULL Experiment error \n";
            exit(0);
        }
        struct result res;
        Experiment shot_result;
        
        while(shots--)
        {
            shot_result = func(handler);
                if(!handler->wr.is_open())
                std::cerr<<"end of shot closed";
            res.m[shot_result.c_reg_value]++;
        }
        
        // last shot results
        res.n_cbits = shot_result.n_cbits;
        
        res.creg.resize(res.n_cbits); 
        for (uint64_t i = 0; i < res.n_cbits; ++i) 
            res.creg[shot_result.n_cbits - 1 - i] = (shot_result.c_reg_value >> i) & 1;
        res.state = shot_result.final_state; // TODO final state doesnt make sense
        
        res.handler = handler;
            return res;
        }

    

#endif

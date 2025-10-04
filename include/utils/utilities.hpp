#ifndef UTILITIES_H
#define UTILITIES_H

#include "typedefs.hpp"
#include "space_efficient_vector.hpp"

// void printRecompressionRLSLP(const RecompressionRLSLP *recompression_rlslp) {

//     cout << "RECOMPRESSION PRINTING STARTED..." << endl;

//     if(!(recompression_rlslp->nonterm).size()) {
//         cout << "RecompressionRLSLP is Empty!!" << endl;
//         return;
//     }

//     c_size_t i = 0;

//     for(const RLSLPNonterm & rlslp_nonterm : recompression_rlslp->nonterm) {
//         cout << i << " --> ";
//         cout << rlslp_nonterm.type << ' ' << rlslp_nonterm.first << ' ' << rlslp_nonterm.second <<  ' ' << rlslp_nonterm.explen << endl;
//         i++;
//     }

//     cout << "RECOMPRESSION PRINTING ENDED!" << endl << endl;

//     return;
// }

// void printSLG(const SLG *slg) {
//     for(int i=0; i<slg->nonterm.size(); i++) {
//         cout << i << " : ";
//         c_size_t start_index = slg->nonterm[i].start_index;
//         c_size_t end_index = (i==slg->nonterm.size()-1) ? slg->rhs.size()-1 : slg->nonterm[i+1].start_index - 1;

//         for(int j=start_index; j<=end_index; j++) {
//             cout << slg->rhs[j] << ' ';
//         }

//         cout << endl;
//     }

//     return;
// }

void expandRLSLP(c_size_t var, const space_efficient_vector<RLSLPNonterm> & rlslp_nonterm_vec, space_efficient_vector<c_size_t> & result) {
    if(rlslp_nonterm_vec[var].type == '0') {
        result.push_back(rlslp_nonterm_vec[var].first);
        return;
    }

    char_t type = rlslp_nonterm_vec[var].type;
    c_size_t first = rlslp_nonterm_vec[var].first;
    c_size_t second = rlslp_nonterm_vec[var].second;

    if(type == '1') {
        expandRLSLP(first, rlslp_nonterm_vec, result);
        expandRLSLP(second, rlslp_nonterm_vec, result);
    }
    // Type == '2'
    else {
        while(second--) {
            expandRLSLP(first, rlslp_nonterm_vec, result);
        }
    }
    return;
}

void expandRLSLP(const RecompressionRLSLP *recompression_rlslp, space_efficient_vector<c_size_t> &result) {

    const space_efficient_vector<RLSLPNonterm> &rlslp_nonterm_vec = recompression_rlslp->nonterm;
    expandRLSLP(rlslp_nonterm_vec.size()-1, rlslp_nonterm_vec, result);

    return;
}

space_efficient_vector<c_size_t> expandRLSLP(const RecompressionRLSLP *recompression_rlslp) {

    space_efficient_vector<c_size_t> result;
    space_efficient_vector<RLSLPNonterm> rlslp_nonterm_vec = recompression_rlslp->nonterm;
    expandRLSLP(rlslp_nonterm_vec.size()-1, rlslp_nonterm_vec, result);


    return result;
}

InputSLP* getSLP(c_size_t grammar_size) {

    InputSLP* slp = new InputSLP();
    space_efficient_vector<SLPNonterm> &nonterm = slp->nonterm;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    srand(seed);



    // Generate random number X < n/2
    c_size_t num_terminals = rand() % ((c_size_t)((3*grammar_size)/4));

    num_terminals++;

    for(c_size_t i=1; i<=num_terminals; i++) {
        nonterm.push_back(SLPNonterm('0', -i, 1));
    }


    for(c_size_t i = num_terminals; i < grammar_size; ++i) {
        c_size_t num1 = rand() % i; // Generate first number less than i
        c_size_t num2 = rand() % i; // Generate second number less than i


        c_size_t bit1 = rand()%2;
        c_size_t bit2 = rand()%2;

        c_size_t num = (bit2 == 1 ? num2 : num1);

        if(bit1 == 0)
        nonterm.push_back(SLPNonterm('1', i-1, num));
        else {
            nonterm.push_back(SLPNonterm('1', num, i-1));
        }
    }

    return slp;
}

void print_current_timestamp() {
    // current date/time based on current system
    time_t now = time(0);

    // convert now to string form
    char* dt = ctime(&now);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);
    cout << endl << "Generated on [UTC Timestamp] "<< dt << endl;;
}

// void expandSLG(c_size_t var, const space_efficient_vector<SLGNonterm> & slg_nonterm_vec, space_efficient_vector<c_size_t> & result) {
//     if(var < 0) {
//         result.push_back(abs(var));
//         return;
//     }

//     const space_efficient_vector<c_size_t> & rhs = slg_nonterm_vec[var].rhs;

//     for(c_size_t rhs_var : rhs) {
//         expandSLG(rhs_var, slg_nonterm_vec, result);
//     }
// }

// space_efficient_vector<c_size_t> expandSLG(const SLG *slg) {
//     space_efficient_vector<c_size_t> result;

//     const space_efficient_vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

//     expandSLG(slg_nonterm_vec.size() - 1, slg_nonterm_vec, result);

//     return result;
// }
#endif

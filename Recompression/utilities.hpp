// #ifndef UTILITIES_H
// #define UTILITIES_H

// #include "typedefs.hpp"

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
//     cout << "SLG PRINTING STARTED..." << endl;

//     if(!(slg->nonterm).size()) {
//         cout << "SLG is Empty!!" << endl;
//         return;
//     }

//     c_size_t i = 0;

//     for(const SLGNonterm & slg_nonterm : slg->nonterm) {
//         cout << i << " --> ";
//         for(c_size_t x : slg_nonterm.rhs) {
//             cout << x << ' ';
//         }
//         cout << endl;
//         i++;
//     }

//     cout << "SLG PRINTING ENDED!" << endl << endl;

//     return;
// }

void expandRLSLP(c_size_t var, const vector<RLSLPNonterm> & rlslp_nonterm_vec, vector<c_size_t> & result) {
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
vector<c_size_t> expandRLSLP(const RecompressionRLSLP *recompression_rlslp) {

    vector<c_size_t> result;
    vector<RLSLPNonterm> rlslp_nonterm_vec = recompression_rlslp->nonterm;
    expandRLSLP(rlslp_nonterm_vec.size()-1, rlslp_nonterm_vec, result);


    return result;
}

// void expandSLG(c_size_t var, const vector<SLGNonterm> & slg_nonterm_vec, vector<c_size_t> & result) {
//     if(var < 0) {
//         result.push_back(abs(var));
//         return;
//     }

//     const vector<c_size_t> & rhs = slg_nonterm_vec[var].rhs;

//     for(c_size_t rhs_var : rhs) {
//         expandSLG(rhs_var, slg_nonterm_vec, result);
//     }
// }

// vector<c_size_t> expandSLG(const SLG *slg) {
//     vector<c_size_t> result;

//     const vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

//     expandSLG(slg_nonterm_vec.size() - 1, slg_nonterm_vec, result);

//     return result;
// }
// #endif

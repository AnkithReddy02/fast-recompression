#include <bits/stdc++.h>
#include "../include/test_queries.hpp"
#include "../../include/recompression_definitions.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    if(argc != 3) {
        cerr << "Usage: ./random_lce_test input_rlslp input_slp" << endl;
        return 1;
    }

    const string& input_rlslp = string(argv[1]);
    const string& input_slp = string(argv[2]);

    RecompressionRLSLP* rlslp = new RecompressionRLSLP();
    rlslp->read_from_file(input_rlslp);

    if (!test_queries(rlslp, input_slp)) {
        cerr << "LCE Tests Failed!" << endl;
        return 1;
    }
    else {
        cerr << "LCE Tests Passed!" << endl;
    }

    delete rlslp;

    return 0;
}
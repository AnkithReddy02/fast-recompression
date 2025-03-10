#include <bits/stdc++.h>

#include "../include/recompression_definitions.hpp"
#include "../include/typedefs.hpp"

#include "../include/utils/utilities.hpp"
#include "../include/utils/hash_pair.hpp"
#include "../include/utils/packed_pair.hpp"
#include "../include/utils/hash_table.hpp"
#include "../include/utils/utils.hpp"

#include "../include/BComp.hpp"
#include "../include/PComp.hpp"

// #define TEST_SAMPLE
// #define DEBUG_LOG

using namespace std;
using namespace utils;

template<>
std::uint64_t get_hash(const packed_pair<c_size_t, c_size_t> &x) {
    return (std::uint64_t)x.first * (std::uint64_t)4972548694736365 +
           (std::uint64_t)x.second * (std::uint64_t)3878547385475748;
}

template<>
std::uint64_t get_hash(const c_size_t &x) {
    return (std::uint64_t)x * (std::uint64_t)1232545666776865;
}

template<>
std::string keyToString(const c_size_t &x) {
  return utils::intToStr(x);
}

template<>
std::string keyToString(const packed_pair<c_size_t,c_size_t> &x) {
  return "std::pair(" + utils::intToStr(x.first) + ", " + utils::intToStr(x.second) + ")";
}

template<>
std::string valToString(const c_size_t &x) {
  return utils::intToStr(x);
}

template<>
std::string valToString(const bool &x) {
  return utils::intToStr((std::uint64_t)x);
}

uint64_t verbosity;
uint64_t current_run;

RecompressionRLSLP* recompression_on_slp(InputSLP* s) {

    #ifdef TEST_SAMPLE
    delete s;

    s = new InputSLP();
    // s->nonterm.push_back(SLPNonterm('0', 97, 0));
    s->nonterm.push_back(SLPNonterm('0', 98, 0));
    s->nonterm.push_back(SLPNonterm('0', 99, 0));
    s->nonterm.push_back(SLPNonterm('1', 1, 0));
    s->nonterm.push_back(SLPNonterm('1', 2, 0));
    s->nonterm.push_back(SLPNonterm('1', 3, 1));
    s->nonterm.push_back(SLPNonterm('1', 4, 3));
    s->nonterm.push_back(SLPNonterm('1', 5, 4));
    #endif


    SLG* slg = new SLG();

    RecompressionRLSLP* recompression_rlslp = new RecompressionRLSLP();

    // For 0.
    recompression_rlslp->nonterm.push_back(RLSLPNonterm());

    // Compute S0 from S and Initialize Recompression
    for(len_t i=0; i<s->nonterm.size(); ++i) {
        const char_t &type = s->nonterm[i].type;
        const c_size_t &first = s->nonterm[i].first;
        const c_size_t &second = s->nonterm[i].second;

        c_size_t global_rhs_size = slg->rhs.size();

        if(type == '0') {
            // For type '0', in SLG, terminals start from -1.
            // Initially, -1 = recompression_rlslp->nonterm size.
            slg->rhs.push_back(-(c_size_t)(recompression_rlslp->nonterm).size());
            // Here first is the ASCII values of the character.
            recompression_rlslp->nonterm.push_back(RLSLPNonterm('0', first, second));
        }
        else {
            slg->rhs.push_back(first);
            slg->rhs.push_back(second);
        }

        slg->nonterm.push_back(global_rhs_size);
    }

    delete s;

    len_t i = 0;

    double max_BComp_time = 0;
    double max_PComp_time = 0;

    while(++i) {
        // map<pair<c_size_t, c_size_t>, c_size_t> m;
        // hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> m;

        typedef packed_pair<c_size_t, c_size_t> pair_type;

        const space_efficient_vector<SLGNonterm> &slg_nonterm_vec = slg->nonterm;
        const space_efficient_vector<c_size_t> &global_rhs = slg->rhs;
        const SLGNonterm &slg_nonterm = slg_nonterm_vec.back();

        c_size_t start_var = slg_nonterm_vec.size()-1;
        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (start_var == slg_nonterm_vec.size()-1) ? c_size_t(global_rhs.size())-1 : c_size_t(slg_nonterm_vec[start_var+1].start_index) - 1;

        c_size_t start_var_rhs_size = end_index - start_index + 1;

        if(start_var_rhs_size == 1 && global_rhs[end_index] < 0) {
            break;
        }

        ++current_run;

        if(i&1) {
            #ifdef DEBUG_LOG
            cout << '\n' << i << " : Starting BComp..." << endl;
            #endif

            if(verbosity >= 1) {
                cout << "\n\nRound " << i << ": Block Compression" << endl;
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            slg = BComp(slg, recompression_rlslp);

            auto end_time = std::chrono::high_resolution_clock::now();

            auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            max_BComp_time = max(max_BComp_time, duration_seconds);

            if(verbosity >= 1) {
                cout << "  Time = " << duration_seconds << 's' << endl;
            }
        }
        else {
            #ifdef DEBUG_LOG
            cout << '\n' << i << " : Starting PComp..." << endl;
            #endif
            if(verbosity >= 1) {
                cout << "\n\n" << "Round " << i << ": Pair Compression" << endl;
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            slg = PComp(slg, recompression_rlslp);

            auto end_time = std::chrono::high_resolution_clock::now();

            auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            max_PComp_time = max(max_PComp_time, duration_seconds);

            if(verbosity >= 1) {
                cout << "  Time = " << duration_seconds << 's' << endl;
            }
        }

        if(verbosity >= 2) {
            cout << "\nStats:" << endl;
            cout << "  Number of nonterms in recompression RLSLP = " << recompression_rlslp->nonterm.size() << endl;
            cout << "  Number of nonterms in SLG = " << slg->nonterm.size() << endl;
            cout << fixed << setprecision(3) <<  "  Avg expansion size in SLG = " << (long double)slg->rhs.size() / slg->nonterm.size() << endl;

            cout << "\nRAM usage:" << endl;
            cout << fixed << setprecision(3) << "  Recompression RLSLP: " << recompression_rlslp->ram_use() / (1024.0 * 1024.0) << "MiB" << endl;
            cout << fixed << setprecision(3) << "  Current SLG: " << slg->ram_use() / (1024.0 * 1024.0) << "MiB" << endl;
            cout << fixed << setprecision(3) << "  Peak RAM usage so far: " << peak_ram_allocation / (1024.0 * 1024.0) << "MiB" << endl;
            cout << "\n\n";
        }

        // m.clear();
        // m.reset();
        // printRecompressionRLSLP(recompression_rlslp);
    }

    cout << endl;
    cout << "Runs: " << i-1 << endl << endl;
    cout << "Max. BComp Time: " << max_BComp_time << " seconds" << endl;
    cout << "Max. PComp TIme: " << max_PComp_time << " seconds" << endl << endl;

    delete slg;

    return recompression_rlslp;
}

void start_compression(const string &input_file, const string &output_file = "", const string &input_text = "") {


    InputSLP *inputSLP = new InputSLP();
    inputSLP->read_from_file(input_file);
    // InputSLP *inputSLP = getSLP(50);

    c_size_t j = 0;
    c_size_t i = -1;

    auto start_time = std::chrono::high_resolution_clock::now();

    RecompressionRLSLP *recompression_rlslp = recompression_on_slp(inputSLP);

    recompression_rlslp->computeExplen();

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Output the duration
    cout << "Time taken for Construction: " << duration_seconds << " seconds" << endl;
    cout << "Total Time taken: " << duration_seconds << " seconds" << endl;
    space_efficient_vector<RLSLPNonterm> & rlslp_nonterm_vec = recompression_rlslp->nonterm;
    // // printRecompressionRLSLP(recompression_rlslp);
    c_size_t text_size = rlslp_nonterm_vec.back().explen;

    cout << "Text Size : " << text_size << endl << endl;

    cout << "****************" << endl;
    cout << "Peak RAM usage for Construction: " << get_peak_ram_allocation()/(1024.0 * 1024.0) << " MB" << endl;
    cout << "****************" << endl << endl;

    // cout << "****************" << endl;
    // cout << "Overall Peak RAM uage: " << get_peak_ram_allocation() << endl;
    // cout << "****************" << endl << endl;

    // Write Recompression to the 'output_file'
    if(!output_file.empty()) {
        cout << "Writing to file: " << output_file << endl;
        recompression_rlslp->write_to_file(output_file);
        cout << "Done writing to file!" << endl;
    }

    // RLSLP expansion test.
    {
        if(!input_text.empty()) {
            cout << "Testing RLSLP expansion..." << endl;
	    space_efficient_vector<c_size_t> expansion;
            expandRLSLP(recompression_rlslp, expansion);

            string input_text_name = utils::absolute_path(input_text);
            c_size_t text_length = utils::file_size(input_text_name);
            cout << "text length: " << text_length << endl;
            char_t *text = allocate_array<char_t>(text_length);
            utils::read_from_file(text, text_length, input_text_name);

            if(expansion.size() != text_length) {
                cerr << "Text and Expansion length didn't match!!" << endl;
            }
	    else {
                bool correct_expansion = true;
                for(uint64_t i = 0; i < text_length; ++i) {
                    if(text[i] != expansion[i]) {
                        cerr << "Mismtach -- RLSLP Expansion and text at Index: " << i << endl;
                        correct_expansion = false;
			break;
                    }
                }

	        if(correct_expansion) {
	            cout << "Correct Extract!" << endl;
                }
	    }

            deallocate(text);
        }
    }

    // Clean up.
    delete recompression_rlslp;

    cout << endl << "****************" << endl;
    cout << "Current RAM usage for Construction: " <<
      get_current_ram_allocation() << " bytes" << endl;
    cout << "****************" << endl << endl;
}

int main(int argc, char *argv[]) {
    utils::initialize_stats();

    if(argc <= 1) {
        cout << "Usage: " + string(argv[0]) + " [-v or -vv] input_slp -o output_file -t input_text" << endl;
        exit(1);
    }

    string input_slp;
    string output_file;
    string input_text;

    verbosity = 0;

    for(len_t i = 1; i < argc; ++i) {
        string curr_arg(argv[i]);

        if(curr_arg == "-o") {
            output_file = string(argv[++i]);
        }
        else if(curr_arg == "-v") {
            verbosity = 1;
        }
        else if(curr_arg == "-vv") {
            verbosity = 2;
        }
        else if(curr_arg == "-t") {
            input_text = string(argv[++i]);
        }
        else {
            input_slp = string(argv[i]);
        }
    }

    if(input_slp.empty()) {
        cerr << "Error: Input SLP file is empty!" << endl;
        exit(1);
    }

    print_current_timestamp();

    cout << "SLP File: " << input_slp << endl;
    cout << "Output File: " << output_file << endl;

    start_compression(input_slp, output_file, input_text);
    return 0;
}

#include<bits/stdc++.h>
#include "../include/utils/bm_compression.hpp"
// #define DEBUG 

int main(int argc, char* argv[]) {
    /*
    
        possible examples;
            
        1. "acxybcbcxy" block_size = 2
            when to extend backward?
            Consider, xy starting from position 8.
            Extend backward only if there is a alone character?
        2. zxQPYxQYPzab block_size = 2
    */

    if(argc < 4) {
        cerr << "Usage: " << argv[0] << " text_file_name -o output_file_name" << endl;
        return 1;
    }

    string text_file_name = argv[1];
    string output_file_name = argv[3];

    cout << "Input text: " << text_file_name << endl;
    cout << "Output: " << output_file_name << endl;

    uint64_t block_size = 50;

    cout << "Block Size: " << block_size << endl;
    
    utils::initialize_stats();

    {
        BM bm(text_file_name, block_size, output_file_name);

        auto compress_start = std::chrono::high_resolution_clock::now();
        // cout << "Computing Optimal Block Size..." << endl;
        // cout << "Optimal Block Size: " << bm.set_opt_block_size() << endl;
        auto compress_end = std::chrono::high_resolution_clock::now();
        cout << "Compressing text using BM..." << endl;
        bm.compress();
        cout << "Compression complete." << endl;
        bm.write_to_file();
        std::chrono::duration<double> elapsed = compress_end - compress_start;
        cout << "time = " << elapsed.count() << " seconds" << endl;
        bm.printStats();
        // bm.test();
    }
    cout << "peak = " << get_peak_ram_allocation() / (1024.0 * 1024.0) << "MiB" << endl;
    cout << "current = " << get_current_ram_allocation() / (1024.0 * 1024.0) << "MiB" << endl;
    return 0;
}

#include<bits/stdc++.h>
#include "../include/utils/hash_table.hpp"
#include "../include/utils/space_efficient_vector.hpp"
#include "../include/utils/packed_pair.hpp"

// #define DEBUG 

using namespace std;

template<>
std::uint64_t get_hash(const packed_pair<uint64_t, uint64_t> &x) {
    return (std::uint64_t)x.first * (std::uint64_t)4972548694736365 +
           (std::uint64_t)x.second * (std::uint64_t)3878547385475748;
}

template<>
std::uint64_t get_hash(const uint64_t &x) {
    return (std::uint64_t)x * (std::uint64_t)1232545666776865;
}


struct Block {
   uint64_t length;
   uint64_t position;

   Block(uint64_t position, uint64_t length) : length(length), position(position) {
    
   }

   Block() {
    length = 0;
    position = 0;
   }
};

class BM {
private:
    hash_table<uint64_t, Block> hash_block;
   
    const char* text;
    uint64_t block_size;
    uint64_t text_length;

public:
    unsigned int stringToHash(string& s) {
        // Lot of Collsions.
        // TODO: Integrate Karp-Rabin
        return 0;
    }


    BM(const char* text, uint64_t block_size): text(text), block_size(block_size) {
        text_length = strlen(text);
    }

    space_efficient_vector<pair<uint64_t, uint64_t>> compress() {
        
        uint64_t current_hash = 0;
        uint64_t position = 0;

        uint64_t buffer_capacity = 1024;
        char* buffer = new char[buffer_capacity];
        uint64_t buffer_length = 0;

        space_efficient_vector<pair<uint64_t, uint64_t>> parsing;

        while(position < text_length) {
            if(text_length - position < block_size) {
                break;
            }

            // PENDING - INTEGRATE KARP-RABIN
            if(current_hash == 0) {
                // THIS IS SAMPLE CODE.
                string curr;
                for(int i = position; i < position + block_size; i++) {
                    curr.push_back(text[i]);
                }

                current_hash = stringToHash(curr);
                // Compute hash from [text + position, text + position + block_size)
            }
            // PENDING - INTEGRATE KARP-RABIN
            else {
                string curr;
                // compute next_hash;
                for(int i = position; i < position + block_size; i++) {
                    curr.push_back(text[i]);
                }

                current_hash = stringToHash(curr);
            }

            if(position % block_size == 0) {
                if(hash_block.find(current_hash) == NULL) {
                    const Block new_block = Block(position, block_size);
                    hash_block.insert(current_hash, new_block);
                }
            }

            Block match;
            uint64_t begin_len = 0;
            uint64_t end_len = 0;
            {
                
                const Block* block = hash_block.find(current_hash);
                // const Block* block = hash_block[current_hash];
                if(block != NULL && block->position != position) {
                    // Note the difference
                    uint64_t block_end_position = block->position + block->length;
                    uint64_t current_end_position = position + block_size;

                    uint64_t block_begin_position = block->position - 1;
                    uint64_t current_begin_position = position - 1;

                    // Forward Match.
                    uint64_t end_index = 0;

                    while(block_end_position + end_index < text_length && current_end_position + end_index < text_length) {
                        if(text[block_end_position + end_index] != text[current_end_position + end_index]) {
                            break;
                        }

                        end_index++;
                    }

                    

                    // Moved end_index steps forward.

                    // Backward Match -- Max block_size - 1 characters

                    uint64_t begin_index = 0;
                    // uint64_t parsingIndex = parsing.size(); 

                    while(block_begin_position - begin_index >= 0 && current_begin_position - begin_index >= 0 && begin_index < block_size - 1
                    && parsing.size() > 0 && parsing.back().second == 0) {
                        if(text[block_begin_position - begin_index] != text[current_begin_position - begin_index]) {
                            break;
                        }
                        parsing.pop_back();
                        begin_index++;
                    }

                    begin_len = begin_index;
                    end_len = end_index;

                    match.position = block->position - begin_index;
                    match.length = block->length + begin_index + end_index;
                }
            }

            if(match.length > 0) {
                parsing.push_back(make_pair(match.position, match.length));
                // cout << "Line 151" <<  position << ' ' << match.position << ' ' << match.length << endl;
                current_hash = 0;
                position = position - begin_len + match.length;
                // cout << "Line 154: " << match.position << ' ' << match.length << endl;
            }
            else {
                parsing.push_back(make_pair(text[position], 0));
                // cout << position << ' ' << position << ' ' << 0 << endl;
                ++position;
            }

            // cout << position << endl;
        }

        // Remaining text.
        while(position < text_length) {
            parsing.push_back(make_pair(text[position], 0));
            ++position;
        }

        #ifdef DEBUG
        for(uint64_t i = 0; i < parsing.size(); ++i) {
            cout << parsing[i].first << ' ' << parsing[i].second << endl;
        }
        #endif

        return parsing;
        
    }

    string decompress(space_efficient_vector<pair<uint64_t, uint64_t>>& parsing) {
        string result;

        for(uint64_t i = 0; i < parsing.size(); ++i) {
            if(parsing[i].second == 0) {
                result.push_back(parsing[i].first);
            }
            else {
                uint64_t start_index = parsing[i].first;
                uint64_t end_index = parsing[i].first + parsing[i].second - 1;

                for(uint64_t j = start_index; j <= end_index; ++j) {
                    result.push_back(result[j]);
                }
            }
        }
        return result;
    }
};

int main() {
    /*
    
        possible examples;
            
        1. "acxybcbcxy" block_size = 2
            when to extend backward?
            Consider, xy starting from position 8.
            Extend backward only if there is a alone character?
        2. zxQPYxQYPzab block_size = 2
    */
    const char* sample_text = "acxybcbcxy";
    uint64_t block_size = 2;
    BM bm(sample_text, block_size);

    space_efficient_vector<pair<uint64_t, uint64_t>> parsing = bm.compress();
    cout << "Compression complete." << endl;
    string decompressed_string = bm.decompress(parsing);
    cout << "Decompression complete." << endl;
    cout << "Original String: " << string(sample_text) << endl;
    cout << "Decompressed String: " << decompressed_string << endl;
    cout << ((string(sample_text) == decompressed_string) ? "Equal" : "Mismatch!!") << endl;
    cout << "Parsing Size: " << parsing.size() << endl;
    cout << "String length: " << strlen(sample_text) << endl;


    return 0;
}

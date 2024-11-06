#ifndef BM_COMPRESSION_HPP
#define BM_COMPRESSION_HPP

#include "hash_table.hpp"
#include "space_efficient_vector.hpp"
#include "packed_pair.hpp"
#include "utils.hpp"
#include "../typedefs.hpp"
#include "karp_rabin_hashing.hpp"

using namespace std;
using namespace karp_rabin_hashing;
using namespace utils;

template<>
std::uint64_t get_hash(const packed_pair<uint64_t, uint64_t> &x) {
    return (std::uint64_t)x.first * (std::uint64_t)4972548694736365 +
           (std::uint64_t)x.second * (std::uint64_t)3878547385475748;
}

template<>
std::uint64_t get_hash(const uint64_t &x) {
    return (std::uint64_t)x * (std::uint64_t)1232545666776865;
}

struct __attribute__((packed)) Block {
   // redudandant.
   uint64_t length;
   uint64_t position;
   uint64_t next;

   Block(uint64_t position, uint64_t length) : length(length), position(position) {
    next = -1;
   }

   Block() {
    length = 0;
    position = 0;
   }
};
// 148067 * 16 + 148160 * 24 + 429265758 + 10000 + 118092513 * 16
struct ListPtr {
    // Index of first element.
    uint64_t head;
    // Index of last element.
    uint64_t tail;

    ListPtr(uint64_t head, uint64_t tail) : head(head), tail(tail) {

    }

    ListPtr() {
        head = -1;
        tail = -1;
    }
};

class BM {
private:
    hash_table<uint64_t, ListPtr> hash_block;
    space_efficient_vector<Block> blocks;
    space_efficient_vector<pair<c_size_t, c_size_t>> parsing;
    string output_file_name;
   
    char_t* text;
    uint64_t block_size;
    uint64_t text_length;

public:
    BM(char_t* text, uint64_t block_size, const string& output_file_name): text(text), block_size(block_size), output_file_name(output_file_name) {
        text_length = strlen(reinterpret_cast<const char*>(text));
        karp_rabin_hashing::init();
    }

    BM(string text_file_name, uint64_t block_size, const string& output_file_name): output_file_name(output_file_name){
        text_file_name = utils::absolute_path(text_file_name); 
        text_length = utils::file_size(text_file_name);
        cout << "text length: " << text_length << endl;
        text = allocate_array<char_t>(text_length);
        utils::read_from_file(text, text_length, text_file_name);
        karp_rabin_hashing::init();
        this->block_size = min(block_size, text_length/2);
    }

    ~BM() {
        deallocate(text);
    }

    void compress() {
        parsing.clear();
        uint64_t current_hash = 0;
        uint64_t position = 0;

        

        while(position < text_length) {
            if(text_length - position < block_size) {
                break;
            }

            // PENDING - INTEGRATE KARP-RABIN
            if(current_hash == 0) {
                // THIS IS SAMPLE CODE.
                // string curr;
                // for(int i = position; i < position + block_size; i++) {
                //     curr.push_back(text[i]);
                // }

                // current_hash = stringToHash(curr);

                char_t* curr_str = allocate_array<char_t>(block_size);
                for(int i = position; i < position + block_size; i++) {
                    curr_str[i - position] = (text[i]);
                }
                current_hash = hash_string<char_t>(curr_str, block_size);

                deallocate(curr_str);
                // Compute hash from [text + position, text + position + block_size)
            }
            // PENDING - INTEGRATE KARP-RABIN
            else {
                // string curr;
                // // compute next_hash;
                // for(int i = position; i < position + block_size; i++) {
                //     curr.push_back(text[i]);
                // }

                // current_hash = stringToHash(curr);

                uint64_t left_hash = hash_char<char_t>(text[position - 1]);
                uint64_t combined_hash = current_hash;
                uint64_t remaining_left_hash = unconcat(left_hash, combined_hash, block_size - 1);

                uint64_t right_hash = hash_char<char_t>(text[position + block_size - 1]);

                current_hash = concat(remaining_left_hash, right_hash, 1);
            }

            if(position % block_size == 0) {
                const Block new_block = Block(position, block_size);
                if(hash_block.find(current_hash) == NULL) {
                    const ListPtr list_ptr = ListPtr(blocks.size(), blocks.size());
                    blocks.push_back(new_block);
                    hash_block.insert(current_hash, list_ptr);
                }
                else {
                    ListPtr* list_ptr = hash_block.find(current_hash);
                    blocks[list_ptr->tail].next = blocks.size();
                    list_ptr->tail = blocks.size();
                    blocks.push_back(new_block);
                }
            }

            Block match;
            uint64_t begin_len = 0;
            uint64_t end_len = 0;
            {
                const ListPtr* list_ptr = hash_block.find(current_hash);

                if(list_ptr != NULL) {
                    uint64_t list_start_index = list_ptr->head;
                    while(list_start_index != -1) {
                        const Block* block = &blocks[list_start_index];
                        if(block->position >= position) {
                            break;
                        }

                        bool stringsEqual = true;

                        uint64_t block_start = block->position;
                        uint64_t curr_start = position;
                        uint64_t len = 0;
                        while(len < block_size) {
                            if(text[block_start + len] != text[curr_start + len]) {
                                stringsEqual = false;
                                break;
                            }
                            len++;
                        }

                        if(!stringsEqual) {
                            list_start_index = block->next;
                            continue;
                        }

                        // const Block* block = hash_block[current_hash];
                        if(block != NULL && block->position != position) {
                            // Note the difference
                            uint64_t block_end_position = (block->position )+ (block->length);
                            uint64_t current_end_position = position + block_size;

                            uint64_t block_begin_position = (block->position) - 1;
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
                            && parsing.size() > 0 && parsing.back().second == (uint40)0) {
                                if(text[block_begin_position - begin_index] != text[current_begin_position - begin_index]) {
                                    break;
                                }
                                parsing.pop_back();
                                ++begin_index;
                            }

                            begin_len = begin_index;
                            end_len = end_index;

                            match.position = block->position - begin_index;
                            match.length = block->length + begin_index + end_index;
                            // ** Missed this **
                            break;
                        }

                        list_start_index = block->next;
                    }
                }

                
            }

            if(match.length > 0) {
                c_size_t first = match.position;
                c_size_t second = match.length;
                parsing.push_back(make_pair(first, second));
                current_hash = 0;
                position = position - begin_len + match.length;
            }
            else {
              
                parsing.push_back(make_pair(text[position], 0));
                ++position;
            }
        }

        // Remaining text.
        while(position < text_length) {
            parsing.push_back(make_pair(text[position], 0));
            ++position;
        }

        #ifdef DEBUG
        cout << endl << endl;
        for(uint64_t i = 0; i < parsing.size(); ++i) {
            cout << parsing[i].first << ' ' << parsing[i].second << endl;
        }
        #endif
        
        return;
    }

    void write_to_file() {
        cout << "Writing to file... " << endl;
        parsing.write_to_file(output_file_name);
        return;
    }

    string decompress(space_efficient_vector<pair<c_size_t, c_size_t>>& parsing) {
        string result;

        for(uint64_t i = 0; i < parsing.size(); ++i) {
            if(parsing[i].second == (uint40)0) {
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
        cout << "Decompression complete." << endl;
        return result;
    }

    void printStats() {
        cout << "Hash table size: " << hash_block.size() << endl;
        cout << "Parsing Size: " << parsing.size() << endl;
        cout << "String length: " << strlen(reinterpret_cast<const char*>(text)) << endl;
        cout << "Num Blocks: " << blocks.size() << endl;
        cout << "Reduction: " << 100 - (double)(parsing.size() * 16) * 100 / (strlen(reinterpret_cast<const char*>(text))) << '%' << endl;

        for(int i = 0; i < parsing.size(); ++i) {
            // cout << parsing[i].first << ' ' << parsing[i].second << endl;
        }
    }

    uint64_t set_opt_block_size() {
        uint64_t low = 1;
        uint64_t high = min((uint64_t)2000, text_length);

        while(high - low > 2) {
            uint64_t block_size_1 = low + (high - low) / 3;
            uint64_t block_size_2 = high - (high - low) / 3;

            block_size = block_size_1;
            compress();
            uint64_t parsing_size_1 = parsing.size();
            block_size = block_size_2;
            compress();
            uint64_t parsing_size_2 = parsing.size();

            if(parsing_size_1 < parsing_size_2) {
                high = block_size_2;
            }
            else {
                low = block_size_1;
            }
        }

        block_size = low + 1;
        return min(low + 1, text_length);
    }

    void test() {
        compress();
        
        string decompressed_string = decompress(parsing);
        // cout << "Original String: " << string(reinterpret_cast<const char*>(text)) << endl;
        // cout << "Decompressed String: " << decompressed_string << endl;
        cout << ((string(reinterpret_cast<const char*>(text)) == decompressed_string) ? "Equal" : "Mismatch!!") << endl;
        cout << "Parsing Size: " << parsing.size() << endl;
        cout << "String length: " << strlen(reinterpret_cast<const char*>(text)) << endl;
        cout << "Num Blocks: " << blocks.size() << endl;
        cout << "Reduction: " << 100 - (double)(parsing.size()) * 100 / (strlen(reinterpret_cast<const char*>(text))) << '%' << endl;
    }
};

#endif

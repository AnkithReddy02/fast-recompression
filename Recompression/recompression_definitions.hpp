#ifndef RECOMPRESSION_DEFINITIONS_HPP
#define RECOMPRESSION_DEFINITIONS_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include "typedefs.hpp"

using namespace std;

struct __attribute__((packed)) RLSLPNonterm {
    char_t type;
    c_size_t first;
    c_size_t second;
    c_size_t explen;

    RLSLPNonterm(const char_t &type, const c_size_t &first, const c_size_t &second) : type(type), first(first), second(second), explen(0) {}
    RLSLPNonterm() : type('0'), first(0), second(0), explen(0) {
    }
    RLSLPNonterm(const char_t &type, const c_size_t &first, const c_size_t &second, const c_size_t &explen) : type(type), first(first), second(second), explen(explen) {}
};

class RecompressionRLSLP {
public:
    vector<RLSLPNonterm> nonterm;
};

struct SLGNonterm {
    // rhs empty represents empty variable.
    vector<c_size_t> rhs;
    //int LMS;
    //int RMS;

    // 0 is empty
    // int LB;
    // int RB;
    
    // {-1, -1} is empty
    // pair<int, int> LR;
    // pair<int, int> RR;

    SLGNonterm(const vector<c_size_t> &rhs) : rhs(rhs) {

    }

    SLGNonterm() {

    }
};

class SLG {

public:
    SLG() {

    }
    SLG(vector<SLGNonterm> & nonterm) : nonterm(nonterm) {

    }
    vector<SLGNonterm> nonterm;
};

struct  __attribute__((packed)) SLPNonterm {
    char_t type;
    c_size_t first;
    c_size_t second;

    SLPNonterm(const char_t &type, const c_size_t &first, const c_size_t &second) : type(type), first(first), second(second) {

    }

    SLPNonterm() : type('0'), first(0), second(0) {

    }
};

class InputSLP {
public:
    vector<SLPNonterm> nonterm;

    InputSLP() {

    }

    InputSLP(const vector<SLPNonterm>& nonterm) : nonterm(nonterm) {

    }

    void read_from_file(const string &file_name) {
        ifstream file(file_name, ios::binary);

        if (!file.is_open()) {
            cout << "Error: Unable to open the file!" << endl;
            return;
        }

        file.seekg(0, ios::end);
        size_t file_size = file.tellg();
        file.seekg(0, ios::beg);
        cout << "File size: " << file_size << " bytes" << endl;

        unsigned char buffer[10]; // Buffer to hold 10 bytes (5 bytes read two times)

        nonterm.resize(file_size/10);

        c_size_t i = 0;

        while (file.read(reinterpret_cast<char*>(buffer), sizeof(buffer))) {
            // Process the first 5-byte integer
            uint64_t value1 = 0;
            for (c_size_t i = 0; i < 5; ++i) { // First 5 bytes
                value1 |= static_cast<uint64_t>(buffer[i]) << (i * 8);
            }

            // Process the second 5-byte integer
            uint64_t value2 = 0;
            for (c_size_t i = 5; i < 10; ++i) { // Next 5 bytes
                value2 |= static_cast<uint64_t>(buffer[i]) << ((i - 5) * 8);
            }

            if (value1 == 0) {
                nonterm[i++] = SLPNonterm('0', value2, 0);
            } else {
                nonterm[i++] = SLPNonterm('1', value1 - 1, value2 - 1);
            }
        }

        // Check if we have read less than 5 bytes in the last chunk
        c_size_t lastChunkSize = file.gcount();
        if (lastChunkSize > 0) {
            cout << "Read " << lastChunkSize << " bytes in the last chunk, incomplete for a 64-bit value." << endl;
        }

        if (file.eof()) {
            cout << "Reached end of the file!" << endl;
        } else {
            cout << "Error in reading file!" << endl;
        }

        file.close();

        order_slp();
    }

private:
    /*
        Space : O(|G|)
        Time  : O(|G|)
    */

    void order_slp() {

        const c_size_t grammar_size = nonterm.size();
        vector<vector<c_size_t>> graph(grammar_size, vector<c_size_t>());
        vector<uint8_t> inorder(grammar_size, 0);
        vector<c_size_t> old_new_map(grammar_size, 0);

        queue<c_size_t> q;

        // Reverse Graph!
        for(c_size_t i = 0; i < nonterm.size(); i++) {
            if(nonterm[i].type == '1') {
                graph[nonterm[i].first].push_back(i);
                graph[nonterm[i].second].push_back(i);

                inorder[i] += 2;
            }
            else {
                q.push(i);
            }
        }

        c_size_t nonterminal_ptr = 0;

        while(!q.empty()) {
            c_size_t u = q.front();
            q.pop();

            old_new_map[u] = nonterminal_ptr++;

            for(c_size_t v : graph[u]) {
                inorder[v]--;

                if(inorder[v] == 0) {
                    q.push(v);
                }
            }
        }

        graph.clear();
        inorder.clear();

        vector<SLPNonterm> ordered_nonterm(grammar_size);

        for(c_size_t i=0; i<nonterm.size(); i++) {

            const char_t &type = nonterm[i].type;
            const c_size_t &first = nonterm[i].first;
            const c_size_t &second = nonterm[i].second;

            if(type == '0') {
                ordered_nonterm[old_new_map[i]] = SLPNonterm('0', first, second);
            }
            else {
                ordered_nonterm[old_new_map[i]] = SLPNonterm('1', old_new_map[first], old_new_map[second]);
            }
        }
        
        nonterm = move(ordered_nonterm);

        return;
    }
};

struct Node {
    // Variable/Non-Terminal
    c_size_t var;
    // [l, r)
    c_size_t l;
    c_size_t r;

    Node() {
        
    }
    Node(const c_size_t &var, const c_size_t &l, const c_size_t &r) : var(var), l(l), r(r) {

    }
};

struct __attribute__((packed)) AdjListElement {
    c_size_t first;
    c_size_t second;
    bool_t swapped;
    c_size_t vOcc;

    AdjListElement() {}

    AdjListElement(const c_size_t &first, const c_size_t &second, const bool_t &swapped, const c_size_t &vOcc) : first(first), second(second), swapped(swapped), vOcc(vOcc) {}
};

#endif

#ifndef RECOMPRESSION_DEFINITIONS_HPP
#define RECOMPRESSION_DEFINITIONS_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <cstdint>

using namespace std;

struct RLSLPNonterm {
    char type;
    int first;
    int second;
    int explen;

    RLSLPNonterm(const char &type, const int &first, const int &second) : type(type), first(first), second(second), explen(0) {}
    RLSLPNonterm() : type('0'), first(0), second(0), explen(0) {
    }
    RLSLPNonterm(const char &type, const int &first, const int &second, const int &explen) : type(type), first(first), second(second), explen(explen) {}
};

class RecompressionRLSLP {
public:
    vector<RLSLPNonterm> nonterm;
};

struct SLGNonterm {
    // rhs empty represents empty variable.
    vector<int> rhs;
    int vOcc;
    //int LMS;
    //int RMS;

    // 0 is empty
    int LB;
    int RB;
    
    // {-1, -1} is empty
    pair<int, int> LR;
    pair<int, int> RR;

    SLGNonterm(const vector<int> &rhs) : vOcc(0), LB(0), RB(0), LR({-1, -1}), RR({-1, -1}), rhs(rhs) {

    }

    SLGNonterm() : vOcc(0), LB(0), RB(0), LR({-1, -1}), RR({-1, -1}) {

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

struct SLPNonterm {
    char type;
    int first;
    int second;

    SLPNonterm(const char &type, const int &first, const int &second) : type(type), first(first), second(second) {

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

        int i = 0;

        while (file.read(reinterpret_cast<char*>(buffer), sizeof(buffer))) {
            // Process the first 5-byte integer
            uint64_t value1 = 0;
            for (int i = 0; i < 5; ++i) { // First 5 bytes
                value1 |= static_cast<uint64_t>(buffer[i]) << (i * 8);
            }

            // Process the second 5-byte integer
            uint64_t value2 = 0;
            for (int i = 5; i < 10; ++i) { // Next 5 bytes
                value2 |= static_cast<uint64_t>(buffer[i]) << ((i - 5) * 8);
            }

            if (value1 == 0) {
                nonterm[i++] = SLPNonterm('0', value2, 0);
            } else {
                nonterm[i++] = SLPNonterm('1', value1 - 1, value2 - 1);
            }
        }

        // Check if we have read less than 5 bytes in the last chunk
        int lastChunkSize = file.gcount();
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

        const int grammar_size = nonterm.size();
        vector<vector<int>> graph(grammar_size, vector<int>());
        vector<uint8_t> inorder(grammar_size, 0);
        vector<int> old_new_map(grammar_size, 0);

        queue<int> q;

        // Reverse Graph!
        for(int i = 0; i < nonterm.size(); i++) {
            if(nonterm[i].type == '1') {
                graph[nonterm[i].first].push_back(i);
                graph[nonterm[i].second].push_back(i);

                inorder[i] += 2;
            }
            else {
                q.push(i);
            }
        }

        int nonterminal_ptr = 0;

        while(!q.empty()) {
            int u = q.front();
            q.pop();

            old_new_map[u] = nonterminal_ptr++;

            for(int v : graph[u]) {
                inorder[v]--;

                if(inorder[v] == 0) {
                    q.push(v);
                }
            }
        }

        graph.clear();
        inorder.clear();

        vector<SLPNonterm> ordered_nonterm(grammar_size);

        for(int i=0; i<nonterm.size(); i++) {

            const char &type = nonterm[i].type;
            const int &first = nonterm[i].first;
            const int &second = nonterm[i].second;

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
    int var;
    // [l, r)
    int l;
    int r;

    Node() {
        
    }
    Node(const int &var, const int &l, const int &r) : var(var), l(l), r(r) {

    }
};

#endif

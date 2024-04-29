#ifndef RECOMPRESSION_DEFINITIONS_HPP
#define RECOMPRESSION_DEFINITIONS_HPP

#include<bits/stdc++.h>
using namespace std;

struct RLSLPNonterm {
    char type;
    int first;
    int second;
    int explen;

    RLSLPNonterm(char type, int first, int second) : type(type), first(first), second(second), explen(0) {}
    RLSLPNonterm() : type('0'), first(0), second(0), explen(0) {
    }
    RLSLPNonterm(char type, int first, int second, int explen) : type(type), first(first), second(second), explen(explen) {}
};

class RecompressionRLSLP {
public:
    vector<RLSLPNonterm> nonterm;
};

struct SLGNonterm {
    // rhs empty represents empty variable.
    vector<int> rhs;
    int vOcc;
    int LMS;
    int RMS;

    // 0 is empty
    int LB;
    int RB;
    
    // {-1, -1} is empty
    pair<int, int> LR;
    pair<int, int> RR;

    SLGNonterm(vector<int> rhs) : vOcc(0), LMS(0), RMS(0), LB(0), RB(0), LR({-1, -1}), RR({-1, -1}), rhs(rhs) {

    }

    SLGNonterm() : vOcc(0), LMS(0), RMS(0), LB(0), RB(0), LR({-1, -1}), RR({-1, -1}) {

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

    SLPNonterm(char type, int first, int second) : type(type), first(first), second(second) {

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

    // https://stackoverflow.com/a/2974659
    // void read_from_file(const string & file_name) {
    //     ifstream file(file_name, ios::binary);

    //     if (!file.is_open()) {
    //         cout << "Error: Unable to open the file!" << endl;
    //         return;
    //     }

    //     // Move the file pointer to the end.
    //     file.seekg(0, ios::end);
    //     // Get the file size.
    //     size_t file_size = file.tellg();
    //     // Move back the file pointer to the beginning.
    //     file.seekg(0, ios::beg);

    //     cout << "File size: " << file_size << endl;

    //     // if (file_size % 4) {
    //     //     cout << "Error: File size is not a multiple of 4!" << endl;
    //     //     return;
    //     // }

    //     int first, second;
    //     char c;

    //     while (file.read(reinterpret_cast<char*>(&first), sizeof(int))) {
    //         if (first == 0) {
    //             // Read a 1-byte integer for second.
    //             file.read(reinterpret_cast<char*>(&c), sizeof(char));
    //             nonterm.push_back(SLPNonterm('0', c, 0));

    //             cout << c << ' ' << 0 << endl;
    //         } else {
    //             // Read a 4-byte integer for second.
    //             file.read(reinterpret_cast<char*>(&second), sizeof(int));
    //             nonterm.push_back(SLPNonterm('1', first - 1, second - 1));

    //             cout << first-1 << ' ' << second-1 << endl;
    //         }
    //     }

    //     if (file.eof()) {
    //         cout << "Reached end of the file!" << endl;
    //     } else {
    //         cout << "Error in reading file!" << endl;
    //     }

    //     file.close();

    //     order_slp();
        
    // }


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
            nonterm.push_back(SLPNonterm('0', value2, 0));
        } else {
            nonterm.push_back(SLPNonterm('1', value1 - 1, value2 - 1));
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
    // void order_slp() {
    //     assert(nonterm.size() > 0);

    //     vector<SLPNonterm> ordered_nonterm(nonterm.size());

    //     // Kind of visited array but stores the newly assigned nonterminal. 
    //     int dp[nonterm.size()];
    //     fill(dp, dp + nonterm.size(), -1);

    //     queue<int> q;

    //     // Initialize the queue with the START non-terminal.
    //     q.push(nonterm.size()-1);

    //     // First Assign Terminals.
    //     int nonterminal_begin = 0;

    //     for(int i=0; i<nonterm.size(); i++) {
    //         if(nonterm[i].type == '0') {
    //             ordered_nonterm[nonterminal_begin] = nonterm[i];
    //             dp[i] = nonterminal_begin++;
    //         }
    //     }

    //     // Pointer of ordered nonterm to reassign SLPNonterm from back. Dry run.
    //     int ordered_nonterm_ptr = ordered_nonterm.size()-1;

    //     // Assign.
    //     int nonterminal_end = nonterm.size()-1;

    //     // Assign the start Non-Terminal.
    //     dp[nonterm.size()-1] = nonterminal_end--;

    //     int j = 0;

    //     while(!q.empty()) {
    //         int sz = q.size();

    //         while(sz--) {
    //             int node = q.front();
    //             q.pop();

    //             char type = nonterm[node].type;
    //             int first = nonterm[node].first;
    //             int second = nonterm[node].second;

    //             // Explore Neighbors.
    //             if(type == '1') {
    //                 if(dp[first] == -1) {
    //                     // Push to queue only if it is not explored!
    //                     q.push(first);
    //                     // Assign Neighbor.
    //                     dp[first] = nonterminal_end--;
    //                 }

    //                 if(dp[second] == -1) {
    //                     q.push(second);
    //                     // ASsign Neighbor.
    //                     dp[second] = nonterminal_end--;
    //                 }

    //                 // Add the current Non-Terminal RHS to the end of ordered_nonterm.
    //                 ordered_nonterm[ordered_nonterm_ptr--] = SLPNonterm('1', dp[first], dp[second]);
                    

    //                 j++;

    //                 // cout << dp[first] << ' ' << dp[second] << endl;
    //             }
    //             // 'else' never executes! -- Terminals are explored initially.
    //             else {
    //                 ordered_nonterm[ordered_nonterm_ptr--] = SLPNonterm('1', first, second);
    //             }
    //         }
    //     }

    //     // Reassign the nonterm.
    //     nonterm = ordered_nonterm;
    // }

    void order_slp() {

        int grammar_size = nonterm.size();
        vector<vector<int>> graph(grammar_size, vector<int>());
        vector<int> inorder(grammar_size, 0);
        vector<int> old_new_map(grammar_size, 0);

        queue<int> q;

        // Reverse Graph!
        for(int i=0; i<nonterm.size(); i++) {
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
            int sz = q.size();

            while(sz--) {

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
        }

        vector<SLPNonterm> ordered_nonterm(grammar_size);

        for(int i=0; i<nonterm.size(); i++) {

            char type = nonterm[i].type;
            int first = nonterm[i].first;
            int second = nonterm[i].second;

            if(type == '0') {
                ordered_nonterm[old_new_map[i]] = SLPNonterm('0', first, second);
            }
            else {
                ordered_nonterm[old_new_map[i]] = SLPNonterm('1', old_new_map[first], old_new_map[second]);
            }
        }
        
        nonterm = ordered_nonterm;

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
    Node(int var, int l, int r) : var(var), l(l), r(r) {

    }
};

#endif
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
};

class InputSLP {
public:
    vector<SLPNonterm> nonterm;

    InputSLP() {

    }
    InputSLP(const vector<SLPNonterm>& nonterm) : nonterm(nonterm) {

    }

    // https://stackoverflow.com/a/2974659
    void read_from_file(const string & file_name) {
        ifstream file;

        file.open(file_name, ios::binary);

        if(!file) {
            cout << "Error in opening file!" << endl;
            return;
        }

        int num_terminals;

        file.read(reinterpret_cast<char*>(&num_terminals), sizeof(int));

        char terminal;

        for(int i=1; i<=num_terminals; i++) {
            file.read(reinterpret_cast<char*>(&terminal), sizeof(terminal));
            nonterm.push_back(SLPNonterm('0', -i, 1));
        }

        int non_terminal1, non_terminal2;

        while(true) {
            if(file.read(reinterpret_cast<char*>(&non_terminal1), sizeof(int)) && file.read(reinterpret_cast<char*>(&non_terminal2), sizeof(int))) {
                nonterm.push_back(SLPNonterm('1', non_terminal1, non_terminal2));
            }
            else {
                if(file.eof()) {
                    cout << "Reached end of the file!" << endl;
                }
                else {
                    cout << "Error in Reading File!!" << endl;
                }
                break;
            }
        }

        file.close();

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
#ifndef RECOMPRESSION_DEFINITIONS_HPP
#define RECOMPRESSION_DEFINITIONS_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <stack>
#include <cassert>

#include "typedefs.hpp"
#include "utils/space_efficient_vector.hpp"

using namespace std;

struct Node;

struct __attribute__((packed)) RLSLPNonterm {
    char_t type;
    c_size_t first;
    c_size_t second;
    c_size_t explen;

    RLSLPNonterm(const char_t& type, const c_size_t& first, const c_size_t& second);
    RLSLPNonterm();
    RLSLPNonterm(const char_t& type, const c_size_t& first, const c_size_t& second, const c_size_t& explen);
};

class RecompressionRLSLP {
public:
    space_efficient_vector<RLSLPNonterm> nonterm;

    uint64_t ram_use() const;
    void write_to_file(const string &filename);
    void read_from_file(const string &filename);
    void computeExplen();
    c_size_t lce(c_size_t i, c_size_t j);

private:
    void initialize_nodes(c_size_t node, const c_size_t& i, c_size_t left, c_size_t right, stack<Node>& ancestors, const space_efficient_vector<RLSLPNonterm>& grammar, Node& v);
    Node getLeftMostChild(Node v, const space_efficient_vector<RLSLPNonterm> & grammar);
    c_size_t getChildCount(const Node &parent, const space_efficient_vector<RLSLPNonterm> &grammar);
    c_size_t getChildIndex(const Node &parent, const Node &v, const space_efficient_vector<RLSLPNonterm> &grammar);
    Node getKthSibling(const Node &parent, const Node &v, c_size_t k);
    Node replaceWithHighestStartingAtPosition(const Node &v, stack<Node> &ancestors, const space_efficient_vector<RLSLPNonterm> &grammar);
    c_size_t LCE(Node &v1, Node &v2, c_size_t i, stack<Node> & v1_ancestors, stack<Node> & v2_ancestors, const space_efficient_vector<RLSLPNonterm> &grammar);
    c_size_t computeExplen(const c_size_t i);
};

struct SLGNonterm {
    c_size_t start_index;

    SLGNonterm(const c_size_t& start_index);
    SLGNonterm();
};

class SLG {
public:
    SLG();
    
    space_efficient_vector<SLGNonterm> nonterm;
    space_efficient_vector<c_size_t> rhs;

    uint64_t ram_use() const;
};

struct __attribute__((packed)) SLPNonterm {
    char_t type;
    c_size_t first;
    c_size_t second;

    SLPNonterm(const char_t& type, const c_size_t& first, const c_size_t& second);
    SLPNonterm();
};

class InputSLP {
public:
    space_efficient_vector<SLPNonterm> nonterm;

    InputSLP();
    InputSLP(const space_efficient_vector<SLPNonterm>& nonterm);

    void read_from_file(const string& file_name);

private:
    void order_slp();
};

struct Node {
    // Variable/Non-Terminal
    c_size_t var;
    // [l, r)
    c_size_t l;
    c_size_t r;

    Node();
    Node(const c_size_t& var, const c_size_t& l, const c_size_t& r);
};

struct __attribute__((packed)) AdjListElement {
    c_size_t first;
    c_size_t second;
    bool_t swapped;
    c_size_t vOcc;

    AdjListElement();
    AdjListElement(const c_size_t& first, const c_size_t& second, const bool_t& swapped, const c_size_t& vOcc);

    bool operator<(const AdjListElement& x) const;
};

#endif

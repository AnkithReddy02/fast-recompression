#include "recompression_definitions.hpp"

uint64_t RecompressionRLSLP::ram_use() const {
    return nonterm.ram_use();
}

void RecompressionRLSLP::write_to_file(const string &filename) {
    ofstream ofs(filename, ios::binary);

    if(!ofs) {
        cerr << "Error opening file for writing: " << filename << endl;
        return;
    }

    for(len_t i = 0; i < nonterm.size(); ++i) {
        const RLSLPNonterm &rlslpNonterm = nonterm[i];

        // Write type
        ofs.write(reinterpret_cast<const char*>(&rlslpNonterm.type), sizeof(rlslpNonterm.type));

        // Write first - 5 bytes
        ofs.write(reinterpret_cast<const char*>(&rlslpNonterm.first), 5);

        // Write second - 5 bytes
        ofs.write(reinterpret_cast<const char*>(&rlslpNonterm.second), 5);
    }

    ofs.close();
}

void RecompressionRLSLP::read_from_file(const string &filename) {
    ifstream ifs(filename, ios::binary);

    if(!ifs) {
        cerr << "Error opening file for reading: " << filename << endl;
        return;
    }

    nonterm.clear();

    while(ifs.peek() != EOF) {
        RLSLPNonterm rlslpNonterm;

        // Read type
        ifs.read(reinterpret_cast<char*>(&rlslpNonterm.type), sizeof(rlslpNonterm.type));

        // Read first - 5 bytes
        ifs.read(reinterpret_cast<char*>(&rlslpNonterm.first), 5);

        // Read second - 5 bytes
        ifs.read(reinterpret_cast<char*>(&rlslpNonterm.second), 5);

        nonterm.push_back(rlslpNonterm);
    }

    ifs.close();

    computeExplen();
}

c_size_t RecompressionRLSLP::computeExplen(const c_size_t i) {
    space_efficient_vector<RLSLPNonterm>& rlslp_nonterm_vec = nonterm;
    // Check if already computed
    if (rlslp_nonterm_vec[i].explen != 0) {
        return rlslp_nonterm_vec[i].explen;
    }

    switch (rlslp_nonterm_vec[i].type) {
        case '0':  // Terminal case
            return rlslp_nonterm_vec[i].explen = 1;

        case '1':  // Binary non-terminal
            return rlslp_nonterm_vec[i].explen = computeExplen(rlslp_nonterm_vec[i].first) +
                                                computeExplen(rlslp_nonterm_vec[i].second);

        default:  // Repetition case or other types
            return rlslp_nonterm_vec[i].explen = rlslp_nonterm_vec[i].second *
                                                computeExplen(rlslp_nonterm_vec[i].first);
    }

    // In case of unexpected type
    return 0; 
}

void RecompressionRLSLP::computeExplen() {
    for(len_t i = nonterm.size() - 1; i >= 1; --i) {
        computeExplen(i);
    }
}

uint64_t SLG::ram_use() const {
    return nonterm.ram_use() + rhs.ram_use();
}

SLG::SLG() {}

void RecompressionRLSLP::initialize_nodes(c_size_t node, const c_size_t& i, c_size_t left, c_size_t right, stack<Node>& ancestors, const space_efficient_vector<RLSLPNonterm>& grammar, Node& v) {
    if (left > right || i < left || i > right)
        return;

    const RLSLPNonterm& nt = grammar[node];

    if (left == i) {
        v = Node(node, left, right + 1);
        return;
    }

    ancestors.push(Node(node, left, right + 1));

    if (nt.type == '0') {
        return;
    } else if (nt.type == '1') {
        const RLSLPNonterm& left_nt = grammar[nt.first];
        const RLSLPNonterm& right_nt = grammar[nt.second];

        initialize_nodes(nt.first, i, left, left + left_nt.explen - 1, ancestors, grammar, v);
        initialize_nodes(nt.second, i, left + left_nt.explen, right, ancestors, grammar, v);

    } else {
        const RLSLPNonterm& child_nt = grammar[nt.first];
        c_size_t child_explen = child_nt.explen;

        c_size_t newLeft = left + ((i - left) / child_explen) * child_explen;
        c_size_t newRight = left + ((i - left) / child_explen) * child_explen + child_explen - 1;

        initialize_nodes(nt.first, i, newLeft, newRight, ancestors, grammar, v);
    }

    return;
}

Node RecompressionRLSLP::getLeftMostChild(Node v, const space_efficient_vector<RLSLPNonterm> & grammar) {
    char type = grammar[v.var].type;

    const RLSLPNonterm& nt = grammar[v.var];

    c_size_t left = v.l;
    c_size_t right = v.r - 1;

    if(type == '1') {
        const RLSLPNonterm& left_nt = grammar[nt.first];
        const RLSLPNonterm& right_nt = grammar[nt.second];
        return Node(nt.first, left, left + left_nt.explen);
    }
    else if(type == '2') {
        const RLSLPNonterm& child_nt = grammar[nt.first];
        c_size_t child_explen = child_nt.explen;
        return Node(nt.first, left, left + child_explen);
    }

    return Node();
}

c_size_t RecompressionRLSLP::getChildCount(const Node &parent, const space_efficient_vector<RLSLPNonterm> &grammar) {
    const RLSLPNonterm &nt = grammar[parent.var];
    switch (nt.type) {
        case '0': return 1;
        case '1': return 2;
        default: return nt.second;
    }
}

c_size_t RecompressionRLSLP::getChildIndex(const Node &parent, const Node &v, const space_efficient_vector<RLSLPNonterm> &grammar) {
    const RLSLPNonterm &nt = grammar[parent.var];

    switch (nt.type) {
        case '0': 
            return 1;

        case '1': 
            return (parent.l == v.l) ? 1 : 2;

        default:
            c_size_t span = v.r - v.l;
            return (v.l - parent.l) / span + 1;
    }
}

Node RecompressionRLSLP::getKthSibling(const Node &parent, const Node &v, c_size_t k) {
    c_size_t segmentLength = v.r - v.l;
    c_size_t newLeft = parent.l + (k - 1) * segmentLength;
    c_size_t newRight = newLeft + segmentLength;

    return Node(v.var, newLeft, newRight);
}

Node RecompressionRLSLP::replaceWithHighestStartingAtPosition(const Node &v, stack<Node> &ancestors, const space_efficient_vector<RLSLPNonterm> &grammar) {
    Node child = v;
    while (!ancestors.empty() && ancestors.top().r == v.r) {
        child = ancestors.top();
        ancestors.pop();
    }

    if (ancestors.empty()) {
        cout << "ANCESTORS is EMPTY -- Exiting" << endl;
        exit(1);
        return Node(); // Or some other appropriate action
    }

    Node ancestor = ancestors.top();
    const RLSLPNonterm &nt = grammar[ancestor.var];

    if (nt.type == '1') {
        return Node(nt.second, child.r, ancestor.r);
    } else {
        c_size_t childIndex = getChildIndex(ancestor, child, grammar);
        return getKthSibling(ancestor, child, childIndex + 1);
    }

    return Node();
}

c_size_t RecompressionRLSLP::LCE(Node &v1, Node &v2, c_size_t i, stack<Node> & v1_ancestors, stack<Node> & v2_ancestors, const space_efficient_vector<RLSLPNonterm> &grammar) {
    c_size_t exp_len_v1 = v1.r - v1.l;
    c_size_t exp_len_v2 = v2.r - v2.l;

    if(exp_len_v1 == 1 and exp_len_v1 == exp_len_v2 and v1.var != v2.var) {
        return v1.l - i;
    }
    else if(exp_len_v1 > exp_len_v2) {
        v1_ancestors.push(v1);
        v1 = getLeftMostChild(v1, grammar);
    }
    else if(exp_len_v1 < exp_len_v2) {
        v2_ancestors.push(v2);
        v2 = getLeftMostChild(v2, grammar);
    }
    else if(v1.var != v2.var) {
        v2_ancestors.push(v2);
        v1_ancestors.push(v1);
        v1 = getLeftMostChild(v1, grammar);
        v2 = getLeftMostChild(v2, grammar);
    }
    else {

        Node v1_parent = v1_ancestors.top();
        Node v2_parent = v2_ancestors.top();

        c_size_t j1 = getChildIndex(v1_parent, v1, grammar);
        c_size_t j2 = getChildIndex(v2_parent, v2, grammar);


        c_size_t d1 = getChildCount(v1_parent, grammar);
        c_size_t d2 = getChildCount(v2_parent, grammar);

        c_size_t lambda = min(d1 - j1, d2 - j2);

        if(lambda <= 1) {
            v1 = replaceWithHighestStartingAtPosition(v1, v1_ancestors, grammar);

            if(v2.r >= grammar.back().explen) {
                return v1.l - i;
            }

            v2 = replaceWithHighestStartingAtPosition(v2, v2_ancestors, grammar);
        }
        else {
            v1 = getKthSibling(v1_parent, v1, j1 + lambda);
            v2 = getKthSibling(v2_parent, v2, j2 + lambda);
        }
    }

    return LCE(v1, v2, i, v1_ancestors, v2_ancestors, grammar);;
}

c_size_t RecompressionRLSLP::lce(c_size_t i, c_size_t j) {

    if(i == j) {
        return nonterm.back().explen - i;
    }
    
    if(i > j) swap(i, j);

    Node v1, v2;
    stack<Node> v1_ancestors, v2_ancestors;
    // v1_ancestors.push(Node(grammar.size()-1, 0, 33));
    // v2_ancestors.push(Node(grammar.size()-1, 0, 33));
    initialize_nodes(nonterm.size() - 1, i, 0, nonterm.back().explen - 1, v1_ancestors, nonterm, v1);
    initialize_nodes(nonterm.size() - 1, j, 0, nonterm.back().explen - 1, v2_ancestors, nonterm, v2);

    return LCE(v1, v2, i, v1_ancestors, v2_ancestors, nonterm);
}

#ifndef LCE_QUERIES_H
#define LCE_QUERIES_H

#include<bits/stdc++.h>
#include "RecompressionDefinitions.hpp"
using namespace std;

void initialize_nodes(int node, const int& i, int left, int right, stack<Node>& ancestors, vector<RLSLPNonterm>& grammar, Node& v) {
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
        int child_explen = child_nt.explen;

        int newLeft = left + ((i - left) / child_explen) * child_explen;
        int newRight = left + ((i - left) / child_explen) * child_explen + child_explen - 1;

        initialize_nodes(nt.first, i, newLeft, newRight, ancestors, grammar, v);
    }

    return;
}

Node getLeftMostChild(Node v, vector<RLSLPNonterm> & grammar) {
    char type = grammar[v.var].type;

    const RLSLPNonterm& nt = grammar[v.var];

    int left = v.l;
    int right = v.r - 1;

    if(type == '1') {
        const RLSLPNonterm& left_nt = grammar[nt.first];
        const RLSLPNonterm& right_nt = grammar[nt.second];
        return Node(nt.first, left, left + left_nt.explen);
    }
    else if(type == '2') {
        const RLSLPNonterm& child_nt = grammar[nt.first];
        int child_explen = child_nt.explen;
        return Node(nt.first, left, left + child_explen);
    }

    return Node();
}


int getChildCount(Node parent, vector<RLSLPNonterm> & grammar) {

    RLSLPNonterm & nt = grammar[parent.var];

    if(nt.type == '0') {
        return 1;
    }
    else if(nt.type == '1') {
        return 2;
    }
    else {
        return nt.second;
    }
}

int getChildIndex(Node parent, Node v, vector<RLSLPNonterm> & grammar) {
    RLSLPNonterm & nt = grammar[parent.var];

    if(nt.type == '0') {
        return 1;
    }
    else if(nt.type == '1') {
        if(parent.l == v.l) {
            return 1;
        }
        else {
            return 2;
        }
    }
    else {
        return (v.l - parent.l)/(v.r - v.l) + 1;
    }

}

Node getKthSibling(Node parent, Node v, int k, vector<RLSLPNonterm> &grammar ) {
    RLSLPNonterm & nt = grammar[parent.var];

    int left = parent.l;
    int right = parent.r;

    // cout << "328: " << parent.var << ' ' << left << ' ' << right << ' ' << k << endl;


    return Node(v.var, left + (k-1)*(v.r - v.l), left + (k-1)*(v.r - v.l) + (v.r - v.l));
}

Node replaceWithHighestStartingAtPosition(Node v, stack<Node> &ancestors, vector<RLSLPNonterm> & grammar) {

    Node child = v;
    while(ancestors.empty() == false and ancestors.top().r == v.r) {
        child = ancestors.top();
        ancestors.pop();
    }

    Node ancestor = ancestors.top();

    // cout << "sz: " << ' ' << ancestor.var << endl;

    RLSLPNonterm nt = grammar[ancestor.var];

    // cout << nt.type << endl;


    if(nt.type == '1') {
        return Node(nt.second, child.r, ancestor.r);
    }
    else {
        int childIndex = getChildIndex(ancestor, child, grammar);

        // cout << "353: " << childIndex << endl;

        return getKthSibling(ancestor, child, childIndex + 1, grammar);
    }

    return Node();

}

int LCE(Node v1, Node v2, int i, stack<Node> & v1_ancestors, stack<Node> & v2_ancestors, vector<RLSLPNonterm> &grammar) {

    int exp_len_v1 = v1.r - v1.l;

    int exp_len_v2 = v2.r - v2.l;

    // cout << exp_len_v1 << ' ' << exp_len_v2 << endl;

    // cout << "NEW" << endl;
    // cout << v1.var << ' ' << v1.l << ' ' << v1.r << endl;
    // cout << v2.var << ' ' << v2.l << ' ' << v2.r << endl;
    // cout << "END" << endl;

    if(exp_len_v1 == 1 and exp_len_v1 == exp_len_v2 and v1.var != v2.var) {
        return v1.l - i;
    }
    else if(exp_len_v1 > exp_len_v2) {
        v1_ancestors.push(v1);
        v1 = getLeftMostChild(v1, grammar);
        // cout << "374" << ' ' << v1.var << ' ' << v1.l << ' ' << v1.r << endl;
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

        int j1 = getChildIndex(v1_parent, v1, grammar);
        int j2 = getChildIndex(v2_parent, v2, grammar);


        int d1 = getChildCount(v1_parent, grammar);
        int d2 = getChildCount(v2_parent, grammar);

        int lambda = min(d1 - j1, d2 - j2);

        // cout << "OK" << endl;

        // cout << j1 << ' ' << d1 << endl;
        // cout << j2 << ' ' << d2 << endl;

        // cout << lambda << endl;

        // cout << "OK END" << endl;


        if(lambda <= 1) {

            // cout << v1.l << ' ' << v1.r << ' ' << v2.l << ' ' << v2.r << endl;


            v1 = replaceWithHighestStartingAtPosition(v1, v1_ancestors, grammar);

            if(v2.r >= grammar.back().explen) {
                return v1.l - i;
            }

            v2 = replaceWithHighestStartingAtPosition(v2, v2_ancestors, grammar);

            // if(v2.r >= grammar.back().explen) {
            //     return v1.l - i;
            // }
        }
        else {
            v1 = getKthSibling(v1_parent, v1, j1 + lambda, grammar);
            v2 = getKthSibling(v2_parent, v2, j2 + lambda, grammar);
        }


    }

    return LCE(v1, v2, i, v1_ancestors, v2_ancestors, grammar);;
}


#endif
#ifndef LCE_QUERIES_H
#define LCE_QUERIES_H

#include <algorithm>
#include <iostream>
#include <stack>
#include <vector>

#include "recompression_definitions.hpp"
#include "typedefs.hpp"
#include "utils/space_efficient_vector.hpp"

using namespace std;

void initialize_nodes(
  c_size_t node, const c_size_t& i, c_size_t left, c_size_t right,
  stack<Node>& ancestors, const space_efficient_vector<RLSLPNonterm>& grammar,
  Node& v) {

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
    c_size_t newRight =
      left + ((i - left) / child_explen) * child_explen + child_explen - 1;

    initialize_nodes(nt.first, i, newLeft, newRight, ancestors, grammar, v);
  }

  return;
}

Node getLeftMostChild(
  Node v, const space_efficient_vector<RLSLPNonterm>& grammar) {

  char type = grammar[v.var].type;

  const RLSLPNonterm& nt = grammar[v.var];

  c_size_t left = v.l;
  c_size_t right = v.r - 1;

  if (type == '1') {
    const RLSLPNonterm& left_nt = grammar[nt.first];
    const RLSLPNonterm& right_nt = grammar[nt.second];
    return Node(nt.first, left, left + left_nt.explen);
  } else if(type == '2') {
    const RLSLPNonterm& child_nt = grammar[nt.first];
    c_size_t child_explen = child_nt.explen;
    return Node(nt.first, left, left + child_explen);
  }

  return Node();
}


c_size_t getChildCount(
  const Node& parent, const space_efficient_vector<RLSLPNonterm>& grammar) {
  const RLSLPNonterm& nt = grammar[parent.var];
  switch (nt.type) {
    case '0': return 1;
    case '1': return 2;
    default: return nt.second;
  }
}


c_size_t getChildIndex(
  const Node& parent,
  const Node& v,
  const space_efficient_vector<RLSLPNonterm>& grammar) {

  const RLSLPNonterm& nt = grammar[parent.var];

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

Node getKthSibling(const Node& parent, const Node& v, c_size_t k) {
  c_size_t segmentLength = v.r - v.l;
  c_size_t newLeft = parent.l + (k - 1) * segmentLength;
  c_size_t newRight = newLeft + segmentLength;

  return Node(v.var, newLeft, newRight);
}


Node replaceWithHighestStartingAtPosition(
  const Node& v,
  stack<Node>& ancestors,
  const space_efficient_vector<RLSLPNonterm>& grammar) {

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
  const RLSLPNonterm& nt = grammar[ancestor.var];

  if (nt.type == '1') {
    return Node(nt.second, child.r, ancestor.r);
  } else {
    c_size_t childIndex = getChildIndex(ancestor, child, grammar);
    return getKthSibling(ancestor, child, childIndex + 1);
  }

  return Node();
}


c_size_t LCE(
  Node& v1,
  Node& v2,
  c_size_t i,
  stack<Node>& v1_ancestors,
  stack<Node>& v2_ancestors,
  const space_efficient_vector<RLSLPNonterm>& grammar) {

  c_size_t exp_len_v1 = v1.r - v1.l;
  c_size_t exp_len_v2 = v2.r - v2.l;

  if (exp_len_v1 == 1 and exp_len_v1 == exp_len_v2 and v1.var != v2.var) {
    return v1.l - i;
  }
  else if (exp_len_v1 > exp_len_v2) {
    v1_ancestors.push(v1);
    v1 = getLeftMostChild(v1, grammar);
  }
  else if (exp_len_v1 < exp_len_v2) {
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

    if (lambda <= 1) {
      v1 = replaceWithHighestStartingAtPosition(v1, v1_ancestors, grammar);

      if (v2.r >= grammar.back().explen) {
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


#endif
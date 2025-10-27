#include "../include/recompression_definitions.hpp"
#include "../include/lce_queries.hpp"

Node::Node() {

}

Node::Node(
  const c_size_t& var,
  const c_size_t& l,
  const c_size_t& r) : 
  var(var),
  l(l),
  r(r) {

}

RLSLPNonterm::RLSLPNonterm(
  const char_t& type,
  const c_size_t& first,
  const c_size_t& second) :
  type(type),
  first(first),
  second(second),
  explen(0) {

}

RLSLPNonterm::RLSLPNonterm() :
  type('0'),
  first(0),
  second(0),
  explen(0) {

}

RLSLPNonterm::RLSLPNonterm(
  const char_t& type,
  const c_size_t& first,
  const c_size_t& second,
  const c_size_t& explen) : 
  type(type),
  first(first),
  second(second),
  explen(explen) {

}

uint64_t RecompressionRLSLP::ram_use() const {
    return nonterm.ram_use();
}

void RecompressionRLSLP::write_to_file(const string& filename) {
  ofstream ofs(filename, ios::binary);

  if (!ofs) {
    cerr << "Error opening file for writing: " << filename << endl;
    return;
  }

  for (len_t i = 0; i < nonterm.size(); ++i) {
    const RLSLPNonterm& rlslpNonterm = nonterm[i];

    // Write type
    ofs.write(reinterpret_cast<const char*>(&rlslpNonterm.type),
              sizeof(rlslpNonterm.type));

    // Write first
    ofs.write(reinterpret_cast<const char*>(&rlslpNonterm.first),
              sizeof(c_size_t));

    // Write second
    ofs.write(reinterpret_cast<const char*>(&rlslpNonterm.second),
              sizeof(c_size_t));
  }

  ofs.close();
}

void RecompressionRLSLP::read_from_file(const string& filename) {
    ifstream ifs(filename, ios::binary);

    if (!ifs) {
      cerr << "Error opening file for reading: " << filename << endl;
      return;
    }

    nonterm.clear();

    cout << "Reading from file: " << filename << endl;

    while (ifs.peek() != EOF) {
      RLSLPNonterm rlslpNonterm;

      // Read type
      ifs.read(reinterpret_cast<char*>(&rlslpNonterm.type),
               sizeof(rlslpNonterm.type));

      // Read first
      ifs.read(reinterpret_cast<char*>(&rlslpNonterm.first),
               sizeof(c_size_t));

      // Read second
      ifs.read(reinterpret_cast<char*>(&rlslpNonterm.second),
               sizeof(c_size_t));

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
        return rlslp_nonterm_vec[i].explen =
          computeExplen(rlslp_nonterm_vec[i].first) +
          computeExplen(rlslp_nonterm_vec[i].second);

      default:  // Repetition case or other types
        return rlslp_nonterm_vec[i].explen = 
          rlslp_nonterm_vec[i].second *
          computeExplen(rlslp_nonterm_vec[i].first);
    }

    // In case of unexpected type
    return 0; 
}

void RecompressionRLSLP::computeExplen() {
  for (len_t i = nonterm.size() - 1; i >= 1; --i) {
    computeExplen(i);
  }
}

void RecompressionRLSLP::initialize_nodes(
  c_size_t node, 
  const c_size_t& i,
  c_size_t left,
  c_size_t right,
  stack<Node>& ancestors,
  const space_efficient_vector<RLSLPNonterm>& grammar,
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

      initialize_nodes(
        nt.first, i, left, left + left_nt.explen - 1, ancestors, grammar, v);
      initialize_nodes(
        nt.second, i, left + left_nt.explen, right, ancestors, grammar, v);

    } else {
      const RLSLPNonterm& child_nt = grammar[nt.first];
      c_size_t child_explen = child_nt.explen;

      c_size_t newLeft =
        left + ((i - left) / child_explen) * child_explen;
      c_size_t newRight =
        left + ((i - left) / child_explen) * child_explen + child_explen - 1;

      initialize_nodes(nt.first, i, newLeft, newRight, ancestors, grammar, v);
    }

    return;
}

Node RecompressionRLSLP::getLeftMostChild(
  Node v, const space_efficient_vector<RLSLPNonterm>& grammar) {
  char type = grammar[v.var].type;

  const RLSLPNonterm& nt = grammar[v.var];

  c_size_t left = v.l;
  c_size_t right = v.r - 1;

  if(type == '1') {
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

c_size_t RecompressionRLSLP::getChildCount(
  const Node& parent, const space_efficient_vector<RLSLPNonterm>& grammar) {
  const RLSLPNonterm& nt = grammar[parent.var];
  switch (nt.type) {
    case '0': return 1;
    case '1': return 2;
    default: return nt.second;
  }
}

c_size_t RecompressionRLSLP::getChildIndex(
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

Node RecompressionRLSLP::getKthSibling(
  const Node& parent, const Node& v, c_size_t k) {
  c_size_t segmentLength = v.r - v.l;
  c_size_t newLeft = parent.l + (k - 1) * segmentLength;
  c_size_t newRight = newLeft + segmentLength;

  return Node(v.var, newLeft, newRight);
}

Node RecompressionRLSLP::replaceWithHighestStartingAtPosition(
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

c_size_t RecompressionRLSLP::LCE(
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
  } else if (exp_len_v1 > exp_len_v2) {
    v1_ancestors.push(v1);
    v1 = getLeftMostChild(v1, grammar);
  } else if (exp_len_v1 < exp_len_v2) {
    v2_ancestors.push(v2);
    v2 = getLeftMostChild(v2, grammar);
  } else if (v1.var != v2.var) {
    v2_ancestors.push(v2);
    v1_ancestors.push(v1);
    v1 = getLeftMostChild(v1, grammar);
    v2 = getLeftMostChild(v2, grammar);
  } else {
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
    } else {
      v1 = getKthSibling(v1_parent, v1, j1 + lambda);
      v2 = getKthSibling(v2_parent, v2, j2 + lambda);
    }
  }

  return LCE(v1, v2, i, v1_ancestors, v2_ancestors, grammar);;
}

c_size_t RecompressionRLSLP::lce(c_size_t i, c_size_t j) {
    if (i == j) {
      return nonterm.back().explen - i;
    }
    
    if (i > j) swap(i, j);

    Node v1, v2;
    stack<Node> v1_ancestors, v2_ancestors;
    // v1_ancestors.push(Node(grammar.size()-1, 0, 33));
    // v2_ancestors.push(Node(grammar.size()-1, 0, 33));
    initialize_nodes(
      nonterm.size() - 1, i, 0, nonterm.back().explen - 1, v1_ancestors,
      nonterm, v1);
    initialize_nodes(
      nonterm.size() - 1, j, 0, nonterm.back().explen - 1, v2_ancestors,
      nonterm, v2);

    return LCE(v1, v2, i, v1_ancestors, v2_ancestors, nonterm);
}


SLGNonterm::SLGNonterm(
  const c_size_t& start_index) :
  start_index(start_index) {

}

SLGNonterm::SLGNonterm() {

}

uint64_t SLG::ram_use() const {
  return nonterm.ram_use() + rhs.ram_use();
}

SLG::SLG() {

}

SLPNonterm::SLPNonterm(
  const char_t& type,
  const c_size_t& first,
  const c_size_t& second) :
  type(type),
  first(first),
  second(second) {

}

SLPNonterm::SLPNonterm() :
  type('0'),
  first(0),
  second(0) {

}

InputSLP::InputSLP() {

}

InputSLP::InputSLP(
  const space_efficient_vector<SLPNonterm>& nonterm) :
  nonterm(nonterm) {

}

void InputSLP::read_from_file(const string& file_name) {
  ifstream file(file_name, ios::binary);

  if (!file.is_open()) {
    cerr << "Error: Unable to open the file - " + file_name << endl;
    exit(1);
  }

  file.seekg(0, ios::end);
  size_t file_size = file.tellg();
  file.seekg(0, ios::beg);
  cout << "Number of Non-Terminals : " << file_size/(2 * sizeof(c_size_t))
       << endl;
  cout << "File size: " << file_size << " bytes" << endl;

  unsigned char buffer[2 * sizeof(c_size_t)];

  // ****
  for (uint64_t i = 0; i < (file_size/(2 * sizeof(c_size_t))); ++i) {
      nonterm.push_back(SLPNonterm());
  }

  c_size_t i = 0;

  while (file.read(reinterpret_cast<char*>(buffer), sizeof(buffer))) {
    // Process the first x-byte integer
    c_size_t value1 = 0;
    for (uint64_t i = 0; i < sizeof(c_size_t); ++i) {
      value1 |= static_cast<uint64_t>(buffer[i]) << (i * 8);
    }

    // Process the second x-byte integer
    c_size_t value2 = 0;
    for (uint64_t i = sizeof(c_size_t); i < 2 * sizeof(c_size_t); ++i) {
      value2 |=
        static_cast<uint64_t>(buffer[i]) << ((i - sizeof(c_size_t)) * 8);
    }

    if (value1 == 0) {
      nonterm[i++] = SLPNonterm('0', value2, 0);
    } else {
      nonterm[i++] = SLPNonterm('1', value1 - 1, value2 - 1);
    }
  }

  // Check if we have read less than c_size_t bytes in the last chunk
  c_size_t lastChunkSize = c_size_t(file.gcount());
  if (lastChunkSize > 0) {
      cout << "Read " << lastChunkSize << " bytes in the last chunk, "
           << "incomplete for a 64-bit value." << endl;
  }

  if (file.eof()) {
    cout << "File successfully read and stored SLP" << endl;
  } else {
    cout << "Error in reading file!" << endl;
  }

  file.close();

  order_slp();
}

void InputSLP::order_slp() {

  cout << "Ordering SLP..." << endl;
  const c_size_t grammar_size = nonterm.size();

  space_efficient_vector<uint8_t> inorder(grammar_size, 0);
  space_efficient_vector<c_size_t> old_new_map(grammar_size, 0);

  cout << "  Creating Reverse Graph..." << endl;
  // Create reversed graph.
  space_efficient_vector<c_size_t> vertex_ptr(grammar_size + 1, (c_size_t)0);
  space_efficient_vector<c_size_t> edges;
  {
    cout << "    Computing Vertex Degress..." << endl;
    // Step 1: Compute vertex degrees.
    for (c_size_t i = 0; i < nonterm.size(); ++i) {
      if (nonterm[i].type == '1') {
        ++vertex_ptr[nonterm[i].first];
        ++vertex_ptr[nonterm[i].second];
      }
    }

    cout << "    Computing Prefix Sum..." << endl;
    // Step 2: turn vertex_ptr into exclusive prefix sum.
    c_size_t degsum = 0;
    for (c_size_t i = 0; i < grammar_size; ++i) {
      c_size_t newsum = degsum + vertex_ptr[i];
      vertex_ptr[i] = degsum;
      degsum = newsum;
    }
    vertex_ptr[grammar_size] = degsum;

    cout << "    Resizing edges vector..." << endl;
    // Step 3: resize `edges' to accomodate all edges.
    edges.resize(degsum);

    cout << "    Computing Edges..." << endl;
    // Step 4: compute edges.
    for (c_size_t i = 0; i < nonterm.size(); ++i) {
      if (nonterm[i].type == '1') {
        edges[vertex_ptr[nonterm[i].first]++] = i;
        edges[vertex_ptr[nonterm[i].second]++] = i;
      }
    }

    cout << "    Restoring pointers to list begin..." << endl;
    // Step 5: restore pointers to adj list begin.
    for (c_size_t i = grammar_size; i > 0; --i)
    vertex_ptr[i] = vertex_ptr[i - 1];
    vertex_ptr[0] = 0;
  }

  cout << "  Computing inorder and Initializing queue..." << endl;
  queue<c_size_t> q;

  // Note: graph is reversed!
  for (c_size_t i = 0; i < nonterm.size(); i++) {
    if (nonterm[i].type == '1') {
      inorder[i] += 2;
    } else {
      q.push(i);
    }
  }

  c_size_t nonterminal_ptr = 0;

  cout << "  Performing BFS..." << endl;
  while (!q.empty()) {
    c_size_t u = q.front();
    q.pop();
    old_new_map[u] = nonterminal_ptr++;

    for (c_size_t i = vertex_ptr[u]; i < vertex_ptr[u+1]; ++i) {
      const c_size_t v = edges[i];
      inorder[v]--;
      if(inorder[v] == 0) {
        q.push(v);
      }
    }
  }

  edges.clear();
  vertex_ptr.clear();
  inorder.clear();

  cout << "  Creating Ordered Nonterm..." << endl;
  space_efficient_vector<SLPNonterm> ordered_nonterm(grammar_size);

  for (c_size_t i = 0; i < nonterm.size(); i++) {
    const char_t& type = nonterm[i].type;
    const c_size_t& first = nonterm[i].first;
    const c_size_t& second = nonterm[i].second;

    if(type == '0') {
      ordered_nonterm[old_new_map[i]] = SLPNonterm('0', first, second);
    }
    else {
      ordered_nonterm[old_new_map[i]] =
        SLPNonterm('1', old_new_map[first], old_new_map[second]);
    }
  }

  // *****
  // nonterm = move(ordered_nonterm);
  assert(nonterm.size() == ordered_nonterm.size());

  cout << "  Assigning ordered nonterm to the existing one..." << endl;
  for(c_size_t i = 0; i < grammar_size; ++i) {
      nonterm[i] = ordered_nonterm[i];
  }

  cout << "Ordered SLP!" << endl << endl;

  return;
}

AdjListElement::AdjListElement() {

}

AdjListElement::AdjListElement(
  const c_size_t& first,
  const c_size_t& second,
  const bool_t& swapped,
  const c_size_t& vOcc) :
  first(first),
  second(second),
  swapped(swapped),
  vOcc(vOcc) {

}

bool AdjListElement::operator<(const AdjListElement& x) const {
  return (first < x.first) || (first == x.first && second < x.second);
}

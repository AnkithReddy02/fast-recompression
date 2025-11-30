#ifndef BCOMP_HPP
#define BCOMP_HPP

extern uint64_t verbosity;

void combineFrequenciesInRange(
  const space_efficient_vector<packed_pair<c_size_t, c_size_t>>& vec,
  const len_t& lr_pointer,
  const len_t& rr_pointer,
  space_efficient_vector<SLGNonterm>& new_slg_nonterm_vec,
  /* map<pair<c_size_t, c_size_t>, c_size_t> &m */
  hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t>& m,
  RecompressionRLSLP *recompression_rlslp,
  space_efficient_vector<c_size_t>& new_rhs) {

  typedef packed_pair<c_size_t, c_size_t> pair_type;

  // Check if vector is empty
  if(vec.empty() || lr_pointer > rr_pointer) {
    new_slg_nonterm_vec.push_back((c_size_t)new_rhs.size());
    return;
  }

  c_size_t curr_new_rhs_size = new_rhs.size();

  // Iterate through the vector within the specified range
  c_size_t currNum = vec[lr_pointer].first;
  c_size_t currFreq = vec[lr_pointer].second;

  for(len_t i = lr_pointer + 1;
      i <= rr_pointer && i < (len_t)vec.size();
      ++i) {
    // Check if the current and previous elements have the same number
    if(vec[i].first == currNum) {
      // Merge frequencies
      currFreq += vec[i].second;
    }
    else {
      {
        const pair_type p(currNum, currFreq);

        // Frequency is >=2 so merge and push to rlslp
        // Merged variable is terminal for a new_slg, so it is marked negative
        if(p.second >= 2) {
          // if(m.find(p) == m.end()) {
          if(!m.find(p)) {
            // m[p] = recompression_rlslp->nonterm.size();
            m.insert(p, recompression_rlslp->nonterm.size());
            recompression_rlslp->nonterm.push_back(
              RLSLPNonterm('2', abs(p.first), p.second));
          }
          //  ** Negative **
          new_rhs.push_back(-m[p]);
        }
        else {
          new_rhs.push_back(p.first);
        }
      }

      // Move to the next number
      currNum = vec[i].first;
      currFreq = vec[i].second;
    }
  }

  const pair_type p(currNum, currFreq);

  // Frequency is >=2 so merge and push to rlslp
  // Merged variable is terminal for a new_slg, so it is marked negative
  if(p.second >= 2) {
    // if(m.find(p) == m.end()) {
    if(!m.find(p)) {
      // m[p] = recompression_rlslp->nonterm.size();
      m.insert(p, recompression_rlslp->nonterm.size());
      recompression_rlslp->nonterm.push_back(
        RLSLPNonterm('2', abs(p.first), p.second));
    }
    //  ** Negative **
    new_rhs.push_back(-m[p]);
  }
  else {
    new_rhs.push_back(p.first);
  }
  // rhs.shrink_to_fit();
  new_slg_nonterm_vec.push_back(curr_new_rhs_size);

  return;
}

// Block Compression
SLG* BComp(SLG *slg, RecompressionRLSLP *recompression_rlslp) {
  hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t> m;

  // Current slg non-term list
  space_efficient_vector<SLGNonterm>& slg_nonterm_vec = slg->nonterm;

  // Global RHS
  space_efficient_vector<c_size_t>& global_rhs = slg->rhs;

  // // New SLG that needs to be created by applying BComp
  SLG *new_slg = new SLG();

  // New SLG non-term list that new_slg needs.
  space_efficient_vector<SLGNonterm>& new_slg_nonterm_vec = new_slg->nonterm;
  new_slg_nonterm_vec.reserve(slg_nonterm_vec.size() + 1);
  space_efficient_vector<c_size_t>& new_rhs = new_slg->rhs;

  const c_size_t& grammar_size = slg_nonterm_vec.size();

  typedef packed_pair<c_size_t, char_t> small_pair_type;
  space_efficient_vector<small_pair_type> LR_vec(
    grammar_size, small_pair_type(-1, 0));
  space_efficient_vector<small_pair_type> RR_vec(
    grammar_size, small_pair_type(-1, 0));

  typedef packed_pair<c_size_t, c_size_t> pair_type;
  space_efficient_vector<pair_type> large_LR_vec;
  space_efficient_vector<pair_type> large_RR_vec;

  space_efficient_vector<pair_type> rhs_expansion;

  // We shall iterate throught each production rule in the increasing order of variable.
  // 'i' --> represents the variable.
  for(len_t i = 0; i < grammar_size; ++i) {
    #ifdef DEBUG_LOG
    if (!(i & ((1 << 19) - 1)))
      cout << fixed << setprecision(2) << "Progress: "
           << ((i+1)*(long double)100)/grammar_size << '\r' << flush;
    #endif

    if(verbosity >= 1) {
      if (!(i & ((1 << 19) - 1))) {
        cout << fixed << setprecision(2) << "  Progress: "
             << ((i+1)*(long double)100)/grammar_size << "%\r" << flush;
      }
    }

    SLGNonterm& slg_nonterm = slg_nonterm_vec[i];

    c_size_t start_index = slg_nonterm.start_index;
    c_size_t end_index =
      (i == grammar_size-1) ?
        c_size_t(global_rhs.size())-1 :
        c_size_t(slg_nonterm_vec[i+1].start_index) - 1;

    if(start_index > end_index) {
      new_slg_nonterm_vec.push_back((c_size_t)new_rhs.size());
      continue;
    }

    // Create the expansion of RHS.
    // space_efficient_vector<pair<c_size_t, c_size_t>> rhs_expansion;
    rhs_expansion.set_empty();

    // Compute the Expansion.
    for(len_t j = start_index; j <= end_index; ++j) {
      const c_size_t& rhs_symbol = global_rhs[j];

      // Single Terminal
      if(rhs_symbol < 0) {
        rhs_expansion.push_back(pair_type(rhs_symbol, 1));
        continue;
      }

      // **LR** of a current variable(rhs_symbol) in current RHS is **not empty**.
      if(static_cast<c_size_t>(LR_vec[rhs_symbol].second) != 0) {
        c_size_t k = static_cast<c_size_t>(LR_vec[rhs_symbol].second);
        if(k == 255) {
          pair_type search_element(rhs_symbol, 0);

          // c_size_t search_index = lower_bound(large_LR_vec, search_element);
          c_size_t search_index = large_LR_vec.lower_bound(search_element);

          assert(large_LR_vec[search_index].first == rhs_symbol);
          assert(search_index != large_LR_vec.size());

          k = large_LR_vec[search_index].second;
        }

        rhs_expansion.push_back(
          pair_type((c_size_t)LR_vec[rhs_symbol].first, (c_size_t)k));
      }

      c_size_t rhs_symbol_start_index =
        new_slg_nonterm_vec[rhs_symbol].start_index;
      c_size_t rhs_symbol_end_index =
        (rhs_symbol == (c_size_t)new_slg_nonterm_vec.size()-1) ?
          (c_size_t)new_rhs.size()-1 :
          new_slg_nonterm_vec[rhs_symbol+1].start_index - 1;

      // Cap is not empty --> in new SLG the variable(rhs_symbol) RHS is
      // not empty --> then Cap is not empty.
      if(rhs_symbol_start_index <= rhs_symbol_end_index) {
        rhs_expansion.push_back(pair_type(rhs_symbol, 0));
      }

      // **RR** of a current variable(rhs_symbol) in current RHS is **not empty**.
      if(static_cast<c_size_t>(RR_vec[rhs_symbol].second) != 0) {
        c_size_t k = static_cast<c_size_t>(RR_vec[rhs_symbol].second);
        if(k == 255) {
          pair_type search_element = pair_type(rhs_symbol, 0);

          c_size_t search_index = large_RR_vec.lower_bound(search_element);

          assert(large_RR_vec[search_index].first == rhs_symbol);
          assert(search_index != large_RR_vec.size());

          k = large_RR_vec[search_index].second;
        }

        rhs_expansion.push_back(
          pair_type((c_size_t)RR_vec[rhs_symbol].first, (c_size_t)k));
      }
    }

    //rhs.clear();
    //vector<int>().swap(rhs);

    // Compute LR
    len_t lr_pointer = 1;

    pair_type LR = rhs_expansion[0];

    while(lr_pointer < rhs_expansion.size()) {
      if(LR.first == rhs_expansion[lr_pointer].first) {
        LR.second += rhs_expansion[lr_pointer].second;
      }
      else {
        break;
      }

      lr_pointer++;
    }

    if(LR.second >= 255) {
      LR_vec[i] = small_pair_type((c_size_t)LR.first, (char_t)255);
      large_LR_vec.push_back(pair_type((c_size_t)i, (c_size_t)LR.second));
    }
    else {
      LR_vec[i] = small_pair_type((c_size_t)LR.first, (char_t)LR.second);
    }

    // Compute RR
    len_t rr_pointer = rhs_expansion.size()-2;

    pair_type RR = rhs_expansion.back();

    // Case 1 : Everything is consumed by lr_pointer
    if(lr_pointer == rhs_expansion.size()) {
      // Cap is empty; set Cap
      new_slg_nonterm_vec.push_back((c_size_t)new_rhs.size());
      // RR is empty; set RR
      RR_vec[i] = packed_pair<c_size_t, char_t>(-1, 0);
    }
    // Case 2 : There is room for RR
    else {
      // rr_pointer can max go upto lr_pointer
      while(lr_pointer <= rr_pointer) {
        // Check whether terminals can be merged
        // Obviously, only terminals!
        if(RR.first == rhs_expansion[rr_pointer].first) {
          RR.second += rhs_expansion[rr_pointer].second;
        }
        else {
          break;
        }

        // rr_pointer--;
        --rr_pointer;
      }

      // set RR
      if(RR.second >= 255) {
        RR_vec[i] = small_pair_type((c_size_t)RR.first, (char_t)255);
        large_RR_vec.push_back(pair_type(i, RR.second));
      }
      else {
        RR_vec[i] = small_pair_type((c_size_t)RR.first, (char_t)RR.second);
      }

      // Compress Cap(middle part)
      combineFrequenciesInRange(
        rhs_expansion, lr_pointer, rr_pointer, new_slg_nonterm_vec, m,
        recompression_rlslp, new_rhs);
    }
  }
  #ifdef DEBUG_LOG
  cout << endl;
  #endif

  if(verbosity >= 1) {
    cout << "  Progress: 100.0%" << endl;
  }

  // Add new Starting Variable to the new Grammar G'(G Prime).
  const c_size_t& start_var = slg_nonterm_vec.size() - 1;

  packed_pair<c_size_t, c_size_t> start_var_LR;
  if(static_cast<c_size_t>(LR_vec[start_var].second) == 255) {

    c_size_t k = static_cast<c_size_t>(LR_vec[start_var].second);
    pair_type search_arr((c_size_t)start_var, (c_size_t)0);

    c_size_t search_index = large_LR_vec.lower_bound(search_arr);

    assert(search_index != large_LR_vec.size());

    k = large_LR_vec[search_index].second;

    start_var_LR = pair_type((c_size_t)LR_vec[start_var].first, (c_size_t)k);
  }
  else {
    start_var_LR = pair_type(
      (c_size_t)LR_vec[start_var].first,
      static_cast<c_size_t>(LR_vec[start_var].second));
  }

  pair_type start_var_RR;
  if(static_cast<c_size_t>(RR_vec[start_var].second) == 255) {
    c_size_t k = static_cast<c_size_t>(RR_vec[start_var].second);
    pair_type search_arr(start_var, 0);

    c_size_t search_index = large_RR_vec.lower_bound(search_arr);

    assert(search_index != large_RR_vec.size());

    k = large_RR_vec[search_index].second;

    start_var_RR = pair_type((c_size_t)RR_vec[start_var].first, (c_size_t)k);
  }
  else {
    start_var_RR = pair_type(
      (c_size_t)RR_vec[start_var].first,
      static_cast<c_size_t>(RR_vec[start_var].second));
  }

  c_size_t curr_new_rhs_size = new_rhs.size();

  // This always holds true.
  if(start_var_LR.second != 0) {
    if(start_var_LR.second >= 2) {
      // if(m.find(start_var_LR) == m.end()) {
      if(!m.find(start_var_LR)) {
        // m[start_var_LR] = recompression_rlslp->nonterm.size();
        m.insert(start_var_LR, recompression_rlslp->nonterm.size());
        recompression_rlslp->nonterm.push_back(
          RLSLPNonterm('2', abs(start_var_LR.first), start_var_LR.second));
      }
      new_rhs.push_back(-m[start_var_LR]);
    }
    else {
      new_rhs.push_back(start_var_LR.first);
    }
  }

  c_size_t start_var_start_index =
    new_slg_nonterm_vec[start_var].start_index;
  c_size_t start_var_end_index =
    (start_var == (c_size_t)new_slg_nonterm_vec.size()-1) ?
      (c_size_t)curr_new_rhs_size-1 :
      new_slg_nonterm_vec[start_var+1].start_index - 1;

  if(start_var_start_index <= start_var_end_index) {
    new_rhs.push_back(start_var);
  }

  if(start_var_RR.second != 0) {
    if(start_var_RR.second >= 2) {
      // if(m.find(start_var_RR) == m.end()) {
      if(!m.find(start_var_RR)) {
        // m[start_var_RR] = recompression_rlslp->nonterm.size();
        m.insert(start_var_RR, recompression_rlslp->nonterm.size());
        recompression_rlslp->nonterm.push_back(
          RLSLPNonterm('2', abs(start_var_RR.first), start_var_RR.second));
      }
      new_rhs.push_back(-m[start_var_RR]);
    }
    else {
      new_rhs.push_back(start_var_RR.first);
    }
  }

  new_slg_nonterm_vec.push_back(curr_new_rhs_size);

  delete slg;

  return new_slg;
  // return new SLG(new_slg_nonterm_vec, new_rhs);
}


#endif // BCOMP_HPP

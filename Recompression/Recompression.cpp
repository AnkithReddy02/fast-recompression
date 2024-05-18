#include <bits/stdc++.h>
#include "recompression_definitions.hpp"
#include "../radix_sort.h"
#include "../io.h"
#include "utilities.hpp"
#include "lce_queries.hpp"
#include "hash_pair.hpp"
#include "typedefs.hpp"
#include "packed_pair.hpp"
#include "hash_table.hpp"

// #define TEST_SAMPLE

using namespace std;

template<>
std::uint64_t get_hash(const std::pair<c_size_t, c_size_t> &x) {
    return (std::uint64_t)x.first * (std::uint64_t)4972548694736365 +
           (std::uint64_t)x.second * (std::uint64_t)3878547385475748;
}

template<>
std::uint64_t get_hash(const c_size_t &x) {
    return (std::uint64_t)x * (std::uint64_t)1232545666776865;
}

template<>
std::string keyToString(const c_size_t &x) {
  return utils::intToStr(x);
}

template<>
std::string keyToString(const std::pair<c_size_t,c_size_t> &x) {
  return "std::pair(" + utils::intToStr(x.first) + ", " + utils::intToStr(x.second) + ")";
}

template<>
std::string valToString(const c_size_t &x) {
  return utils::intToStr(x);
}

template<>
std::string valToString(const bool &x) {
  return utils::intToStr((std::uint64_t)x);
}

#if 0
void combineFrequenciesInRange(const vector<pair<c_size_t, c_size_t>>& vec, const c_size_t &lr_pointer, const c_size_t &rr_pointer, vector<pair<c_size_t, c_size_t>> &result) {
    // Check if vector is empty
    if (vec.empty() || lr_pointer > rr_pointer)
        return;

    // Iterate through the vector within the specified range
    c_size_t currNum = vec[lr_pointer].first;
    c_size_t currFreq = vec[lr_pointer].second;

    for (size_t i = lr_pointer + 1; i <= rr_pointer && i < vec.size(); ++i) {
        // Check if the current and previous elements have the same number
        if (vec[i].first == currNum) {
            // Merge frequencies
            currFreq += vec[i].second;
        } else {
            // Add current pair to result vector
            result.emplace_back(currNum, currFreq);
            // Move to the next number
            currNum = vec[i].first;
            currFreq = vec[i].second;
        }
    }

    // Add the last pair to result vector
    result.emplace_back(currNum, currFreq);

    return;
}

vector<pair<c_size_t, c_size_t>> combineFrequenciesInRange(const vector<pair<c_size_t, c_size_t>>& vec, const c_size_t &lr_pointer, const c_size_t &rr_pointer) {
    vector<pair<c_size_t, c_size_t>> result;

    // Check if vector is empty
    if (vec.empty() || lr_pointer > rr_pointer)
        return result;

    // Iterate through the vector within the specified range
    c_size_t currNum = vec[lr_pointer].first;
    c_size_t currFreq = vec[lr_pointer].second;

    for (size_t i = lr_pointer + 1; i <= rr_pointer && i < vec.size(); ++i) {
        // Check if the current and previous elements have the same number
        if (vec[i].first == currNum) {
            // Merge frequencies
            currFreq += vec[i].second;
        } else {
            // Add current pair to result vector
            result.emplace_back(currNum, currFreq);
            // Move to the next number
            currNum = vec[i].first;
            currFreq = vec[i].second;
        }
    }

    // Add the last pair to result vector
    result.emplace_back(currNum, currFreq);

    return result;
}
#endif

void combineFrequenciesInRange(const vector<pair<c_size_t, c_size_t>>& vec, const c_size_t &lr_pointer, const c_size_t &rr_pointer, vector<SLGNonterm> &new_slg_nonterm_vec, /* map<pair<c_size_t, c_size_t>, c_size_t> &m */hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m, RecompressionRLSLP *recompression_rlslp, vector<c_size_t> &new_rhs) {
    // Check if vector is empty
    if (vec.empty() || lr_pointer > rr_pointer) {
        new_slg_nonterm_vec.push_back((c_size_t)new_rhs.size());
        return;
    }

    c_size_t curr_new_rhs_size = new_rhs.size();

    // Iterate through the vector within the specified range
    c_size_t currNum = vec[lr_pointer].first;
    c_size_t currFreq = vec[lr_pointer].second;

    for (size_t i = lr_pointer + 1; i <= rr_pointer && i < vec.size(); ++i) {
        // Check if the current and previous elements have the same number
        if (vec[i].first == currNum) {
            // Merge frequencies
            currFreq += vec[i].second;
        } else {
            
            {
                const pair<c_size_t, c_size_t> p = {currNum, currFreq};

                // Frequency is >=2 so merge and push to rlslp
                // Merged variable is terminal for a new_slg, so it is marked negative
                if(p.second >= 2) {
                    // if(m.find(p) == m.end()) {
                    if(!m.find(p)) {
                        // m[p] = recompression_rlslp->nonterm.size();
                        m.insert(p, recompression_rlslp->nonterm.size());
                        recompression_rlslp->nonterm.emplace_back('2', abs(p.first), p.second);
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

    const pair<c_size_t, c_size_t> p = {currNum, currFreq};

    // Frequency is >=2 so merge and push to rlslp
    // Merged variable is terminal for a new_slg, so it is marked negative
    if(p.second >= 2) {
        // if(m.find(p) == m.end()) {
        if(!m.find(p)) {    
            // m[p] = recompression_rlslp->nonterm.size();
            m.insert(p, recompression_rlslp->nonterm.size());
            recompression_rlslp->nonterm.emplace_back('2', abs(p.first), p.second);
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

template <typename T>
c_size_t lower_bound(const vector<T> &vec, const T &search_element) {
    c_size_t low = 0;
    c_size_t high = vec.size() - 1;

    while(high - low > 1) {
        c_size_t mid = low + (high - low)/2;

        if(vec[mid] < search_element) {
            low = mid + 1;
        }
        else {
            high = mid;
        }
    }

    if(search_element <= vec[low]) return low;
    if(search_element <= vec[high]) return high;


    return vec.size();
}

// Block Compression
SLG* BComp(SLG *slg, RecompressionRLSLP *recompression_rlslp, //map<pair<c_size_t, c_size_t>, c_size_t> & m) {
    hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> & m) {

    // Current slg non-term list
    vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

    // Global RHS
    vector<c_size_t> &global_rhs = slg->rhs;

    // // New SLG that needs to be created by applying BComp
    // SLG *new_slg = new SLG();

    // New SLG non-term list that new_slg needs.
    vector<SLGNonterm> new_slg_nonterm_vec; // = new_slg->nonterm;
    vector<c_size_t> new_rhs;

    const c_size_t &grammar_size = slg_nonterm_vec.size();

    vector<packed_pair<c_size_t, char_t>> LR_vec(grammar_size, packed_pair<c_size_t, char_t>(-1,0));
    vector<packed_pair<c_size_t, char_t>> RR_vec(grammar_size, packed_pair<c_size_t, char_t>(-1,0));

    vector<pair<c_size_t, c_size_t>> large_LR_vec;
    vector<pair<c_size_t, c_size_t>> large_RR_vec;

    // We shall iterate throught each production rule in the increasing order of variable.
    // 'i' --> represents the variable.
    for(c_size_t i = 0; i < grammar_size; i++) {
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (i == grammar_size-1) ? global_rhs.size()-1 : slg_nonterm_vec[i+1].start_index - 1;

        if(start_index > end_index) {
            new_slg_nonterm_vec.push_back((c_size_t)new_rhs.size());
            continue;
        }

        // Create the expansion of RHS.
        vector<pair<c_size_t, c_size_t>> rhs_expansion;
     
        // Compute the Expansion.
        for(c_size_t j = start_index; j <= end_index; j++) {

            const c_size_t &rhs_symbol = global_rhs[j];

            // Single Terminal
            if(rhs_symbol < 0) {
                rhs_expansion.emplace_back(rhs_symbol, 1);
                continue;
            }

            // **LR** of a current variable(rhs_symbol) in current RHS is **not empty**.
            if(static_cast<c_size_t>(LR_vec[rhs_symbol].second) != 0) {
                c_size_t k = static_cast<c_size_t>(LR_vec[rhs_symbol].second);
                if(k == 255) {
                    pair<c_size_t, c_size_t> search_element = {rhs_symbol, 0};

                    c_size_t search_index = lower_bound(large_LR_vec, search_element);

                    assert(large_LR_vec[search_index].first == rhs_symbol);
                    assert(search_index != large_LR_vec.size());

                    k = large_LR_vec[search_index].second;
                }

                rhs_expansion.push_back(make_pair((c_size_t)LR_vec[rhs_symbol].first, (c_size_t)k));
            }

            c_size_t rhs_symbol_start_index = new_slg_nonterm_vec[rhs_symbol].start_index;
            c_size_t rhs_symbol_end_index = (rhs_symbol == (c_size_t)new_slg_nonterm_vec.size()-1) ? new_rhs.size()-1 : new_slg_nonterm_vec[rhs_symbol+1].start_index - 1;

           
            // Cap is not empty --> in new SLG the variable(rhs_symbol) RHS is not empty --> then Cap is not empty.
            if(rhs_symbol_start_index <= rhs_symbol_end_index) {
                rhs_expansion.emplace_back(rhs_symbol, 0);
            }

            // **RR** of a current variable(rhs_symbol) in current RHS is **not empty**.
            if(static_cast<c_size_t>(RR_vec[rhs_symbol].second) != 0) {
                c_size_t k = static_cast<c_size_t>(RR_vec[rhs_symbol].second);
                if(k == 255) {
                    pair<c_size_t, c_size_t> search_element = {rhs_symbol, 0};

                    c_size_t search_index = lower_bound(large_RR_vec, search_element);

                    assert(large_RR_vec[search_index].first == rhs_symbol);
                    assert(search_index != large_RR_vec.size());

                    k = large_RR_vec[search_index].second;
                }

                rhs_expansion.push_back(make_pair((c_size_t)RR_vec[rhs_symbol].first, (c_size_t)k));
            }
        }

        //rhs.clear();
        //vector<int>().swap(rhs);

        // Compute LR
        c_size_t lr_pointer = 1;

        pair<c_size_t, c_size_t> LR = rhs_expansion[0];

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
            LR_vec[i] = packed_pair<c_size_t, char_t>(LR.first, 255);
            large_LR_vec.push_back(make_pair((c_size_t)i, (c_size_t)LR.second));
        }
        else {
            LR_vec[i] = packed_pair<c_size_t, char_t>(LR.first, LR.second);
        }

        // Compute RR
        c_size_t rr_pointer = rhs_expansion.size()-2;

        pair<c_size_t, c_size_t> RR = rhs_expansion.back();

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

                rr_pointer--;
            }

            // set RR
            if(RR.second >= 255) {
                RR_vec[i] = packed_pair<c_size_t, char_t>(RR.first, 255);
                large_RR_vec.push_back(make_pair(i, RR.second));
            }
            else {
                RR_vec[i] = packed_pair<c_size_t, char_t>(RR.first, RR.second);
            }

            // Compress Cap(middle part)
            combineFrequenciesInRange(rhs_expansion, lr_pointer, rr_pointer, new_slg_nonterm_vec, m, recompression_rlslp, new_rhs);
        }
    }

    // Add new Starting Variable to the new Grammar G'(G Prime).
    const c_size_t &start_var = slg_nonterm_vec.size()-1;

    pair<c_size_t, c_size_t> start_var_LR;
    if(static_cast<c_size_t>(LR_vec[start_var].second) == 255) {

        c_size_t k = static_cast<c_size_t>(LR_vec[start_var].second);
        pair<c_size_t, c_size_t> search_arr = make_pair((c_size_t)start_var, (c_size_t)0);
        
        c_size_t search_index = lower_bound(large_LR_vec, search_arr);

        assert(search_index != large_LR_vec.size());

        k = large_LR_vec[search_index].second;

        start_var_LR = make_pair((c_size_t)LR_vec[start_var].first, (c_size_t)k);
    }
    else {
        start_var_LR = make_pair((c_size_t)LR_vec[start_var].first, static_cast<c_size_t>(LR_vec[start_var].second));
    }

    pair<c_size_t, c_size_t> start_var_RR;
    if(static_cast<c_size_t>(RR_vec[start_var].second) == 255) {

        c_size_t k = static_cast<c_size_t>(RR_vec[start_var].second);
        pair<c_size_t, c_size_t> search_arr = make_pair(start_var, 0);
        
        c_size_t search_index = lower_bound(large_RR_vec, search_arr);

        assert(search_index != large_RR_vec.size());

        k = large_RR_vec[search_index].second;

        start_var_RR = make_pair((c_size_t)RR_vec[start_var].first, (c_size_t)k);
    }
    else {
        start_var_RR = make_pair((c_size_t)RR_vec[start_var].first, static_cast<c_size_t>(RR_vec[start_var].second));
    }

    c_size_t curr_new_rhs_size = new_rhs.size();

    // This always holds true.
    if(start_var_LR.second != 0) {
        if(start_var_LR.second >= 2) {
            // if(m.find(start_var_LR) == m.end()) {
            if(!m.find(start_var_LR)) {
                // m[start_var_LR] = recompression_rlslp->nonterm.size();
                m.insert(start_var_LR, recompression_rlslp->nonterm.size());
                recompression_rlslp->nonterm.emplace_back('2', abs(start_var_LR.first), start_var_LR.second);
            }
            new_rhs.push_back(-m[start_var_LR]);
        }
        else {
            new_rhs.push_back(start_var_LR.first);
        }
    }

    c_size_t start_var_start_index = new_slg_nonterm_vec[start_var].start_index;
    c_size_t start_var_end_index = (start_var == (c_size_t)new_slg_nonterm_vec.size()-1) ? curr_new_rhs_size-1 : new_slg_nonterm_vec[start_var+1].start_index - 1;

    if(start_var_start_index <= start_var_end_index) {
        new_rhs.push_back(start_var);
    }

    if(start_var_RR.second != 0) {
        if(start_var_RR.second >= 2) {
            // if(m.find(start_var_RR) == m.end()) {
            if(!m.find(start_var_RR)) {
                // m[start_var_RR] = recompression_rlslp->nonterm.size();
                m.insert(start_var_RR, recompression_rlslp->nonterm.size());
                recompression_rlslp->nonterm.emplace_back('2', abs(start_var_RR.first), start_var_RR.second);
            }
            new_rhs.push_back(-m[start_var_RR]);
        }
        else {
            new_rhs.push_back(start_var_RR.first);
        }
    }

    new_slg_nonterm_vec.push_back(curr_new_rhs_size);

    delete slg;

    // return new_slg;
    return new SLG(new_slg_nonterm_vec, new_rhs);
}

#if 0
pair<c_size_t, c_size_t> computeAdjListHelper(c_size_t var, SLG *slg, vector<AdjListElement> & adjList, vector<pair<c_size_t, c_size_t>> & dp, vector<c_size_t> &vOcc) {

    if(var < 0) {
        return {var, var};
    }

    if(dp[var].first != 1 && dp[var].second != 1) {
        return dp[var];
    }

    // 1. Modifying LMS and RMS
    // 2. Access rhs
    SLGNonterm &slg_nonterm = slg->nonterm[var];

    const vector<c_size_t> &global_rhs = slg->rhs;
    const vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

    c_size_t start_index = slg_nonterm.start_index;
    c_size_t end_index = (var == slg_nonterm_vec.size()-1) ? global_rhs.size()-1 : slg_nonterm_vec[var+1].start_index - 1;

    c_size_t curr_rhs_size = end_index - start_index + 1;

    // 
    if(curr_rhs_size == 0) {
        return {0, 0};
    }
    else if(curr_rhs_size == 1 && global_rhs[start_index] < 0) {
        // Hopefully this should be always true when there is only 1 character to right
        // assert(rhs[0] < 0);
        return dp[var] = {global_rhs[start_index], global_rhs[start_index]};

        // return dp[var] = {rhs[0], rhs[0]};
    }
    
    vector<pair<c_size_t, c_size_t>> lms_rms_list;

    for(c_size_t j=start_index; j<=end_index; j++) {
        const c_size_t &rhs_symbol = global_rhs[j];
        pair<c_size_t, c_size_t> lms_rms = computeAdjListHelper(rhs_symbol, slg, adjList, dp, vOcc);
        lms_rms_list.push_back(lms_rms);
    }

    for(c_size_t i = 0; i < lms_rms_list.size() - 1; i++) {

        c_size_t f = lms_rms_list[i].second;
        c_size_t s = lms_rms_list[i+1].first;

        bool_t swapped = false;

        if(abs(f) < abs(s)) {
            swap(f, s);
            swapped = true;
        }
        adjList.push_back(AdjListElement(f, s, swapped, vOcc[var]));
    }

    //slg_nonterm.LMS = lms_rms_list.front().first;
    //slg_nonterm.RMS = lms_rms_list.back().second;

    return dp[var] = {lms_rms_list.front().first, lms_rms_list.back().second};
}

void computeAdjList(SLG *slg, vector<AdjListElement> &adjList, vector<c_size_t> &vOcc) {
    vector<pair<c_size_t, c_size_t>> dp(slg->nonterm.size(), make_pair(1, 1));

    computeAdjListHelper(slg->nonterm.size()-1, slg, adjList, dp, vOcc);
    return;
}

vector<AdjListElement> computeAdjList(SLG *slg, vector<c_size_t> &vOcc) {

    vector<AdjListElement> adjList;

    vector<pair<c_size_t, c_size_t>> dp(slg->nonterm.size(), make_pair(1, 1));

    for(c_size_t i = 0; i < slg->nonterm.size(); ++i) {
        computeAdjListHelper(i, slg, adjList, dp, vOcc);
    }

    return adjList;
}
#endif

c_size_t computeVOccHelper(vector<pair<c_size_t, c_size_t>> & edges, vector<c_size_t> &curr_index, vector<bool_t> &have_edges, c_size_t u, vector<c_size_t> & dp) {

    // Base Case : Target is Reached / Target is Same as the current node.
    if(u == curr_index.size()-1) {
        return dp[u] = 1;
    }

    // Already Computed
    if(dp[u] != -1) {
        return dp[u];
    }

    c_size_t num_paths = 0;

    c_size_t u_curr_index = curr_index[u];

    while(true) {

       
        if(u_curr_index >= edges.size() || have_edges[u] == false || (u < curr_index.size()-1 && u_curr_index == curr_index[u+1])) {
            break;
        }

        // cout << u_curr_index << ' ' << curr_index[u+1] << endl;


        const auto &edge = edges[u_curr_index++];
        const c_size_t &v = edge.first;
        const c_size_t &weight = edge.second;

        
        num_paths += weight * computeVOccHelper(edges, curr_index, have_edges, v, dp);
        // cout << u << ' ' << v << ' ' << weight << endl;
        // cout << num_paths << endl;
    }

    return dp[u] = num_paths;
}

void computeVOcc(SLG *slg, vector<c_size_t> &dp) {

    // Here, any list size is slg_nonterm_vec.size()
    vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;
    vector<c_size_t> &global_rhs = slg->rhs;

    // Create Reverse Graph.
    // vector<vector<pair<c_size_t,c_size_t>>> graph(slg_nonterm_vec.size(), vector<pair<c_size_t, c_size_t>>());
    vector<c_size_t> curr_index(slg_nonterm_vec.size(), -1);
    vector<pair<c_size_t, c_size_t>> edges;

    for(c_size_t i=0; i<slg_nonterm_vec.size(); i++) {

        // Current SLGNonterm
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (i == slg_nonterm_vec.size()-1) ? global_rhs.size()-1 : slg_nonterm_vec[i+1].start_index - 1;

        c_size_t curr_rhs_size = end_index - start_index + 1;

        set<c_size_t> unique_var;
        // hash_table<c_size_t, bool> unique_var;

        // Enumerate each character of RHS
        // Frequency calculation for reverse of the graph.
        for(c_size_t j=start_index; j<=end_index; j++) {
            const c_size_t & var = global_rhs[j];
            // Only Non-Terminals
            if(var >= 0) {
                unique_var.insert(var);
                // unique_var.insert(var, true);
            }
        }

        // Construct Edges.
        for(const auto & v : unique_var) {
        // for(c_size_t i = 0; i < unique_var.size(); ++i) {
            // c_size_t v = unique_var.get(i).first;
            // Weighted Edge from v to u , v --> u
            if(curr_index[v]==-1) curr_index[v] = 0;
            curr_index[v]++;
            
        }
    }

    c_size_t prefix_sum = 0;
    for(c_size_t i=0; i<curr_index.size(); i++) {
        if(curr_index[i] > 0) {
            c_size_t temp = curr_index[i];
            curr_index[i] += prefix_sum;
            prefix_sum += temp;
        }
    }

    edges.resize(prefix_sum);

    // Enumerate Each Production Rule.
    for(c_size_t i=0; i<slg_nonterm_vec.size(); i++) {
        // Current SLGNonterm
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (i == slg_nonterm_vec.size()-1) ? global_rhs.size()-1 : slg_nonterm_vec[i+1].start_index - 1;

        c_size_t curr_rhs_size = end_index - start_index + 1;

        unordered_map<c_size_t, c_size_t> var_freq;
        // hash_table<c_size_t, c_size_t, c_size_t> var_freq;

        // Enumerate each character of RHS
        // Frequency calculation for reverse of the graph.
        for(c_size_t j=start_index; j<=end_index; j++) {

            const c_size_t & var = global_rhs[j];
            // Only Non-Terminals
            if(var >= 0) {
                var_freq[var]++;
                // if(!var_freq.find(var)) {
                //     var_freq.insert(var, 0);
                // }
                // var_freq.insert(var, var_freq[var] + 1);
            }
        }

        // Construct Edges.
        for(const auto & x : var_freq) {
        // for(c_size_t j = 0; j < var_freq.size(); ++j) {
            // auto x = var_freq.get(j); 

            // edge from v to u, u <-- v
            const c_size_t &u = i;
            const c_size_t &v = x.first;
            const c_size_t &weight = x.second;

            // Weighted Edge from v to u , v --> u
            // graph[v].emplace_back(u, weight);
            edges[--curr_index[v]] = make_pair(u, weight);
            // if(curr_index[v]==-1) curr_index[v] = 0;
            // curr_index[v]++;
            
        }
    }

    vector<bool_t> have_edges(curr_index.size(), false);

    for(c_size_t i=curr_index.size()-1; i>=0; i--) {

        if(curr_index[i] == -1) {
            have_edges[i] = false;
        }
        else {
            have_edges[i] = true;
        }
        if(i <= curr_index.size()-2) {
            if(curr_index[i] == -1 && curr_index[i+1] != -1) {
                curr_index[i] = curr_index[i+1];
            }
        }
    }

    // Store VOcc.
    dp.resize(slg_nonterm_vec.size(), -1);

    // Compute vOcc.
    for(c_size_t i = slg_nonterm_vec.size() - 1; i >= 0; i--) {
        computeVOccHelper(edges, curr_index, have_edges, i, dp);
    }

    return;
}

bool AdjListElementCompare(const AdjListElement &a, const AdjListElement &b) {
    if(a.first != b.first) {
        return a.first < b.first;
    }

    return a.second < b.second;
}

void sortAdjList(vector<AdjListElement> & adjList) {

    // Make Positive
    // for(array<int, 4> & a : adjList) {
    //     a[0] = -a[0];
    //     a[1] = -a[1];
    // }

    for(AdjListElement & a : adjList) {
        a.first = -a.first;
        a.second = -a.second;
    }

    // Sort It.
    // radixSort(adjList);
    sort(adjList.begin(), adjList.end(), AdjListElementCompare);


    // Revert to Negative
    for(AdjListElement & a : adjList) {
        a.first = -a.first;
        a.second = -a.second;
    }
}

// void createPartition(const vector<AdjListElement> & adjList, array<set<c_size_t>, 2> &partition_set) {
pair<hash_table<c_size_t, bool> *, hash_table<c_size_t, bool> *> createPartition(const vector<AdjListElement> & adjList) {

    // // Make Positive
    // for(array<int, 4> & a : adjList) {
    //  a[0] = -a[0];
    //  a[1] = -a[1];
    // }
    
    typedef c_size_t key_type;
    typedef bool value_type;
    typedef hash_table<key_type, value_type> hash_table_type;
    // set<c_size_t>& leftSet = partition_set[0];
    hash_table_type *leftSet_ptr = new hash_table_type();
    // set<c_size_t>& rightSet = partition_set[1];
    hash_table_type *rightSet_ptr = new hash_table_type();

    hash_table_type &leftSet = *leftSet_ptr;
    hash_table_type &rightSet = *rightSet_ptr;

    c_size_t currentIndex = 0;
    size_t n = adjList.size();

    c_size_t c = adjList[currentIndex].first;

    while(currentIndex < n) {
        c_size_t leftSetFreq = 0;
        c_size_t rightSetFreq = 0;
        while (currentIndex < n && adjList[currentIndex].first == c) {
            if(!rightSet.find(adjList[currentIndex].second) /*== rightSet.end()*/) {
                leftSet.insert(adjList[currentIndex].second, true);
            }

            if (leftSet.find(adjList[currentIndex].second) /*!= leftSet.end()*/) {
                leftSetFreq += adjList[currentIndex].vOcc;
            } else {
                rightSetFreq += adjList[currentIndex].vOcc;
            }
            currentIndex++;
        }


        if (leftSetFreq >= rightSetFreq) {
            rightSet.insert(c, true);
        } else {
            leftSet.insert(c, true);
        }

        if(currentIndex < n) {
            c = adjList[currentIndex].first;
        }
        
    }
    
    c_size_t LRPairsCount = 0;
    c_size_t RLPairsCount = 0;

    /*
        for (int i = 0; i < arr.size() - 1; i++) {


            LRPairsCount += (leftSet.find(arr[i]) != leftSet.end()) && (rightSet.find(arr[i + 1]) != rightSet.end());
            RLPairsCount += (rightSet.find(arr[i]) != rightSet.end()) && (leftSet.find(arr[i + 1]) != leftSet.end());
        }
    */

    for(const AdjListElement & arr : adjList) {
        c_size_t f = arr.first;
        c_size_t s = arr.second;

        if(arr.swapped == true) {
            swap(f, s);
        }

        LRPairsCount += ((leftSet.find(f) /*!= leftSet.end()*/) && (rightSet.find(s) /*!= rightSet.end()*/)) ? arr.vOcc : 0;
        RLPairsCount += ((rightSet.find(f) /*!= rightSet.end()*/) && (leftSet.find(s) /*!= leftSet.end()*/)) ? arr.vOcc : 0;

    }

    if (RLPairsCount < LRPairsCount) {
        // swap(leftSet, rightSet);
        swap(leftSet_ptr, rightSet_ptr);
    }

    // swap(rightSet, leftSet);
    swap(rightSet_ptr, leftSet_ptr);

 //    // Revert to Negative
    // for(array<int, 4> & a : adjList) {
    //  a[0] = -a[0];
    //  a[1] = -a[1];
    // }

    

    return make_pair(leftSet_ptr, rightSet_ptr);
}

#if 0
array<set<c_size_t>, 2> createPartition(const vector<array<c_size_t, 4>> & adjList) {

    // // Make Positive
    // for(array<int, 4> & a : adjList) {
    //  a[0] = -a[0];
    //  a[1] = -a[1];
    // }
    
    set<c_size_t> leftSet, rightSet;
    c_size_t currentIndex = 0;
    size_t n = adjList.size();

    c_size_t c = adjList[currentIndex][0];

    while(currentIndex < n) {
        c_size_t leftSetFreq = 0;
        c_size_t rightSetFreq = 0;
        while (currentIndex < n && adjList[currentIndex][0] == c) {
            if(rightSet.find(adjList[currentIndex][1]) == rightSet.end()) {
                leftSet.insert(adjList[currentIndex][1]);
            }

            if (leftSet.find(adjList[currentIndex][1]) != leftSet.end()) {
                leftSetFreq += adjList[currentIndex][3];
            } else {
                rightSetFreq += adjList[currentIndex][3];
            }
            currentIndex++;
        }


        if (leftSetFreq >= rightSetFreq) {
            rightSet.insert(c);
        } else {
            leftSet.insert(c);
        }

        if(currentIndex < n) {
            c = adjList[currentIndex][0];
        }
        
    }
    
    c_size_t LRPairsCount = 0;
    c_size_t RLPairsCount = 0;

    /*
        for (c_size_t i = 0; i < arr.size() - 1; i++) {


            LRPairsCount += (leftSet.find(arr[i]) != leftSet.end()) && (rightSet.find(arr[i + 1]) != rightSet.end());
            RLPairsCount += (rightSet.find(arr[i]) != rightSet.end()) && (leftSet.find(arr[i + 1]) != leftSet.end());
        }
    */

    for(const array<c_size_t, 4> & arr : adjList) {
        c_size_t f = arr[0];
        c_size_t s = arr[1];

        if(arr[2] == 1) {
            swap(f, s);
        }

        LRPairsCount += ((leftSet.find(f) != leftSet.end()) && (rightSet.find(s) != rightSet.end())) ? arr[3] : 0;
        RLPairsCount += ((rightSet.find(f) != rightSet.end()) && (leftSet.find(s) != leftSet.end())) ? arr[3] : 0;

    }

    if (RLPairsCount < LRPairsCount) {
        swap(leftSet, rightSet);
    }

 //    // Revert to Negative
    // for(array<int, 4> & a : adjList) {
    //  a[0] = -a[0];
    //  a[1] = -a[1];
    // }

    

    return { rightSet, leftSet };
}
#endif

// MAP USAGE
pair<c_size_t, c_size_t> computeAdjListHelper(c_size_t var, SLG *slg, //map<pair<c_size_t, c_size_t>, c_size_t> &m0, map<pair<c_size_t, c_size_t>, c_size_t> &m1, vector<pair<c_size_t, c_size_t>> & dp, vector<c_size_t> &vOcc) {
hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m0, hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m1, vector<pair<c_size_t, c_size_t>> & dp, vector<c_size_t> &vOcc) {

    if(var < 0) {
        return {var, var};
    }

    if(dp[var].first != 1 && dp[var].second != 1) {
        return dp[var];
    }

    /*
        1. Modifying LMS and RMS
        2. Access rhs
    */
    SLGNonterm &slg_nonterm = slg->nonterm[var];

    const vector<c_size_t> & global_rhs = slg->rhs;
    vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

    c_size_t start_index = slg_nonterm.start_index;
    c_size_t end_index = (var == slg_nonterm_vec.size()-1) ? (c_size_t)global_rhs.size()-1 : slg_nonterm_vec[var+1].start_index - 1;

    c_size_t curr_rhs_size = end_index - start_index + 1;


    // 
    if(curr_rhs_size == 0) {
        return {0, 0};
    }
    else if(curr_rhs_size == 1 && global_rhs[start_index] < 0) {
        // Hopefully this should be always true when there is only 1 character to right
        // assert(rhs[0] < 0);
        return dp[var] = {global_rhs[start_index], global_rhs[start_index]};

        // return dp[var] = {rhs[0], rhs[0]};
    }
    
    vector<pair<c_size_t, c_size_t>> lms_rms_list;

    

    for(c_size_t j = start_index; j <= end_index; j++) {
        const c_size_t &rhs_symbol = global_rhs[j];
        pair<c_size_t, c_size_t> lms_rms = computeAdjListHelper(rhs_symbol, slg, m0, m1, dp, vOcc);
        lms_rms_list.push_back(lms_rms);
    }

    for(c_size_t i = 0; i < lms_rms_list.size() - 1; i++) {

        c_size_t f = lms_rms_list[i].second;
        c_size_t s = lms_rms_list[i+1].first;

        bool_t swapped = false;

        if(abs(f) < abs(s)) {
            swap(f, s);
            swapped = true;
        }

        if(swapped) {
            // m1[{f, s}] += vOcc[var];
            pair<c_size_t, c_size_t> p = {f, s};

            if(!m1.find(p)) {
                m1.insert(p, 0);
            }

            m1.insert(p, m1[p] + vOcc[var]);
        }
        else {
            // m0[{f, s}] += vOcc[var];
            pair<c_size_t, c_size_t> p = {f, s};

            if(!m0.find(p)) {
                m0.insert(p, 0);
            }

            m0.insert(p, m0[p] + vOcc[var]);
        }
        // adjList.push_back(AdjListElement(f, s, swapped, vOcc[var]));
    }

    //slg_nonterm.LMS = lms_rms_list.front().first;
    //slg_nonterm.RMS = lms_rms_list.back().second;

    return dp[var] = {lms_rms_list.front().first, lms_rms_list.back().second};
}

void computeAdjList(SLG *slg, //map<pair<c_size_t, c_size_t>, c_size_t> &m0, map<pair<c_size_t, c_size_t>, c_size_t> &m1, vector<c_size_t> &vOcc) {
    hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m0, hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m1, vector<c_size_t> &vOcc) {

    vector<pair<c_size_t, c_size_t>> dp(slg->nonterm.size(), make_pair(1, 1));

    computeAdjListHelper(slg->nonterm.size()-1, slg, m0, m1, dp, vOcc);
    return;
}

// Pair-Wise Compression
SLG * PComp(SLG *slg, RecompressionRLSLP *recompression_rlslp,  //map<pair<c_size_t, c_size_t>, c_size_t> & m) { 
    hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> & m) { 
    vector<c_size_t> vOcc;
    // Compute vOcc
    computeVOcc(slg, vOcc);

    // Compute AdjList
    hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> m0;
    hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> m1;
    computeAdjList(slg, m0, m1, vOcc);
    // Avoid resizing.
    vector<AdjListElement> adjList(m0.size() + m1.size());

    assert(adjList.size() != 0);

    c_size_t index = 0;
    // for(auto & x : m0) {
    //     adjList[index++] = AdjListElement(x.first.first, x.first.second, false, x.second);
    // }
    for(c_size_t i = 0; i < m0.size(); i++) {
        pair<pair<c_size_t, c_size_t>, c_size_t> x = m0.get(i);
        adjList[index++] = AdjListElement(x.first.first, x.first.second, false, x.second);
    }

    // m0.clear();

    // for(auto & x : m1) {
    //     adjList[index++] = AdjListElement(x.first.first, x.first.second, true, x.second);
    // }
    for(c_size_t i = 0; i < m1.size(); i++) {
        pair<pair<c_size_t, c_size_t>, c_size_t> x = m1.get(i);
        adjList[index++] = AdjListElement(x.first.first, x.first.second, true, x.second);
    }

    // m1.clear();

    sortAdjList(adjList);

    // Create Partition.
    // array<set<c_size_t>,2> arr;
    // createPartition(adjList, arr);
    typedef c_size_t key_type;
    typedef bool value_type;
    typedef hash_table<key_type, value_type> hash_table_type;
    pair<hash_table_type*, hash_table_type*> partition = createPartition(adjList);

    adjList.clear();  // Clear immediately
    vector<AdjListElement>().swap(adjList);  // Ensure memory is released

    // const set<c_size_t> &left_set = arr[0], &right_set = arr[1];
    const hash_table_type &left_set = *partition.first;
    const hash_table_type &right_set = *partition.second;

    // Current slg non-term list
    vector<SLGNonterm> &slg_nonterm_vec = slg->nonterm;

    // Current rhs
    vector<c_size_t> &global_rhs = slg->rhs;

    // New SLG that needs to be created by applying BComp
    // SLG *new_slg = new SLG();

    // New SLG non-term list that new_slg needs.
    vector<SLGNonterm> new_slg_nonterm_vec; // = new_slg->nonterm;

    // New RHS
    vector<c_size_t> new_rhs;

    vector<c_size_t> LB(slg_nonterm_vec.size(), 0), RB(slg_nonterm_vec.size(), 0);

    // We shall iterate throught each production rule in the increasing order of variable.
    // 'i' --> represents the variable.
    for(c_size_t i=0; i<slg_nonterm_vec.size(); i++) {

        // cout << i << " : ";
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        // const vector<c_size_t> & rhs = slg_nonterm.rhs;
        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (i == slg_nonterm_vec.size()-1) ? global_rhs.size()-1 : slg_nonterm_vec[i+1].start_index - 1;

        c_size_t curr_rhs_size = end_index - start_index + 1;

        c_size_t curr_new_rhs_size = new_rhs.size();

        if(curr_rhs_size == 0) {
            new_slg_nonterm_vec.push_back((c_size_t)new_rhs.size());
            // cout << '\n';
            continue;
        }
        if((c_size_t)curr_rhs_size >= 2) {

            vector<c_size_t> rhs_expansion;

            // Expanding RHS.
            for(c_size_t j=start_index; j<=end_index; j++) {

                // cout << global_rhs[j] << ' ';
                // RHS Terminal
                if(global_rhs[j] < 0) {

                    rhs_expansion.push_back(global_rhs[j]);

                    if(j==start_index) {
                        if(left_set.find(global_rhs[j]) /*!= left_set.end()*/) {
                            LB[i] = 0;
                        }
                        else if(right_set.find(global_rhs[j]) /*!= right_set.end()*/) {
                            LB[i] = global_rhs[j];

                            // Remove the pushed element
                            rhs_expansion.pop_back();
                        }
                    }
                    else if(j == end_index) {
                        if(left_set.find(global_rhs[j]) /*!= left_set.end()*/) {
                            RB[i] = global_rhs[j];

                            // Remove the pushed element
                            rhs_expansion.pop_back();
                        }
                        else if(right_set.find(global_rhs[j]) /*!= right_set.end()*/) {
                            RB[i] = 0;
                        }
                    }

                }

                // RHS Non-Terminal
                else {

                    if(j==start_index) {
                        LB[i] = LB[global_rhs[j]];

                    }
                    else if(j == end_index) {
                        RB[i] = RB[global_rhs[j]];
                    }
                    // LB --> not equal to 0
                    // If j == 0 we set it to LB in slg_nonterm but not in to the expansion.
                    if((j != start_index) and LB[global_rhs[j]]) {
                        rhs_expansion.push_back(LB[global_rhs[j]]);
                    }

                    c_size_t rhs_symbol_start_index = new_slg_nonterm_vec[global_rhs[j]].start_index;
                    c_size_t rhs_symbol_end_index = (global_rhs[j] == (c_size_t)new_slg_nonterm_vec.size()-1) ? new_rhs.size()-1 : new_slg_nonterm_vec[global_rhs[j]+1].start_index - 1;

                    // Check whether Cap is Empty in new_slg_nonterm_vec.
                    if(rhs_symbol_start_index <= rhs_symbol_end_index){
                        rhs_expansion.push_back(global_rhs[j]);
                    }



                    // RB --> not equal to 0.
                    // If j is last element, we set it to RB in slg_nonterm but not in to the expansion
                    if((j != end_index) and RB[global_rhs[j]]) {
                        rhs_expansion.push_back(RB[global_rhs[j]]);
                    }
                }
            }

            // cout << '\n';

            // if(i==5) {
            //     cout << global_rhs[start_index] << ' ' << global_rhs[start_index+1]  << endl;
            //     cout << "RHS Expansion: " << endl;
            //     for(c_size_t x : rhs_expansion) {
            //         cout << x << ' ';
            //     }
            //     cout << endl;
            // }

            // To handle last character/corner case
            // if(rhs_expansion.size()>=2)
            rhs_expansion.push_back(rhs_expansion.back());

            // vector<c_size_t> cap_rhs;

            for(c_size_t j=0; j<(c_size_t)rhs_expansion.size()-1; j++) {
                if(rhs_expansion[j] < 0 && rhs_expansion[j+1] < 0) {
                    if(left_set.find(rhs_expansion[j]) /*!= left_set.end()*/ && right_set.find(rhs_expansion[j+1]) /*!= right_set.end()*/) {
                        // if(m.find({rhs_expansion[j], rhs_expansion[j+1]}) == m.end()) {
                        if(!m.find(make_pair(rhs_expansion[j], rhs_expansion[j+1]))) {
                            // m[{rhs_expansion[j], rhs_expansion[j+1]}] = recompression_rlslp->nonterm.size();
                            m.insert(make_pair(rhs_expansion[j], rhs_expansion[j+1]), recompression_rlslp->nonterm.size());
                            recompression_rlslp->nonterm.emplace_back('1', abs(rhs_expansion[j]), abs(rhs_expansion[j+1]));
                        }

                        new_rhs.push_back(-m[{rhs_expansion[j], rhs_expansion[j+1]}]);
                        j++;
                    }
                    else {
                        new_rhs.push_back(rhs_expansion[j]);
                    }
                }
                else if(rhs_expansion[j] >= 0) {
                    new_rhs.push_back(rhs_expansion[j]);
                }
                else {
                    new_rhs.push_back(rhs_expansion[j]);
                }
            }

	        // cap_rhs.shrink_to_fit();
            new_slg_nonterm_vec.push_back(curr_new_rhs_size);
        }
        else if(curr_rhs_size == 1) {

            // cout << global_rhs[start_index] << '\n';

            // Terminal;
            // A --> a
            // B --> b
            if(global_rhs[start_index] < 0) {
                if(left_set.find(global_rhs[start_index]) /*!= left_set.end()*/) {
                    LB[i] = 0;
                    new_slg_nonterm_vec.push_back(curr_new_rhs_size);
                    RB[i] = global_rhs[start_index];
                }
                else if(right_set.find(global_rhs[start_index]) /*!= right_set.end()*/) {
                    LB[i] = global_rhs[start_index];
                    new_slg_nonterm_vec.push_back(curr_new_rhs_size);
                    RB[i] = 0;
                }
                else {
                    // Non-reachable Non-Terminals.
                    // A --> a
                    // 'a' is not found in adjacency list so push empty.
                    new_slg_nonterm_vec.push_back(curr_new_rhs_size);
                    cerr << "Error: Not Found in Left and Right." << endl;
                }
            }
            else {

                // A --> B (LB(B) . B^ . RB(B))
                // LB(A) = LB(B)
                // RB(A) = RB(B)
                // A^ = B^ (only if B^ is not empty)
                // vector<c_size_t> cap_rhs;

                c_size_t rhs_symbol = global_rhs[start_index];

                LB[i] = LB[rhs_symbol];

                c_size_t rhs_symbol_start_index = new_slg_nonterm_vec[rhs_symbol].start_index;
                c_size_t rhs_symbol_end_index = (rhs_symbol == (c_size_t)new_slg_nonterm_vec.size()-1) ? new_rhs.size()-1 : new_slg_nonterm_vec[rhs_symbol+1].start_index - 1;


                if(rhs_symbol_start_index <= rhs_symbol_end_index) {
                    new_rhs.push_back(rhs_symbol);
                }

                RB[i] = RB[global_rhs[start_index]];
		        // new_rhs.shrink_to_fit();
                new_slg_nonterm_vec.push_back(curr_new_rhs_size);
            }
            
        }
    }

    // arr[0].clear();
    // arr[1].clear();
    // set<c_size_t>().swap(arr[0]);
    // set<c_size_t>().swap(arr[1]);

    delete partition.first;
    delete partition.second;

    const c_size_t &start_var = slg_nonterm_vec.size()-1;

    const c_size_t &start_var_LB = LB[start_var];
    const c_size_t &start_var_RB = RB[start_var];


    // vector<c_size_t> new_start_rhs;
    c_size_t curr_new_rhs_size = new_rhs.size();
    

    if(start_var_LB != 0) {
        new_rhs.push_back(start_var_LB);
    }

    c_size_t rhs_symbol_start_index = new_slg_nonterm_vec[start_var].start_index;
    c_size_t rhs_symbol_end_index = (start_var == (c_size_t)new_slg_nonterm_vec.size()-1) ? curr_new_rhs_size-1 : new_slg_nonterm_vec[start_var+1].start_index - 1;

    if(rhs_symbol_start_index <= rhs_symbol_end_index) {
        new_rhs.push_back(start_var);
    }
    
    if(start_var_RB != 0) {
        new_rhs.push_back(start_var_RB);
    }

    new_slg_nonterm_vec.push_back(curr_new_rhs_size);

    delete slg;

    // return new_slg;
    return new SLG(new_slg_nonterm_vec, new_rhs);
}

RecompressionRLSLP* recompression_on_slp(InputSLP* s) {

    #ifdef TEST_SAMPLE
    delete s;

    s = new InputSLP();
    // s->nonterm.push_back(SLPNonterm('0', 97, 0));
    s->nonterm.push_back(SLPNonterm('0', 98, 0));
    s->nonterm.push_back(SLPNonterm('0', 99, 0));
    s->nonterm.push_back(SLPNonterm('1', 1, 0));
    s->nonterm.push_back(SLPNonterm('1', 2, 0));
    s->nonterm.push_back(SLPNonterm('1', 3, 1));
    s->nonterm.push_back(SLPNonterm('1', 4, 3));
    s->nonterm.push_back(SLPNonterm('1', 5, 4));
    #endif


    SLG* slg = new SLG();

    RecompressionRLSLP* recompression_rlslp = new RecompressionRLSLP();

    // For 0.
    recompression_rlslp->nonterm.emplace_back();

    // Compute S0 from S and Initialize Recompression
    for(c_size_t i=0; i<s->nonterm.size(); i++) {
        const char_t &type = s->nonterm[i].type;
        const c_size_t &first = s->nonterm[i].first;
        const c_size_t &second = s->nonterm[i].second;

        c_size_t global_rhs_size = slg->rhs.size();

        if(type == '0') {
            // For type '0', in SLG, terminals start from -1.
            // Initially, -1 = recompression_rlslp->nonterm size.
            slg->rhs.push_back(-(c_size_t)(recompression_rlslp->nonterm).size());
            // Here first is the ASCII values of the character.
            recompression_rlslp->nonterm.emplace_back('0', first, second);
        }
        else {
            slg->rhs.push_back(first);
            slg->rhs.push_back(second);
        }

        slg->nonterm.push_back(global_rhs_size);

    }

    delete s;

    c_size_t i = 0;

    double max_BComp_time = 0;
    double max_PComp_time = 0;
    while(++i) {
        // map<pair<c_size_t, c_size_t>, c_size_t> m;
        hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> m;

        const vector<SLGNonterm> &slg_nonterm_vec = slg->nonterm;

        const vector<c_size_t> &global_rhs = slg->rhs;

        const SLGNonterm &slg_nonterm = slg_nonterm_vec.back();

        c_size_t start_var = slg_nonterm_vec.size()-1;

        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (start_var == slg_nonterm_vec.size()-1) ? global_rhs.size()-1 : slg_nonterm_vec[start_var+1].start_index - 1;

        c_size_t start_var_rhs_size = end_index - start_index + 1;

        if(start_var_rhs_size == 1 && global_rhs[end_index] < 0) {
            break;
        }

        if(i&1) {
            // cout << i << " BComp \n" << flush;
            auto start_time = std::chrono::high_resolution_clock::now();

            slg = BComp(slg, recompression_rlslp, m);

            auto end_time = std::chrono::high_resolution_clock::now();

            auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            max_BComp_time = max(max_BComp_time, duration_seconds);
        }
        else {
            // cout << i << " PComp \n" << flush;
            auto start_time = std::chrono::high_resolution_clock::now();

            slg = PComp(slg, recompression_rlslp, m);

            auto end_time = std::chrono::high_resolution_clock::now();

            auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            max_PComp_time = max(max_PComp_time, duration_seconds);
        }

        // m.clear();
        m.reset();
        // printRecompressionRLSLP(recompression_rlslp);
    }

    cout << "Runs: " << i << endl << endl;
    cout << "Max. BComp Time: " << max_BComp_time << endl;
    cout << "Max. PComp TIme: " << max_PComp_time << endl;

    delete slg;

    return recompression_rlslp;
}

c_size_t computeExplen(const c_size_t i, vector<RLSLPNonterm>& rlslp_nonterm_vec) {
    // Check if already computed
    if (rlslp_nonterm_vec[i].explen != 0) {
        return rlslp_nonterm_vec[i].explen;
    }

    switch (rlslp_nonterm_vec[i].type) {
        case '0':  // Terminal case
            return rlslp_nonterm_vec[i].explen = 1;

        case '1':  // Binary non-terminal
            return rlslp_nonterm_vec[i].explen = computeExplen(rlslp_nonterm_vec[i].first, rlslp_nonterm_vec) +
                                                 computeExplen(rlslp_nonterm_vec[i].second, rlslp_nonterm_vec);

        default:  // Repetition case or other types
            return rlslp_nonterm_vec[i].explen = rlslp_nonterm_vec[i].second *
                                                 computeExplen(rlslp_nonterm_vec[i].first, rlslp_nonterm_vec);
    }

    // In case of unexpected type
    return 0; 
}

InputSLP* getSLP(c_size_t grammar_size) {

    vector<SLPNonterm> nonterm;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    srand(seed);



    // Generate random number X < n/2
    c_size_t num_terminals = rand() % ((c_size_t)((3*grammar_size)/4));

    num_terminals++;

    for(c_size_t i=1; i<=num_terminals; i++) {
        nonterm.push_back(SLPNonterm('0', -i, 1));
    }


    for (c_size_t i = num_terminals; i < grammar_size; ++i) {
        c_size_t num1 = rand() % i; // Generate first number less than i
        c_size_t num2 = rand() % i; // Generate second number less than i


        c_size_t bit1 = rand()%2;
        c_size_t bit2 = rand()%2;

        c_size_t num = (bit2 == 1 ? num2 : num1);

        if(bit1 == 0)
        nonterm.push_back(SLPNonterm('1', i-1, num));
        else {
            nonterm.push_back(SLPNonterm('1', num, i-1));
        }
    }

    InputSLP* slp = new InputSLP(nonterm);

    return slp;
}

vector<pair<c_size_t, c_size_t>> get_random_queries(c_size_t text_size) {
    // Create a vector to store the pairs
    vector<pair<c_size_t, c_size_t>> pairs;
    pairs.reserve(1000000); // Reserve space for 4096 pairs to avoid reallocation

    // Random number generation setup
    random_device rd;  // Obtain a random number from hardware
    mt19937 gen(rd()); // Seed the generator
    uniform_int_distribution<> distrib(0, text_size - 1);

    // Generate 4096 pairs
    for (c_size_t i = 0; i < 1000000; ++i) {
        c_size_t first = distrib(gen);
        c_size_t second = distrib(gen);
        pairs.emplace_back(first, second); // Emplace pair into the vector
    }

    return pairs;
}


void test(c_size_t text_size, RecompressionRLSLP *recompression_rlslp, vector<RLSLPNonterm> & rlslp_nonterm_vec, const string &raw_input_text) {

    vector<c_size_t> arr = expandRLSLP(recompression_rlslp);

    ifstream file(raw_input_text);

    if (!file.is_open()) {
        cerr << "Error opening file - " + raw_input_text << endl;
        exit(1);
    }

    vector<c_size_t> text;

    char ch;
    while (file.get(ch)) {
        text.push_back(static_cast<unsigned char>(ch));
    }

    // Close the file
    file.close();

    assert(text.size()==arr.size());

    cout << "Text Size Matched!" << endl;

    for(c_size_t i=0; i<text.size(); i++) {
        if(text[i] != arr[i]) {
            cout << "Text didn't match : " << i << ' ' << text[i] << ' ' << arr[i] << endl;
            exit(1);
        }
    }

    cout << "Text Matched!" << endl;

    vector<pair<c_size_t, c_size_t>> random_queries = get_random_queries(text_size);

    auto start_time = std::chrono::high_resolution_clock::now();

    for(auto x : random_queries) {

        c_size_t i = x.first;
        c_size_t j = x.second;
        // if(i!=7 or j!=11) continue;

        // cout << i << ' ' << j << endl;

    
        Node v1, v2;
        stack<Node> v1_ancestors, v2_ancestors;
        // v1_ancestors.push(Node(grammar.size()-1, 0, 33));
        // v2_ancestors.push(Node(grammar.size()-1, 0, 33));
        initialize_nodes(rlslp_nonterm_vec.size() - 1, i, 0, rlslp_nonterm_vec.back().explen - 1, v1_ancestors, rlslp_nonterm_vec, v1);
        initialize_nodes(rlslp_nonterm_vec.size() - 1, j, 0, rlslp_nonterm_vec.back().explen - 1, v2_ancestors, rlslp_nonterm_vec, v2);


        // cout << v1.var << ' ' << v1.l << ' ' << v1.r << endl;
        // cout << v2.var << ' ' << v2.l << ' ' << v2.r << endl;

        // cout << i << ' ' << j << endl;

        c_size_t res1 = LCE(v1, v2, i, v1_ancestors, v2_ancestors, rlslp_nonterm_vec);

        c_size_t res2 = 0;

        c_size_t ii = i;
        c_size_t jj = j;

        while(jj < text_size && arr[ii] == arr[jj]) {
            res2++;

            jj++;
            ii++;
        }

        // cout << i << ' ' << j << ' ' << res1 << ' ' << res2 << endl;
        // assert(res1==res2);

        if(res1 != res2) {
            cerr << "Error: LCE Query didn't match with Naive!!" << endl;
            exit(1);
        }


    }

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    cout << "LCE Queries Passed!" << endl;
    // Output the duration
    cout << "Time taken for LCE Queries: " << duration_seconds << " seconds" << endl;
}

void start_compression(const string &input_file, const string &raw_input_text) {
    InputSLP *inputSLP = new InputSLP();
    inputSLP->read_from_file(input_file);

    c_size_t j = 0;
    cout << "Number of Non-Terminals : " << inputSLP->nonterm.size() << endl;

    c_size_t i = -1;

    auto start_time = std::chrono::high_resolution_clock::now();

    RecompressionRLSLP *recompression_rlslp = recompression_on_slp(inputSLP);

    vector<RLSLPNonterm> & rlslp_nonterm_vec = recompression_rlslp->nonterm;

    for(c_size_t i=rlslp_nonterm_vec.size()-1; i>=1; i--) {
        computeExplen(i, rlslp_nonterm_vec);
    }

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Output the duration
    cout << "Time taken for Construction: " << duration_seconds << " seconds" << endl;

    // // printRecompressionRLSLP(recompression_rlslp);
    c_size_t text_size = rlslp_nonterm_vec.back().explen;

    cout << "Text Size : " << text_size << endl;

    // TEST
    if(!raw_input_text.empty())
    test(text_size, recompression_rlslp, rlslp_nonterm_vec, raw_input_text);
}

int main(int argc, char *argv[]) {

    if(!(argc == 2 or argc == 4)) {
        cout << "Usage: " + string(argv[0]) + " [-t raw_input_text] input_slp" << endl;
        exit(1);
    }

    string raw_input_text;
    string input_slp;

    for(c_size_t i=1; i<argc; i++) {
        string curr_arg(argv[i]);

        if(curr_arg == "-t") {
            if(i+1 == argc) {
                cout << "Usage: " + string(argv[0]) + " [-t raw_input_text] input_slp" << endl;
                exit(1);
            }

            raw_input_text = string(argv[i+1]);
            i++;
        }
        else {
            input_slp = string(argv[i]);
        }
    }

    if(!raw_input_text.empty()) {
        ifstream file(raw_input_text);

        if (!file.is_open()) {
            cerr << "Error: Unable to open the file - " + raw_input_text << endl;
            exit(1);
        }

        file.close();
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    start_compression(input_slp, raw_input_text);

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Output the duration
    cout << "Total Time taken: " << duration_seconds << " seconds" << endl;

    return 0;
}

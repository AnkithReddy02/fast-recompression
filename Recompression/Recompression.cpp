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
#include "utils.hpp"

// #define TEST_SAMPLE
// #define DEBUG_LOG

using namespace std;
using namespace utils;

template<>
std::uint64_t get_hash(const packed_pair<c_size_t, c_size_t> &x) {
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
std::string keyToString(const packed_pair<c_size_t,c_size_t> &x) {
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

void combineFrequenciesInRange(
    const space_efficient_vector<packed_pair<c_size_t, c_size_t>>& vec,
    const c_size_t &lr_pointer,
    const c_size_t &rr_pointer,
    space_efficient_vector<SLGNonterm> &new_slg_nonterm_vec,
    /* map<pair<c_size_t, c_size_t>, c_size_t> &m */
    hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m,
    RecompressionRLSLP *recompression_rlslp,
    space_efficient_vector<c_size_t> &new_rhs) {

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

    for(size_t i = lr_pointer + 1; i <= rr_pointer && i < vec.size(); ++i) {
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
                        recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', abs(p.first), p.second));
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
            recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', abs(p.first), p.second));
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
c_size_t lower_bound(const space_efficient_vector<T> &vec, const T &search_element) {
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
    hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t> & m) {

    // Current slg non-term list
    space_efficient_vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

    // Global RHS
    space_efficient_vector<c_size_t> &global_rhs = slg->rhs;

    // // New SLG that needs to be created by applying BComp
    SLG *new_slg = new SLG();

    // New SLG non-term list that new_slg needs.
    space_efficient_vector<SLGNonterm> &new_slg_nonterm_vec = new_slg->nonterm;
    new_slg_nonterm_vec.reserve(slg_nonterm_vec.size() + 1);
    space_efficient_vector<c_size_t> &new_rhs = new_slg->rhs;

    const c_size_t &grammar_size = slg_nonterm_vec.size();

    typedef packed_pair<c_size_t, char_t> small_pair_type;
    space_efficient_vector<small_pair_type> LR_vec(grammar_size, small_pair_type(-1, 0));
    space_efficient_vector<small_pair_type> RR_vec(grammar_size, small_pair_type(-1, 0));

    typedef packed_pair<c_size_t, c_size_t> pair_type;
    space_efficient_vector<pair_type> large_LR_vec;
    space_efficient_vector<pair_type> large_RR_vec;

    space_efficient_vector<pair_type> rhs_expansion;

    // We shall iterate throught each production rule in the increasing order of variable.
    // 'i' --> represents the variable.
    for(c_size_t i = 0; i < grammar_size; i++) {
        #ifdef DEBUG_LOG
        if (!(i & ((1 << 19) - 1)))
          cout << fixed << setprecision(2) << "Progress: " << ((i+1)*(long double)100)/grammar_size << '\r' << flush;
        #endif

        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (i == grammar_size-1) ? c_size_t(global_rhs.size())-1 : c_size_t(slg_nonterm_vec[i+1].start_index) - 1;

        if(start_index > end_index) {
            new_slg_nonterm_vec.push_back((c_size_t)new_rhs.size());
            continue;
        }

        // Create the expansion of RHS.
        // space_efficient_vector<pair<c_size_t, c_size_t>> rhs_expansion;
        rhs_expansion.set_empty();
     
        // Compute the Expansion.
        for(c_size_t j = start_index; j <= end_index; j++) {

            const c_size_t &rhs_symbol = global_rhs[j];

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

                    c_size_t search_index = lower_bound(large_LR_vec, search_element);

                    assert(large_LR_vec[search_index].first == rhs_symbol);
                    assert(search_index != large_LR_vec.size());

                    k = large_LR_vec[search_index].second;
                }

                rhs_expansion.push_back(pair_type((c_size_t)LR_vec[rhs_symbol].first, (c_size_t)k));
            }

            c_size_t rhs_symbol_start_index = new_slg_nonterm_vec[rhs_symbol].start_index;
            c_size_t rhs_symbol_end_index = (rhs_symbol == (c_size_t)new_slg_nonterm_vec.size()-1) ? (c_size_t)new_rhs.size()-1 : new_slg_nonterm_vec[rhs_symbol+1].start_index - 1;

           
            // Cap is not empty --> in new SLG the variable(rhs_symbol) RHS is not empty --> then Cap is not empty.
            if(rhs_symbol_start_index <= rhs_symbol_end_index) {
                rhs_expansion.push_back(pair_type(rhs_symbol, 0));
            }

            // **RR** of a current variable(rhs_symbol) in current RHS is **not empty**.
            if(static_cast<c_size_t>(RR_vec[rhs_symbol].second) != 0) {
                c_size_t k = static_cast<c_size_t>(RR_vec[rhs_symbol].second);
                if(k == 255) {
                    pair_type search_element = pair_type(rhs_symbol, 0);

                    c_size_t search_index = lower_bound(large_RR_vec, search_element);

                    assert(large_RR_vec[search_index].first == rhs_symbol);
                    assert(search_index != large_RR_vec.size());

                    k = large_RR_vec[search_index].second;
                }

                rhs_expansion.push_back(pair_type((c_size_t)RR_vec[rhs_symbol].first, (c_size_t)k));
            }
        }

        //rhs.clear();
        //vector<int>().swap(rhs);

        // Compute LR
        c_size_t lr_pointer = 1;

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
        c_size_t rr_pointer = rhs_expansion.size()-2;

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

                rr_pointer--;
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
            combineFrequenciesInRange(rhs_expansion, lr_pointer, rr_pointer, new_slg_nonterm_vec, m, recompression_rlslp, new_rhs);
        }
    }
    #ifdef DEBUG_LOG
    cout << endl;
    #endif

    // Add new Starting Variable to the new Grammar G'(G Prime).
    const c_size_t &start_var = slg_nonterm_vec.size()-1;

    packed_pair<c_size_t, c_size_t> start_var_LR;
    if(static_cast<c_size_t>(LR_vec[start_var].second) == 255) {

        c_size_t k = static_cast<c_size_t>(LR_vec[start_var].second);
        pair_type search_arr((c_size_t)start_var, (c_size_t)0);
        
        c_size_t search_index = lower_bound(large_LR_vec, search_arr);

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
        
        c_size_t search_index = lower_bound(large_RR_vec, search_arr);

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
                recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', abs(start_var_LR.first), start_var_LR.second));
            }
            new_rhs.push_back(-m[start_var_LR]);
        }
        else {
            new_rhs.push_back(start_var_LR.first);
        }
    }

    c_size_t start_var_start_index = new_slg_nonterm_vec[start_var].start_index;
    c_size_t start_var_end_index = (start_var == (c_size_t)new_slg_nonterm_vec.size()-1) ? (c_size_t)curr_new_rhs_size-1 : new_slg_nonterm_vec[start_var+1].start_index - 1;

    if(start_var_start_index <= start_var_end_index) {
        new_rhs.push_back(start_var);
    }

    if(start_var_RR.second != 0) {
        if(start_var_RR.second >= 2) {
            // if(m.find(start_var_RR) == m.end()) {
            if(!m.find(start_var_RR)) {
                // m[start_var_RR] = recompression_rlslp->nonterm.size();
                m.insert(start_var_RR, recompression_rlslp->nonterm.size());
                recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', abs(start_var_RR.first), start_var_RR.second));
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

c_size_t computeVOccHelper(space_efficient_vector<packed_pair<c_size_t, c_size_t>> & edges, space_efficient_vector<c_size_t> &curr_index, space_efficient_vector<bool_t> &have_edges, c_size_t u, space_efficient_vector<c_size_t> & dp) {

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

void computeVOcc(SLG *slg, space_efficient_vector<c_size_t> &dp) {

    // Here, any list size is slg_nonterm_vec.size()
    space_efficient_vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;
    space_efficient_vector<c_size_t> &global_rhs = slg->rhs;

    // Create Reverse Graph.
    // vector<vector<pair<c_size_t,c_size_t>>> graph(slg_nonterm_vec.size(), vector<pair<c_size_t, c_size_t>>());
    space_efficient_vector<c_size_t> curr_index(slg_nonterm_vec.size(), -1);
    // space_efficient_vector<pair<c_size_t, c_size_t>> edges;

    hash_table<c_size_t, bool> unique_var;

    for(c_size_t i=0; i<slg_nonterm_vec.size(); i++) {

        // Current SLGNonterm
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (i == slg_nonterm_vec.size()-1) ? (c_size_t)global_rhs.size()-1 : slg_nonterm_vec[i+1].start_index - 1;

        c_size_t curr_rhs_size = end_index - start_index + 1;

        // set<c_size_t> unique_var;
        // hash_table<c_size_t, bool> unique_var;
        unique_var.reset();

        // Enumerate each character of RHS
        // Frequency calculation for reverse of the graph.
        for(c_size_t j=start_index; j<=end_index; j++) {
            const c_size_t & var = global_rhs[j];
            // Only Non-Terminals
            if(var >= 0) {
                // unique_var.insert(var);
                unique_var.insert(var, true);
            }
        }

        // Construct Edges.
        // for(const auto & v : unique_var) {
        for(c_size_t i = 0; i < unique_var.size(); ++i) {
            c_size_t v = unique_var.get(i).first;
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

    // edges.resize(prefix_sum);
    typedef packed_pair<c_size_t, c_size_t> pair_type;
    space_efficient_vector<pair_type> edges(prefix_sum);

    hash_table<c_size_t, c_size_t, c_size_t> var_freq;

    // Enumerate Each Production Rule.
    for(c_size_t i=0; i<slg_nonterm_vec.size(); i++) {
        // Current SLGNonterm
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (i == slg_nonterm_vec.size()-1) ? (c_size_t)global_rhs.size()-1 : slg_nonterm_vec[i+1].start_index - 1;

        c_size_t curr_rhs_size = end_index - start_index + 1;

        // unordered_map<c_size_t, c_size_t> var_freq;
        // hash_table<c_size_t, c_size_t, c_size_t> var_freq;
        var_freq.reset();

        // Enumerate each character of RHS
        // Frequency calculation for reverse of the graph.
        for(c_size_t j=start_index; j<=end_index; j++) {

            const c_size_t & var = global_rhs[j];
            // Only Non-Terminals
            if(var >= 0) {
                // var_freq[var]++;
                /*if(!var_freq.find(var)) {
                    var_freq.insert(var, 0);
                }
                var_freq.insert(var, var_freq[var] + 1);*/
                c_size_t *freq = var_freq.find(var);
                if (!freq) var_freq.insert(var, 1);
                else *freq += 1;
            }
        }

        // Construct Edges.
        // for(const auto & x : var_freq) {
        for(c_size_t j = 0; j < var_freq.size(); ++j) {
            auto x = var_freq.get(j); 

            // edge from v to u, u <-- v
            const c_size_t u = i;
            const c_size_t v = x.first;
            const c_size_t weight = x.second;

            // Weighted Edge from v to u , v --> u
            // graph[v].emplace_back(u, weight);
            edges[--curr_index[v]] = pair_type(u, weight);
            // if(curr_index[v]==-1) curr_index[v] = 0;
            // curr_index[v]++;
            
        }
    }

    space_efficient_vector<bool_t> have_edges(curr_index.size(), false);

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
    // dp.resize(slg_nonterm_vec.size(), -1);

    // Compute vOcc.
    for(c_size_t i = slg_nonterm_vec.size() - 1; i >= 0; i--) {
        #ifdef DEBUG_LOG
        if (!(i & ((1 << 19) - 1)))
          cout << fixed << setprecision(2) << "  VOcc Progress: " << (slg_nonterm_vec.size() - i) * (double)100/(slg_nonterm_vec.size()) << '\r';
        #endif
        computeVOccHelper(edges, curr_index, have_edges, i, dp);
    }

    #ifdef DEBUG_LOG
        cout << endl;
    #endif

    return;
}

bool AdjListElementCompare(const AdjListElement &a, const AdjListElement &b) {
    if(a.first != b.first) {
        return a.first < b.first;
    }

    return a.second < b.second;
}

void sortAdjList(space_efficient_vector<AdjListElement> & adjList) {

    // Make Positive
    // for(array<int, 4> & a : adjList) {
    //     a[0] = -a[0];
    //     a[1] = -a[1];
    // }

    // for(AdjListElement & a : adjList) {
    //     a.first = -a.first;
    //     a.second = -a.second;
    // }

    for(c_size_t i = 0; i < adjList.size(); ++i) {
        adjList[i].first = -adjList[i].first;
        adjList[i].second = -adjList[i].second;
    }

    // Sort It.
    //radixSort(adjList);
    adjList.sort();
    // sort(adjList.begin(), adjList.end(), AdjListElementCompare);


    // Revert to Negative
    // for(AdjListElement & a : adjList) {
    //     a.first = -a.first;
    //     a.second = -a.second;
    // }

    for(c_size_t i = 0; i < adjList.size(); ++i) {
        adjList[i].first = -adjList[i].first;
        adjList[i].second = -adjList[i].second;
    }
}

// MAP USAGE
packed_pair<c_size_t, c_size_t> computeAdjListHelper(
    c_size_t var,
    SLG *slg,
    //map<pair<c_size_t, c_size_t>, c_size_t> &m0, map<pair<c_size_t, c_size_t>, c_size_t> &m1, vector<pair<c_size_t, c_size_t>> & dp, vector<c_size_t> &vOcc) {
    hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m0,
    hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m1,
    space_efficient_vector<packed_pair<c_size_t, c_size_t>> &dp,
    space_efficient_vector<c_size_t> &vOcc) {

    typedef packed_pair<c_size_t, c_size_t> pair_type;

    if(var < 0) {
        return pair_type(var, var);
    }

    if(dp[var].first != 1 && dp[var].second != 1) {
        return dp[var];
    }

    /*
        1. Modifying LMS and RMS
        2. Access rhs
    */
    SLGNonterm &slg_nonterm = slg->nonterm[var];

    const space_efficient_vector<c_size_t> & global_rhs = slg->rhs;
    space_efficient_vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

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
        return dp[var] = pair_type(global_rhs[start_index], global_rhs[start_index]);

        // return dp[var] = {rhs[0], rhs[0]};
    }

    pair_type first_lms_rms;
    pair_type last_lms_rms;
    pair_type prev_lms_rms;

    for(c_size_t j = start_index; j <= end_index; ++j) {
        const c_size_t &rhs_symbol = global_rhs[j];
        pair_type lms_rms = computeAdjListHelper(rhs_symbol, slg, m0, m1, dp, vOcc);

        if(j == start_index) {
            first_lms_rms = lms_rms;
        }

        if(j == end_index) {
            last_lms_rms = lms_rms;
        }

        if(j > start_index) {
            c_size_t f = prev_lms_rms.second;
            c_size_t s = lms_rms.first;

            bool_t swapped = false;

            if(abs(f) < abs(s)) {
                swap(f, s);
                swapped = true;
            }

            if(swapped) {
                pair_type p(f, s);

                c_size_t *val = m1.find(p);
                if(!val) {
                    m1.insert(p, vOcc[var]);
                }
                else {
                    *val += vOcc[var];
                }
            }
            else {
                pair_type p(f, s);

                c_size_t *val = m0.find(p);
                if(!val) {
                    m0.insert(p, vOcc[var]);
                }
                else {
                    *val += vOcc[var];
                }
            }
        }

        prev_lms_rms = lms_rms;
    }

    return dp[var] = pair_type(first_lms_rms.first, last_lms_rms.second);

    /*
    space_efficient_vector<pair_type> lms_rms_list;

    

    for(c_size_t j = start_index; j <= end_index; j++) {
        const c_size_t &rhs_symbol = global_rhs[j];
        pair_type lms_rms = computeAdjListHelper(rhs_symbol, slg, m0, m1, dp, vOcc);
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
            pair_type p(f, s);

            if(!m1.find(p)) {
                m1.insert(p, 0);
            }

            m1.insert(p, m1[p] + vOcc[var]);
        }
        else {
            // m0[{f, s}] += vOcc[var];
            pair_type p(f, s);

            if(!m0.find(p)) {
                m0.insert(p, 0);
            }

            m0.insert(p, m0[p] + vOcc[var]);
        }
        // adjList.push_back(AdjListElement(f, s, swapped, vOcc[var]));
    }

    //slg_nonterm.LMS = lms_rms_list.front().first;
    //slg_nonterm.RMS = lms_rms_list.back().second;

    return dp[var] = pair_type(lms_rms_list.front().first, lms_rms_list.back().second);
    */


}

void computeAdjList(
    SLG *slg,
    //map<pair<c_size_t, c_size_t>, c_size_t> &m0, map<pair<c_size_t, c_size_t>, c_size_t> &m1, vector<c_size_t> &vOcc) {
    hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m0,
    hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m1,
    space_efficient_vector<c_size_t> &vOcc) {

    typedef packed_pair<c_size_t, c_size_t> pair_type;
    space_efficient_vector<pair_type> dp(slg->nonterm.size(), pair_type((c_size_t)1, (c_size_t)1));

    computeAdjListHelper(slg->nonterm.size()-1, slg, m0, m1, dp, vOcc);
    return;
}

// void createPartition(const vector<AdjListElement> & adjList, array<set<c_size_t>, 2> &partition_set) {
std::pair<hash_table<c_size_t, bool> *, hash_table<c_size_t, bool> *>
createPartition(const space_efficient_vector<AdjListElement> & adjList) {
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

    // for(const AdjListElement & arr : adjList) {
    for(c_size_t i = 0; i < adjList.size(); ++i) {
        const AdjListElement & arr = adjList[i];

        c_size_t f = arr.first;
        c_size_t s = arr.second;

        if(arr.swapped == true) {
            swap(f, s);
        }

        LRPairsCount += ((leftSet.find(f) /*!= leftSet.end()*/) && (rightSet.find(s) /*!= rightSet.end()*/)) ? arr.vOcc : (c_size_t)0;
        RLPairsCount += ((rightSet.find(f) /*!= rightSet.end()*/) && (leftSet.find(s) /*!= leftSet.end()*/)) ? arr.vOcc : (c_size_t)0;
    }

    if (RLPairsCount < LRPairsCount) {
        // swap(leftSet, rightSet);
        swap(leftSet_ptr, rightSet_ptr);
    }

    // swap(rightSet, leftSet);
    swap(rightSet_ptr, leftSet_ptr);

    return make_pair(leftSet_ptr, rightSet_ptr);
}

std::pair<hash_table<c_size_t, bool> *, hash_table<c_size_t, bool> *>
createPartition(SLG *slg) {
    space_efficient_vector<c_size_t> vOcc(slg->nonterm.size(), -1);

    #ifdef DEBUG_LOG
    cout << "Computing VOcc..." << endl;
    #endif

    // Compute vOcc
    computeVOcc(slg, vOcc);

    #ifdef DEBUG_LOG
    cout << "Computing AdjList..." << endl;
    #endif

    // Compute AdjList
    typedef packed_pair<c_size_t, c_size_t> pair_type;
    hash_table<pair_type, c_size_t, c_size_t> m0;
    hash_table<pair_type, c_size_t, c_size_t> m1;
    computeAdjList(slg, m0, m1, vOcc);
    // Avoid resizing.
    space_efficient_vector<AdjListElement> adjList(m0.size() + m1.size());

    assert(adjList.size() != 0);

    c_size_t index = 0;
    // for(auto & x : m0) {
    //     adjList[index++] = AdjListElement(x.first.first, x.first.second, false, x.second);
    // }
    for(c_size_t i = 0; i < m0.size(); i++) {
        std::pair<pair_type, c_size_t> x = m0.get(i);
        adjList[index++] = AdjListElement(x.first.first, x.first.second, false, x.second);
    }

    // m0.clear();

    // for(auto & x : m1) {
    //     adjList[index++] = AdjListElement(x.first.first, x.first.second, true, x.second);
    // }
    for(c_size_t i = 0; i < m1.size(); i++) {
        std::pair<pair_type, c_size_t> x = m1.get(i);
        adjList[index++] = AdjListElement(x.first.first, x.first.second, true, x.second);
    }

    // m1.clear();
    #ifdef DEBUG_LOG
    cout << "Sorting AdjList..." << endl;
    #endif
    sortAdjList(adjList);

    #ifdef DEBUG_LOG
    cout << "Computing Partition..." << endl;
    #endif

    // Create Partition.
    // array<set<c_size_t>,2> arr;
    // createPartition(adjList, arr);
    typedef c_size_t key_type;
    typedef bool value_type;
    typedef hash_table<key_type, value_type> hash_table_type;

    std::pair<hash_table_type*, hash_table_type*> partition = createPartition(adjList);
    
    return partition;
}

// Pair-Wise Compression
SLG * PComp(SLG *slg, RecompressionRLSLP *recompression_rlslp,  //map<pair<c_size_t, c_size_t>, c_size_t> & m) { 
    hash_table<packed_pair<c_size_t, c_size_t>, c_size_t, c_size_t> &m) { 
    
    typedef c_size_t key_type;
    typedef bool value_type;
    typedef hash_table<key_type, value_type> hash_table_type;

    typedef packed_pair<c_size_t, c_size_t> pair_type;

    std::pair<hash_table_type*, hash_table_type*> partition = createPartition(slg);

    #ifdef DEBUG_LOG
    cout << "Performing PComp..." << endl;
    #endif

    // adjList.clear();  // Clear immediately
    // vector<AdjListElement>().swap(adjList);  // Ensure memory is released
    // adjList.clear();

    // const set<c_size_t> &left_set = arr[0], &right_set = arr[1];
    
    const hash_table_type &left_set = *partition.first;
    const hash_table_type &right_set = *partition.second;

    // Current slg non-term list
    space_efficient_vector<SLGNonterm> &slg_nonterm_vec = slg->nonterm;

    // Current rhs
    space_efficient_vector<c_size_t> &global_rhs = slg->rhs;

    // New SLG that needs to be created by applying BComp
    SLG *new_slg = new SLG();

    // New SLG non-term list that new_slg needs.
    space_efficient_vector<SLGNonterm> &new_slg_nonterm_vec = new_slg->nonterm;
    new_slg_nonterm_vec.reserve(slg_nonterm_vec.size() + 1);
    // New RHS
    space_efficient_vector<c_size_t> &new_rhs = new_slg->rhs;

    space_efficient_vector<c_size_t> LB(slg_nonterm_vec.size(), 0), RB(slg_nonterm_vec.size(), 0);

    space_efficient_vector<c_size_t> rhs_expansion;

    // We shall iterate throught each production rule in the increasing order of variable.
    // 'i' --> represents the variable.
    for(c_size_t i=0; i<slg_nonterm_vec.size(); i++) {

        // cout << i << " : ";
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        // const vector<c_size_t> & rhs = slg_nonterm.rhs;
        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (i == slg_nonterm_vec.size()-1) ? (c_size_t)global_rhs.size()-1 : slg_nonterm_vec[i+1].start_index - 1;

        c_size_t curr_rhs_size = end_index - start_index + 1;

        c_size_t curr_new_rhs_size = new_rhs.size();

        if(curr_rhs_size == 0) {
            new_slg_nonterm_vec.push_back((c_size_t)new_rhs.size());
            // cout << '\n';
            continue;
        }
        if((c_size_t)curr_rhs_size >= 2) {

            // space_efficient_vector<c_size_t> rhs_expansion;
            rhs_expansion.set_empty();

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
                    c_size_t rhs_symbol_end_index = (global_rhs[j] == (c_size_t)new_slg_nonterm_vec.size()-1) ? (c_size_t)new_rhs.size()-1 : new_slg_nonterm_vec[global_rhs[j]+1].start_index - 1;

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
                        if(!m.find(pair_type(rhs_expansion[j], rhs_expansion[j+1]))) {
                            // m[{rhs_expansion[j], rhs_expansion[j+1]}] = recompression_rlslp->nonterm.size();
                            m.insert(pair_type(rhs_expansion[j], rhs_expansion[j+1]), recompression_rlslp->nonterm.size());
                            recompression_rlslp->nonterm.push_back(RLSLPNonterm('1', abs(rhs_expansion[j]), abs(rhs_expansion[j+1])));
                        }

                        new_rhs.push_back(-m[pair_type(rhs_expansion[j], rhs_expansion[j+1])]);
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
                c_size_t rhs_symbol_end_index = (rhs_symbol == (c_size_t)new_slg_nonterm_vec.size()-1) ? (c_size_t)new_rhs.size()-1 : new_slg_nonterm_vec[rhs_symbol+1].start_index - 1;


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
    c_size_t rhs_symbol_end_index = (start_var == (c_size_t)new_slg_nonterm_vec.size()-1) ? (c_size_t)curr_new_rhs_size-1 : new_slg_nonterm_vec[start_var+1].start_index - 1;

    if(rhs_symbol_start_index <= rhs_symbol_end_index) {
        new_rhs.push_back(start_var);
    }
    
    if(start_var_RB != 0) {
        new_rhs.push_back(start_var_RB);
    }

    new_slg_nonterm_vec.push_back(curr_new_rhs_size);

    delete slg;

    return new_slg;
    // return new SLG(new_slg_nonterm_vec, new_rhs);
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
    recompression_rlslp->nonterm.push_back(RLSLPNonterm());

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
            recompression_rlslp->nonterm.push_back(RLSLPNonterm('0', first, second));
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

    typedef packed_pair<c_size_t, c_size_t> pair_type;
    hash_table<pair_type, c_size_t, c_size_t> m;

    while(++i) {
        // map<pair<c_size_t, c_size_t>, c_size_t> m;
        // hash_table<pair<c_size_t, c_size_t>, c_size_t, c_size_t> m;

        const space_efficient_vector<SLGNonterm> &slg_nonterm_vec = slg->nonterm;

        const space_efficient_vector<c_size_t> &global_rhs = slg->rhs;

        const SLGNonterm &slg_nonterm = slg_nonterm_vec.back();

        c_size_t start_var = slg_nonterm_vec.size()-1;

        c_size_t start_index = slg_nonterm.start_index;
        c_size_t end_index = (start_var == slg_nonterm_vec.size()-1) ? c_size_t(global_rhs.size())-1 : c_size_t(slg_nonterm_vec[start_var+1].start_index) - 1;

        c_size_t start_var_rhs_size = end_index - start_index + 1;

        if(start_var_rhs_size == 1 && global_rhs[end_index] < 0) {
            break;
        }

        if(i&1) {
            #ifdef DEBUG_LOG
            cout << '\n' << i << " : Starting BComp..." << endl;
            #endif
            auto start_time = std::chrono::high_resolution_clock::now();

            slg = BComp(slg, recompression_rlslp, m);

            auto end_time = std::chrono::high_resolution_clock::now();

            auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            max_BComp_time = max(max_BComp_time, duration_seconds);
        }
        else {
            #ifdef DEBUG_LOG
            cout << '\n' << i << " : Starting PComp..." << endl;
            #endif

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

    cout << endl;
    cout << "Runs: " << i << endl << endl;
    cout << "Max. BComp Time: " << max_BComp_time << " seconds" << endl;
    cout << "Max. PComp TIme: " << max_PComp_time << " seconds" << endl << endl;

    delete slg;

    return recompression_rlslp;
}

c_size_t computeExplen(const c_size_t i, space_efficient_vector<RLSLPNonterm>& rlslp_nonterm_vec) {
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

    space_efficient_vector<SLPNonterm> nonterm;

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

void get_random_queries(c_size_t text_size,
    space_efficient_vector<pair<c_size_t, c_size_t>> &pairs) {
    // pairs.reserve(1000000); 

    random_device rd; 
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(0, text_size - 1);

    // Generate 4096 pairs
    for (c_size_t i = 0; i < 1000000; ++i) {
        c_size_t first = distrib(gen);
        c_size_t second = distrib(gen);
        pairs.push_back(make_pair(first, second));
    }
}

void test(c_size_t text_size, RecompressionRLSLP *recompression_rlslp, space_efficient_vector<RLSLPNonterm> & rlslp_nonterm_vec, const string &raw_input_text) {

    space_efficient_vector<c_size_t> arr;
    expandRLSLP(recompression_rlslp, arr);

    ifstream file(raw_input_text);

    if (!file.is_open()) {
        cerr << "Error opening file - " + raw_input_text << endl;
        exit(1);
    }

    space_efficient_vector<c_size_t> text;

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

    space_efficient_vector<pair<c_size_t, c_size_t>> random_queries;
    get_random_queries(text_size, random_queries);

    auto start_time = std::chrono::high_resolution_clock::now();

    // for(auto x : random_queries) {
    for(c_size_t k = 0; k < random_queries.size(); ++k) {
        const auto& x = random_queries[k];

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

    space_efficient_vector<RLSLPNonterm> & rlslp_nonterm_vec = recompression_rlslp->nonterm;

    for(c_size_t i=rlslp_nonterm_vec.size()-1; i>=1; i--) {
        computeExplen(i, rlslp_nonterm_vec);
    }

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Output the duration
    cout << "Time taken for Construction: " << duration_seconds << " seconds" << endl;

    // // printRecompressionRLSLP(recompression_rlslp);
    c_size_t text_size = rlslp_nonterm_vec.back().explen;

    cout << "Text Size : " << text_size << endl << endl;

    cout << "****************" << endl;
    cout << "Peak RAM usage for Construction: " << get_peak_ram_allocation()/(1024.0 * 1024.0) << " MB" << endl;
    cout << "****************" << endl << endl;

    // TEST
    if(!raw_input_text.empty())
    test(text_size, recompression_rlslp, rlslp_nonterm_vec, raw_input_text);

    // cout << "****************" << endl;
    // cout << "Overall Peak RAM uage: " << get_peak_ram_allocation() << endl;
    // cout << "****************" << endl << endl;

    // Clean up.
    delete recompression_rlslp;

    cout << "****************" << endl;
    cout << "Current RAM usage for Construction: " <<
      get_current_ram_allocation() << " bytes" << endl;
    cout << "****************" << endl << endl;
}

int main(int argc, char *argv[]) {
    utils::initialize_stats();

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
#include<bits/stdc++.h>
#include "recompression_definitions.hpp"
#include "../radix_sort.h"
#include "../io.h"
#include "utilities.hpp"
#include "lce_queries.hpp"
#include "hash_pair.hpp"
using namespace std;

void combineFrequenciesInRange(const vector<pair<int, int>>& vec, const int &lr_pointer, const int &rr_pointer, vector<pair<int, int>> &result) {
    // Check if vector is empty
    if (vec.empty() || lr_pointer > rr_pointer)
        return;

    // Iterate through the vector within the specified range
    int currNum = vec[lr_pointer].first;
    int currFreq = vec[lr_pointer].second;

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

vector<pair<int, int>> combineFrequenciesInRange(const vector<pair<int, int>>& vec, const int &lr_pointer, const int &rr_pointer) {
    vector<pair<int, int>> result;

    // Check if vector is empty
    if (vec.empty() || lr_pointer > rr_pointer)
        return result;

    // Iterate through the vector within the specified range
    int currNum = vec[lr_pointer].first;
    int currFreq = vec[lr_pointer].second;

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

void combineFrequenciesInRange(const vector<pair<int, int>>& vec, const int &lr_pointer, const int &rr_pointer, vector<SLGNonterm> &new_slg_nonterm_vec, unordered_map<pair<int, int>, int, hash_pair> &m, unique_ptr<RecompressionRLSLP> &recompression_rlslp) {
    // Check if vector is empty
    if (vec.empty() || lr_pointer > rr_pointer) {
        new_slg_nonterm_vec.emplace_back();
        return;
    }

    vector<int> cap_rhs;

    // Iterate through the vector within the specified range
    int currNum = vec[lr_pointer].first;
    int currFreq = vec[lr_pointer].second;

    for (size_t i = lr_pointer + 1; i <= rr_pointer && i < vec.size(); ++i) {
        // Check if the current and previous elements have the same number
        if (vec[i].first == currNum) {
            // Merge frequencies
            currFreq += vec[i].second;
        } else {
            
            {
                const pair<int, int> p = {currNum, currFreq};

                // Frequency is >=2 so merge and push to rlslp
                // Merged variable is terminal for a new_slg, so it is marked negative
                if(p.second >= 2) {
                    if(m.find(p) == m.end()) {
                        m[p] = recompression_rlslp->nonterm.size();
                        recompression_rlslp->nonterm.emplace_back('2', abs(p.first), p.second);
                    }
                    
                    //  ** Negative **
                    cap_rhs.push_back(-m[p]);
                }
                else {
                    cap_rhs.push_back(p.first);
                }
            }

            // Move to the next number
            currNum = vec[i].first;
            currFreq = vec[i].second;
        }
    }

    const pair<int, int> p = {currNum, currFreq};

    // Frequency is >=2 so merge and push to rlslp
    // Merged variable is terminal for a new_slg, so it is marked negative
    if(p.second >= 2) {
        if(m.find(p) == m.end()) {
            m[p] = recompression_rlslp->nonterm.size();
            recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', abs(p.first), p.second));
        }
        
        //  ** Negative **
        cap_rhs.push_back(-m[p]);
    }
    else {
        cap_rhs.push_back(p.first);
    }

    new_slg_nonterm_vec.emplace_back(cap_rhs);

    return;
}

// Block Compression
unique_ptr<SLG> BComp(unique_ptr<SLG> & slg, unique_ptr<RecompressionRLSLP> & recompression_rlslp, unordered_map<pair<int, int>, int, hash_pair> & m) {

    // Current slg non-term list
    vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

    // // New SLG that needs to be created by applying BComp
    // unique_ptr<SLG> new_slg = make_unique<SLG>();

    // New SLG non-term list that new_slg needs.
    vector<SLGNonterm> new_slg_nonterm_vec;

    const int &grammar_size = slg_nonterm_vec.size();

    // We shall iterate throught each production rule in the increasing order of variable.
    // 'i' --> represents the variable.
    for(SLGNonterm & slg_nonterm : slg_nonterm_vec) {
        // Take the RHS of the production.
        vector<int> &rhs = slg_nonterm.rhs;

        if(rhs.empty()) {
            new_slg_nonterm_vec.emplace_back();
            continue;
        }

        // Create the expansion of RHS.
        vector<pair<int, int>> rhs_expansion;

        // Compute the Expansion.
        for(const int &rhs_symbol: rhs) {

            // Single Terminal
            if(rhs_symbol < 0) {
                rhs_expansion.emplace_back(rhs_symbol, 1);
                continue;
            }

            // **LR** of a current variable(rhs_symbol) in current RHS is **not empty**.
            if(slg_nonterm_vec[rhs_symbol].LR.second != -1) {
                rhs_expansion.push_back(slg_nonterm_vec[rhs_symbol].LR);
            }

            // Cap is not empty --> in new SLG the variable(rhs_symbol) RHS is not empty --> then Cap is not empty.
            if(new_slg_nonterm_vec[rhs_symbol].rhs.size() != 0) {
                rhs_expansion.emplace_back(rhs_symbol, 0);
            }

            // **RR** of a current variable(rhs_symbol) in current RHS is **not empty**.
            if(slg_nonterm_vec[rhs_symbol].RR.second != -1) {
                rhs_expansion.push_back(slg_nonterm_vec[rhs_symbol].RR);
            }
        }

        rhs.clear();
        vector<int>().swap(rhs);

        // Compute LR
        int lr_pointer = 1;

        pair<int, int> LR = rhs_expansion[0];

        while(lr_pointer < rhs_expansion.size()) {
            if(LR.first == rhs_expansion[lr_pointer].first) {
                LR.second += rhs_expansion[lr_pointer].second;
            }
            else {
                break;
            }

            lr_pointer++;
        }

        // set LR
        slg_nonterm.LR = LR;


        // Compute RR
        int rr_pointer = rhs_expansion.size()-2;

        pair<int, int> RR = rhs_expansion.back();

        // Case 1 : Everything is consumed by lr_pointer
        if(lr_pointer == rhs_expansion.size()) {
            // Cap is empty; set Cap
            new_slg_nonterm_vec.emplace_back();
            // RR is empty; set RR
            slg_nonterm.RR = {-1, -1};
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
            slg_nonterm.RR = RR;

            // Compress Cap(middle part)
            combineFrequenciesInRange(rhs_expansion, lr_pointer, rr_pointer, new_slg_nonterm_vec, m, recompression_rlslp);

            rhs_expansion.clear();
            vector<pair<int, int>>().swap(rhs_expansion);
        }
    }

    // Add new Starting Variable to the new Grammar G'(G Prime).

    const int &start_var = slg_nonterm_vec.size()-1;

    const pair<int, int> &start_var_LR = slg_nonterm_vec[start_var].LR;
    const pair<int, int> &start_var_RR = slg_nonterm_vec[start_var].RR;

    vector<int> new_start_rhs;

    // This always holds true.
    if(start_var_LR.second != -1) {
        if(start_var_LR.second >= 2) {
            if(m.find(start_var_LR) == m.end()) {
                m[start_var_LR] = recompression_rlslp->nonterm.size();
                recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', abs(start_var_LR.first), start_var_LR.second));
            }
            new_start_rhs.push_back(-m[start_var_LR]);
        }
        else {
            new_start_rhs.push_back(start_var_LR.first);
        }
    }

    if(!new_slg_nonterm_vec[start_var].rhs.empty()) {
        new_start_rhs.push_back(start_var);
    }

    if(start_var_RR.second != -1) {
        if(start_var_RR.second >= 2) {
            if(m.find(start_var_RR) == m.end()) {
                m[start_var_RR] = recompression_rlslp->nonterm.size();
                recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', abs(start_var_RR.first), start_var_RR.second));
            }
            new_start_rhs.push_back(-m[start_var_RR]);
        }
        else {
            new_start_rhs.push_back(start_var_RR.first);
        }
    }

    new_slg_nonterm_vec.push_back(SLGNonterm(new_start_rhs));

    slg.reset();

    return make_unique<SLG>(new_slg_nonterm_vec);
}

pair<int, int> computeAdjListHelper(int var, unique_ptr<SLG> & slg, vector<array<int, 4>> & adjList, vector<pair<int, int>> & dp) {

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

    const vector<int> &rhs = slg_nonterm.rhs;

    // 
    if(rhs.empty()) {
        return {0, 0};
    }
    else if(rhs.size() == 1 && rhs[0] < 0) {
        // Hopefully this should be always true when there is only 1 character to right
        // assert(rhs[0] < 0);

        return dp[var] = {rhs[0], rhs[0]};

        // return dp[var] = {rhs[0], rhs[0]};
    }
    
    vector<pair<int, int>> lms_rms_list;

    for(const int &rhs_symbol : rhs) {
        pair<int, int> lms_rms = computeAdjListHelper(rhs_symbol, slg, adjList, dp);
        lms_rms_list.push_back(lms_rms);
    }

    for(int i = 0; i < lms_rms_list.size() - 1; i++) {

        int f = lms_rms_list[i].second;
        int s = lms_rms_list[i+1].first;

        bool swapped = false;

        if(abs(f) < abs(s)) {
            swap(f, s);
            swapped = true;
        }
        adjList.push_back({f, s, swapped ? 1 : 0, slg_nonterm.vOcc});
    }

    slg_nonterm.LMS = lms_rms_list.front().first;
    slg_nonterm.RMS = lms_rms_list.back().second;

    return dp[var] = {slg_nonterm.LMS, slg_nonterm.RMS};
}

void computeAdjList(unique_ptr<SLG> & slg, vector<array<int, 4>> &adjList) {
    vector<pair<int, int>> dp(slg->nonterm.size(), make_pair(1, 1));

    computeAdjListHelper(slg->nonterm.size()-1, slg, adjList, dp);

    return;
}

vector<array<int, 4>> computeAdjList(unique_ptr<SLG> & slg) {

    vector<array<int, 4>> adjList;

    vector<pair<int, int>> dp(slg->nonterm.size(), make_pair(1, 1));

    computeAdjListHelper(slg->nonterm.size()-1, slg, adjList, dp);

    return adjList;
}

int computeVOccHelper(vector<vector<pair<int,int>>> & graph, int u, vector<int> & dp) {

    // Base Case : Target is Reached / Target is Same as the current node.
    if(u == graph.size()-1) {
        return dp[u] = 1;
    }

    // Already Computed
    if(dp[u] != -1) {
        return dp[u];
    }

    int num_paths = 0;

    for(const auto & edge : graph[u]) {
        int v = edge.first;
        int weight = edge.second;

        num_paths += weight * computeVOccHelper(graph, v, dp);
    }

    // u is explored.
    graph[u].clear();
    vector<pair<int, int>>().swap(graph[u]);

    return dp[u] = num_paths;
}
void computeVOcc(unique_ptr<SLG> & slg) {

    // Here, any list size is slg_nonterm_vec.size()

    vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

    // Create Reverse Graph.
    vector<vector<pair<int,int>>> graph(slg_nonterm_vec.size(), vector<pair<int, int>>());


    // Enumerate Each Production Rule.
    for(int i=0; i<slg_nonterm_vec.size(); i++) {

        // Current SLGNonterm
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        unordered_map<int, int> var_freq;

        const vector<int> & rhs = slg_nonterm.rhs;

        // Enumerate each character of RHS
        // Frequency calculation for reverse of the graph.
        for(const int & var : rhs) {

            // Only Non-Terminals
            if(var >= 0) {
                var_freq[var]++;
            }
        }

        // Construct Edges.
        for(const auto & x : var_freq) {
            // edge from v to u, u <-- v
            int u = i;
            int v = x.first;
            int weight = x.second;

            // Weighted Edge from v to u , v --> u
            graph[v].emplace_back(u, weight);
        }
    }

    // Store VOcc.
    vector<int> dp(slg_nonterm_vec.size(), -1);

    // Compute vOcc.
    for(int i = 0; i < slg_nonterm_vec.size(); i++) {
        computeVOccHelper(graph, i, dp);
    }

    // Update SLGNonterm.
    for(int i=0; i<dp.size(); i++) {
        slg_nonterm_vec[i].vOcc = dp[i];
    }

    return;
}

void sortAdjList(vector<array<int, 4>> & adjList) {

    // Make Positive
    for(array<int, 4> & a : adjList) {
        a[0] = -a[0];
        a[1] = -a[1];
    }

    // Sort It.
    radixSort(adjList);


    // Revert to Negative
    for(array<int, 4> & a : adjList) {
        a[0] = -a[0];
        a[1] = -a[1];
    }
}

void createPartition(const vector<array<int, 4>> & adjList, array<unordered_set<int>, 2> &partition_set) {

    // // Make Positive
    // for(array<int, 4> & a : adjList) {
    //  a[0] = -a[0];
    //  a[1] = -a[1];
    // }
    
    unordered_set<int>& leftSet = partition_set[0];
    unordered_set<int>& rightSet = partition_set[1];
    int currentIndex = 0;
    size_t n = adjList.size();

    int c = adjList[currentIndex][0];

    while(currentIndex < n) {
        int leftSetFreq = 0;
        int rightSetFreq = 0;
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
    
    int LRPairsCount = 0;
    int RLPairsCount = 0;

    /*
        for (int i = 0; i < arr.size() - 1; i++) {


            LRPairsCount += (leftSet.find(arr[i]) != leftSet.end()) && (rightSet.find(arr[i + 1]) != rightSet.end());
            RLPairsCount += (rightSet.find(arr[i]) != rightSet.end()) && (leftSet.find(arr[i + 1]) != leftSet.end());
        }
    */

    for(const array<int, 4> & arr : adjList) {
        int f = arr[0];
        int s = arr[1];

        if(arr[2] == 1) {
            swap(f, s);
        }

        LRPairsCount += ((leftSet.find(f) != leftSet.end()) && (rightSet.find(s) != rightSet.end())) ? arr[3] : 0;
        RLPairsCount += ((rightSet.find(f) != rightSet.end()) && (leftSet.find(s) != leftSet.end())) ? arr[3] : 0;

    }

    if (RLPairsCount < LRPairsCount) {
        swap(leftSet, rightSet);
    }

    swap(rightSet, leftSet);

 //    // Revert to Negative
    // for(array<int, 4> & a : adjList) {
    //  a[0] = -a[0];
    //  a[1] = -a[1];
    // }

    

    return;
}

array<unordered_set<int>, 2> createPartition(const vector<array<int, 4>> & adjList) {

    // // Make Positive
    // for(array<int, 4> & a : adjList) {
    //  a[0] = -a[0];
    //  a[1] = -a[1];
    // }
    
    unordered_set<int> leftSet, rightSet;
    int currentIndex = 0;
    size_t n = adjList.size();

    int c = adjList[currentIndex][0];

    while(currentIndex < n) {
        int leftSetFreq = 0;
        int rightSetFreq = 0;
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
    
    int LRPairsCount = 0;
    int RLPairsCount = 0;

    /*
        for (int i = 0; i < arr.size() - 1; i++) {


            LRPairsCount += (leftSet.find(arr[i]) != leftSet.end()) && (rightSet.find(arr[i + 1]) != rightSet.end());
            RLPairsCount += (rightSet.find(arr[i]) != rightSet.end()) && (leftSet.find(arr[i + 1]) != leftSet.end());
        }
    */

    for(const array<int, 4> & arr : adjList) {
        int f = arr[0];
        int s = arr[1];

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

// Pair-Wise Compression
unique_ptr<SLG> PComp(unique_ptr<SLG> & slg, unique_ptr<RecompressionRLSLP> & recompression_rlslp,  unordered_map<pair<int, int>, int, hash_pair> & m) {
   // Compute vOcc
    computeVOcc(slg);

    // Compute AdjList
    vector<array<int, 4>> adjList;
    computeAdjList(slg, adjList);

    sortAdjList(adjList);

    // Create Partition.
    array<unordered_set<int>,2> arr;
    createPartition(adjList, arr);

    adjList.clear();  // Clear immediately
    vector<array<int, 4>>().swap(adjList);  // Ensure memory is released

    const unordered_set<int> &left_set = arr[0], &right_set = arr[1];

    // Current slg non-term list
    vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

    // New SLG that needs to be created by applying BComp
    // unique_ptr<SLG> new_slg = make_unique<SLG>();

    // New SLG non-term list that new_slg needs.
    vector<SLGNonterm> new_slg_nonterm_vec;


    // We shall iterate throught each production rule in the increasing order of variable.
    // 'i' --> represents the variable.
    for(int i=0; i<slg_nonterm_vec.size(); i++) {
        SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

        const vector<int> & rhs = slg_nonterm.rhs;

        if(rhs.empty()) {
            new_slg_nonterm_vec.push_back(SLGNonterm());
            continue;
        }
        if((int)rhs.size() >= 2) {

            vector<int> rhs_expansion;

            // Expanding RHS.
            for(int j=0; j<rhs.size(); j++) {


                // RHS Terminal
                if(rhs[j] < 0) {

                    rhs_expansion.push_back(rhs[j]);

                    if(j==0) {
                        if(left_set.find(rhs[j]) != left_set.end()) {
                            slg_nonterm.LB = 0;
                        }
                        else if(right_set.find(rhs[j]) != right_set.end()) {
                            slg_nonterm.LB = rhs[j];

                            // Remove the pushed element
                            rhs_expansion.pop_back();
                        }
                    }
                    else if(j == rhs.size() - 1) {
                        if(left_set.find(rhs[j]) != left_set.end()) {
                            slg_nonterm.RB = rhs[j];

                            // Remove the pushed element
                            rhs_expansion.pop_back();
                        }
                        else if(right_set.find(rhs[j]) != right_set.end()) {
                            slg_nonterm.RB = 0;
                        }
                    }

                }

                // RHS Non-Terminal
                else {

                    if(j==0) {
                        slg_nonterm.LB = slg_nonterm_vec[rhs[j]].LB;

                    }
                    else if(j == rhs.size() - 1) {
                        slg_nonterm.RB = slg_nonterm_vec[rhs[j]].RB;
                    }
                    // LB --> not equal to 0
                    // If j == 0 we set it to LB in slg_nonterm but not in to the expansion.
                    if((j != 0) and slg_nonterm_vec[rhs[j]].LB) {
                        rhs_expansion.push_back(slg_nonterm_vec[rhs[j]].LB);
                    }

                    // Check whether Cap is Empty in new_slg_nonterm_vec.
                    if(new_slg_nonterm_vec[rhs[j]].rhs.empty() == false){
                        rhs_expansion.push_back(rhs[j]);
                    }



                    // RB --> not equal to 0.
                    // If j is last element, we set it to RB in slg_nonterm but not in to the expansion
                    if((j != rhs.size() - 1) and slg_nonterm_vec[rhs[j]].RB) {
                        rhs_expansion.push_back(slg_nonterm_vec[rhs[j]].RB);
                    }
                }
            }

            // To handle last character/corner case
            // if(rhs_expansion.size()>=2)
            rhs_expansion.push_back(rhs_expansion.back());



            vector<int> cap_rhs;

            for(int j=0; j<(int)rhs_expansion.size()-1; j++) {
                if(rhs_expansion[j] < 0 && rhs_expansion[j+1] < 0) {
                    if(left_set.find(rhs_expansion[j]) != left_set.end() && right_set.find(rhs_expansion[j+1]) != right_set.end()) {
                        if(m.find({rhs_expansion[j], rhs_expansion[j+1]}) == m.end()) {
                            m[{rhs_expansion[j], rhs_expansion[j+1]}] = recompression_rlslp->nonterm.size();
                            recompression_rlslp->nonterm.push_back(RLSLPNonterm('1', abs(rhs_expansion[j]), abs(rhs_expansion[j+1])));
                        }

                        cap_rhs.push_back(-m[{rhs_expansion[j], rhs_expansion[j+1]}]);
                        j++;
                    }
                    else {
                        cap_rhs.push_back(rhs_expansion[j]);
                    }
                }
                else if(rhs_expansion[j] >= 0) {
                    cap_rhs.push_back(rhs_expansion[j]);
                }
                else {
                    cap_rhs.push_back(rhs_expansion[j]);
                }
            }


            new_slg_nonterm_vec.push_back(SLGNonterm(cap_rhs));
        }
        else if((int)rhs.size() == 1) {

            // Terminal;
            // A --> a
            // B --> b
            if(rhs[0] < 0) {
                if(left_set.find(rhs[0]) != left_set.end()) {
                    slg_nonterm.LB = 0;
                    new_slg_nonterm_vec.push_back(SLGNonterm());
                    slg_nonterm.RB = rhs[0];
                }
                else if(right_set.find(rhs[0]) != right_set.end()) {
                    slg_nonterm.LB = rhs[0];
                    new_slg_nonterm_vec.push_back(SLGNonterm());
                    slg_nonterm.RB = 0;
                }
                else {
                    // Non-reachable Non-Terminals.
                    // A --> a
                    // 'a' is not found in adjacency list so push empty.
                    new_slg_nonterm_vec.push_back(SLGNonterm());
                    cout << "Error: Not Found in Left and Right." << endl;
                }
            }
            else {

                // A --> B (LB(B) . B^ . RB(B))
                // LB(A) = LB(B)
                // RB(A) = RB(B)
                // A^ = B^ (only if B^ is not empty)
                vector<int> cap_rhs;

                slg_nonterm.LB = slg_nonterm_vec[rhs[0]].LB;

                if(new_slg_nonterm_vec[rhs[0]].rhs.empty() == false) {
                    cap_rhs.push_back(rhs[0]);
                }

                slg_nonterm.RB = slg_nonterm_vec[rhs[0]].RB;

                new_slg_nonterm_vec.push_back(SLGNonterm(cap_rhs));
            }
            
        }
    }

    arr[0].clear();
    arr[1].clear();
    unordered_set<int>().swap(arr[0]);
    unordered_set<int>().swap(arr[1]);



    int start_var = slg_nonterm_vec.size()-1;

    int start_var_LB = slg_nonterm_vec[start_var].LB;
    int start_var_RB = slg_nonterm_vec[start_var].RB;


    vector<int> new_start_rhs;

    if(start_var_LB != 0) {
        new_start_rhs.push_back(start_var_LB);
    }

    if(new_slg_nonterm_vec[start_var].rhs.empty() == false) {
        new_start_rhs.push_back(start_var);
    }
    
    if(start_var_RB != 0) {
        new_start_rhs.push_back(start_var_RB);
    }

    new_slg_nonterm_vec.push_back(SLGNonterm(new_start_rhs));

    slg.reset();

    return make_unique<SLG>(new_slg_nonterm_vec);
}

unique_ptr<RecompressionRLSLP> recompression_on_slp(unique_ptr<InputSLP>& s) {

    unique_ptr<SLG> slg = make_unique<SLG>();

    unique_ptr<RecompressionRLSLP> recompression_rlslp = make_unique<RecompressionRLSLP>();

    unordered_map<pair<int, int>, int, hash_pair> m;

    // For 0.
    recompression_rlslp->nonterm.push_back(RLSLPNonterm());

    // cout << "OK" << endl;


    // Compute S0 from S and Initialize Recompression
    for(int i=0; i<s->nonterm.size(); i++) {
        char type = s->nonterm[i].type;
        int first = s->nonterm[i].first;
        int second = s->nonterm[i].second;



        vector<int> rhs;

        if(type == '0') {
            // For type '0', in SLG, terminals start from -1.
            // Initially, -1 = recompression_rlslp->nonterm size.
            rhs = {-(int)(recompression_rlslp->nonterm).size()};
            // Here first is the ASCII values of the character.
            recompression_rlslp->nonterm.push_back(RLSLPNonterm('0', first, second));
        }
        else {
            rhs = {first, second};
        }

        slg->nonterm.push_back(SLGNonterm(rhs));

    }

    s.reset();

    // vector<int> arr1 = expandSLG(slg);

    // cout << arr << endl;


    int i = 0;

    // printSLG(slg);

    while(++i) {

        const vector<SLGNonterm> &slg_nonterm_vec = slg->nonterm;

        if(slg_nonterm_vec.back().rhs.size() == 1 && slg_nonterm_vec.back().rhs.back() < 0) {
            break;
        }

        if(i&1) {
            slg = BComp(slg, recompression_rlslp, m);
            // cout << i << ' ' << "BComp" << endl;
        }
        else {
            slg = PComp(slg, recompression_rlslp, m);
            // cout << i << ' ' << "PComp" << endl;
        }

        // printSLG(slg);
        // printRecompressionRLSLP(recompression_rlslp);
    }

    cout << "Runs: " << i << endl;

    // vector<int> arr2 = expandRLSLP(recompression_rlslp);

    // if(arr1 != arr2) {
    //     cout << "SLG and RLSLP expansion didn't match!" << endl;
    //     exit(1);
    // }


    

    return recompression_rlslp;
}

int computeExplen(const int i, vector<RLSLPNonterm>& rlslp_nonterm_vec) {
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

unique_ptr<InputSLP> getSLP(int grammar_size) {

    vector<SLPNonterm> nonterm;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    srand(seed);



    // Generate random number X < n/2
    int num_terminals = rand() % ((int)((3*grammar_size)/4));

    num_terminals++;

    for(int i=1; i<=num_terminals; i++) {
        nonterm.push_back(SLPNonterm('0', -i, 1));
    }


    for (int i = num_terminals; i < grammar_size; ++i) {
        int num1 = rand() % i; // Generate first number less than i
        int num2 = rand() % i; // Generate second number less than i


        int bit1 = rand()%2;
        int bit2 = rand()%2;

        int num = (bit2 == 1 ? num2 : num1);

        if(bit1 == 0)
        nonterm.push_back(SLPNonterm('1', i-1, num));
        else {
            nonterm.push_back(SLPNonterm('1', num, i-1));
        }
    }

    unique_ptr<InputSLP> slp = make_unique<InputSLP>(nonterm);

    return slp;
}

vector<pair<int, int>> get_random_queries(int text_size) {
    // Create a vector to store the pairs
    vector<pair<int, int>> pairs;
    pairs.reserve(1000000); // Reserve space for 4096 pairs to avoid reallocation

    // Random number generation setup
    random_device rd;  // Obtain a random number from hardware
    mt19937 gen(rd()); // Seed the generator
    uniform_int_distribution<> distrib(0, text_size - 1);

    // Generate 4096 pairs
    for (int i = 0; i < 1000000; ++i) {
        int first = distrib(gen);
        int second = distrib(gen);
        pairs.emplace_back(first, second); // Emplace pair into the vector
    }

    return pairs;
}


void test(int text_size, unique_ptr<RecompressionRLSLP> &recompression_rlslp, vector<RLSLPNonterm> & rlslp_nonterm_vec) {

    vector<int> arr = expandRLSLP(recompression_rlslp);

    ifstream file("einstein.en.txt");

    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        exit(1);
    }

    vector<int> text;

    char ch;
    while (file.get(ch)) {
        text.push_back(static_cast<unsigned char>(ch));
    }

    // Close the file
    file.close();

    assert(text.size()==arr.size());

    cout << "Text Size Matched!" << endl;

    for(int i=0; i<text.size(); i++) {
        if(text[i] != arr[i]) {
            cout << "Text didn't match : " << i << ' ' << text[i] << ' ' << arr[i] << endl;
            exit(1);
        }
    }

    cout << "Text Matched!" << endl;

    vector<pair<int, int>> random_queries = get_random_queries(text_size);

    auto start_time = std::chrono::high_resolution_clock::now();

    for(auto x : random_queries) {

        int i = x.first;
        int j = x.second;
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

        int res1 = LCE(v1, v2, i, v1_ancestors, v2_ancestors, rlslp_nonterm_vec);

        int res2 = 0;

        int ii = i;
        int jj = j;

        while(jj < text_size && arr[ii] == arr[jj]) {
            res2++;

            jj++;
            ii++;
        }

        // cout << i << ' ' << j << ' ' << res1 << ' ' << res2 << endl;
        // assert(res1==res2);

        if(res1 != res2) {
            cout << "ERROR" << endl;
            exit(1);
        }


    }

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Output the duration
    cout << "Time taken for LCE Queries: " << duration_seconds << " seconds" << endl;
}

void start_compression(string input_file) {
    unique_ptr<InputSLP> inputSLP = make_unique<InputSLP>();
    inputSLP->read_from_file(input_file);

    int j = 0;
    cout << "Number of Non-Terminals : " << inputSLP->nonterm.size() << endl;

    int i = -1;

    auto start_time = std::chrono::high_resolution_clock::now();

    unique_ptr<RecompressionRLSLP> recompression_rlslp = recompression_on_slp(inputSLP);

    vector<RLSLPNonterm> & rlslp_nonterm_vec = recompression_rlslp->nonterm;

    for(int i=rlslp_nonterm_vec.size()-1; i>=1; i--) {
        computeExplen(i, rlslp_nonterm_vec);
    }

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Output the duration
    cout << "Time taken for Construction: " << duration_seconds << " seconds" << endl;

    // // printRecompressionRLSLP(recompression_rlslp);
    int text_size = rlslp_nonterm_vec.back().explen;

    cout << "Text Size : " << text_size << endl;

    // TEST
    // test(text_size, recompression_rlslp, rlslp_nonterm_vec);
}

int main(int argc, char *argv[]) {

    if(argc < 2) {
        cout << "Please provided the Input SLP file!" << endl;
        exit(1);
    }

    string input_file(argv[1]);

    auto start_time = std::chrono::high_resolution_clock::now();

    start_compression(input_file);

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    // Output the duration
    cout << "Total Time taken: " << duration_seconds << " seconds" << endl;

    return 0;
}
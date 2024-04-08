#include<bits/stdc++.h>
#include "RecompressionDefinitions.hpp"
#include "radix_sort.h"
#include "io.h"
using namespace std;

void printRecompressionRLSLP(const unique_ptr<RecompressionRLSLP> & recompression_rlslp) {

	cout << "RECOMPRESSION PRINTING STARTED..." << endl;

	if(!(recompression_rlslp->nonterm).size()) {
		cout << "RecompressionRLSLP is Empty!!" << endl;
		return;
	}

	int i = 0;

	for(const RLSLPNonterm & rlslp_nonterm : recompression_rlslp->nonterm) {
		cout << i << " --> ";
		cout << rlslp_nonterm.type << ' ' << rlslp_nonterm.first << ' ' << rlslp_nonterm.second << endl;
		i++;
	}

	cout << "RECOMPRESSION PRINTING ENDED!" << endl << endl;

	return;
}

void printSLG(const unique_ptr<SLG> & slg) {
	cout << "SLG PRINTING STARTED..." << endl;

	if(!(slg->nonterm).size()) {
		cout << "SLG is Empty!!" << endl;
		return;
	}

	int i = 0;

	for(const SLGNonterm & slg_nonterm : slg->nonterm) {
		cout << i << " --> ";
		for(int x : slg_nonterm.rhs) {
			cout << x << ' ';
		}
		cout << endl;
		i++;
	}

	cout << "SLG PRINTING ENDED!" << endl << endl;

	return;
}


vector<pair<int, int>> combineFrequenciesInRange(vector<pair<int, int>>& vec, int lr_pointer, int rr_pointer) {
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
            result.push_back({currNum, currFreq});
            // Move to the next number
            currNum = vec[i].first;
            currFreq = vec[i].second;
        }
    }

    // Add the last pair to result vector
    result.push_back(make_pair(currNum, currFreq));

    return result;
}

// Block Compression
unique_ptr<SLG> BComp(unique_ptr<SLG> & slg, unique_ptr<RecompressionRLSLP> & recompression_rlslp, map<pair<int, int>, int> & m) {

	// Current slg non-term list
	vector<SLGNonterm> & slg_nonterm_vec = slg->nonterm;

	// // New SLG that needs to be created by applying BComp
	// unique_ptr<SLG> new_slg = make_unique<SLG>();

	// New SLG non-term list that new_slg needs.
	vector<SLGNonterm> new_slg_nonterm_vec;


	// We shall iterate throught each production rule in the increasing order of variable.
	// 'i' --> represents the variable.
	for(int i=0; i<slg_nonterm_vec.size(); i++) {

		// Take the current production rule.
		SLGNonterm & slg_nonterm = slg_nonterm_vec[i];

		// Take the RHS of the production.
		vector<int> rhs = slg_nonterm.rhs;

		// Create the expansion of RHS.
		vector<pair<int, int>> rhs_expansion;

		// Compute the Expansion.
		for(int j=0; j<rhs.size(); j++) {

			// Single Terminal
			if(rhs[j] < 0) {
				rhs_expansion.push_back({rhs[j], 1});
				continue;
			}

			// **LR** of a current variable(rhs[j]) in current RHS is **not empty**.
			if(slg_nonterm_vec[rhs[j]].LR.second != -1) {
				rhs_expansion.push_back(slg_nonterm_vec[rhs[j]].LR);
			}

			// Cap is not empty --> in new SLG the variable(rhs[j]) RHS is not empty --> then Cap is not empty.
			if(new_slg_nonterm_vec[rhs[j]].rhs.size() != 0) {
				rhs_expansion.push_back({rhs[j], 0});
			}

			// **RR** of a current variable(rhs[j]) in current RHS is **not empty**.
			if(slg_nonterm_vec[rhs[j]].RR.second != -1) {
				rhs_expansion.push_back(slg_nonterm_vec[rhs[j]].RR);
			}
		}

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
			new_slg_nonterm_vec.push_back(SLGNonterm());
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
			vector<pair<int, int>> cap_compressed_slg_nonterm_vec = combineFrequenciesInRange(rhs_expansion, lr_pointer, rr_pointer);

			// There are no pairs presents for the cap.
			if(cap_compressed_slg_nonterm_vec.size() == 0) {
				// Cap is empty; set Cap
				new_slg_nonterm_vec.push_back(SLGNonterm());
			}
			else {

				// RHS of new variable CAP.
				vector<int> cap_rhs;

				for(const pair<int, int> & p : cap_compressed_slg_nonterm_vec) {

					// Frequency is >=2 so merge and push to rlslp
					// Merged variable is terminal for a new_slg, so it is marked negative
					if(p.second >= 2) {
						if(m.find(p) == m.end()) {
							m[p] = recompression_rlslp->nonterm.size();
							recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', p.first, p.second));
						}
						
						//  ** Negative **
						cap_rhs.push_back(-m[p]);
					}
					else {
						cap_rhs.push_back(p.first);
					}
					
				}

				new_slg_nonterm_vec.push_back(SLGNonterm(cap_rhs));
			}
		}
	}

	// Add new Starting Variable to the new Grammar G'(G Prime).

	int start_var = slg_nonterm_vec.size()-1;

	pair<int, int> start_var_LR = slg_nonterm_vec[start_var].LR;
	pair<int, int> start_var_RR = slg_nonterm_vec[start_var].RR;

	vector<int> new_start_rhs;

	// This always holds true.
	if(start_var_LR.second != -1) {
		if(start_var_LR.second >= 2) {
			if(m.find(start_var_LR) == m.end()) {
				m[start_var_LR] = recompression_rlslp->nonterm.size();
				recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', start_var_LR.first, start_var_LR.second));
			}

			new_start_rhs.push_back(m[start_var_LR]);
		}
		else {
			new_start_rhs.push_back(start_var_LR.first);
		}
	}

	new_start_rhs.push_back(start_var);

	if(start_var_RR.second != -1) {
		if(start_var_RR.second >= 2) {
			if(m.find(start_var_RR) == m.end()) {
				m[start_var_RR] = recompression_rlslp->nonterm.size();
				recompression_rlslp->nonterm.push_back(RLSLPNonterm('2', start_var_RR.first, start_var_LR.second));
			}

			new_start_rhs.push_back(m[start_var_RR]);
		}
		else {
			new_start_rhs.push_back(start_var_RR.first);
		}
	}

	new_slg_nonterm_vec.push_back(SLGNonterm(new_start_rhs));


	return make_unique<SLG>(new_slg_nonterm_vec);
}


// Return {RMS, LMS} DP? ** doubt **
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
	SLGNonterm & slg_nonterm = slg->nonterm[var];

	vector<int> & rhs = slg_nonterm.rhs;

	// 
	if(rhs.size() == 0) {
		return {0, 0};
	}
	else if(rhs.size() == 1) {
		// Hopefully this should be always true when there is only 1 character to right
		assert(rhs[0] < 0);

		return dp[var] = {rhs[0], rhs[0]};
	}
	
	vector<pair<int, int>> lms_rms_list;

	for(int i=0; i<rhs.size(); i++) {
		pair<int, int> lms_rms = computeAdjListHelper(rhs[i], slg, adjList, dp);
		lms_rms_list.push_back(lms_rms);
	}

	for(int i=0; i<lms_rms_list.size()-1; i++) {

		int f = lms_rms_list[i].second;
		int s = lms_rms_list[i+1].first;

		int swapped = 0;

		if(abs(f) < abs(s)) {
			swap(f, s);
			swapped = 1;
		}
		// cout << var << ' ' << f << ' ' << s << ' ' << swapped << ' ' << slg_nonterm.vOcc << endl;
		adjList.push_back({f, s, swapped, slg_nonterm.vOcc});
	}

	slg_nonterm.LMS = lms_rms_list[0].first;
	slg_nonterm.RMS = lms_rms_list.back().second;

	return dp[var] = {slg_nonterm.LMS, slg_nonterm.RMS};
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

	for(auto & x : graph[u]) {
		int v = x.first;
		int weight = x.second;

		num_paths += weight * computeVOccHelper(graph, v, dp);
	}

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

		vector<int> & rhs = slg_nonterm.rhs;

		// Enumerate each character of RHS
		// Frequency calculation for reverse of the graph.
		for(int & var : rhs) {
			var_freq[var]++;
		}

		// Construct Edges.
		for(const auto & x : var_freq) {
			// edge from v to u, u <-- v
			int u = i;
			int v = x.first;
			int weight = x.second;

			// Weighted Edge from v to u , v --> u
			graph[v].push_back({u, weight});
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

// only adjList as input. *
array<unordered_set<int>, 2> createPartition(vector<array<int, 4>> & adjList) {
    
    unordered_set<int> leftSet, rightSet;
    int currentIndex = 0;
    size_t n = adjList.size();

    int c = adjList[currentIndex][0];

    for(int i=1;i<=c-1;i++) {
    	leftSet.insert(i);
    }

    while(currentIndex < n) {
        int leftSetFreq = 0;
        int rightSetFreq = 0;
        while (currentIndex < n && adjList[currentIndex][0] == c) {
            if (leftSet.find(adjList[currentIndex][1]) != leftSet.end()) {
                leftSetFreq++;
            } else {
                rightSetFreq++;
            }
            currentIndex++;
        }


        if (leftSetFreq >= rightSetFreq) {
            rightSet.insert(c);
        } else {
            leftSet.insert(c);
        }

        if(currentIndex < n) {
        	
        	for(int i=c+1; i<=adjList[currentIndex][0]-1;i++) {
        		leftSet.insert(i);
        	}
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

	for(array<int, 4> & arr : adjList) {
		int f = arr[0];
		int s = arr[1];

		if(arr[2] == 1) {
			swap(f, s);
		}

		LRPairsCount += (leftSet.find(f) != leftSet.end()) && (rightSet.find(s) != rightSet.end());
	    RLPairsCount += (rightSet.find(f) != rightSet.end()) && (leftSet.find(s) != leftSet.end());

	}

    if (RLPairsCount < LRPairsCount) {
        swap(leftSet, rightSet);
    }

    return { rightSet, leftSet };
}

// Pair-Wise Compression
unique_ptr<SLG> PComp(unique_ptr<SLG> & slg, unique_ptr<RecompressionRLSLP> & recompression_rlslp,  map<pair<int, int>, int> & m) {
	
	// Compute vOcc - Sample Passed
	computeVOcc(slg);

	


	// Compute AdjList
	// ** doubt ** --> I hope just call to the largest varaible/number would suffice.
	vector<array<int, 4>> adjList = computeAdjList(slg);

	sortAdjList(adjList);

	cout << adjList << endl;


	



	

	// Create Partition.

	unordered_set<int> left_set, right_set;

	array<unordered_set<int>,2> arr = createPartition(adjList);

	left_set = arr[0];
	right_set = arr[1];

	cout << "Left set : " << left_set << endl;
	cout << "Right set : " << right_set << endl;


	return make_unique<SLG>();

	



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

		vector<int> & rhs = slg_nonterm.rhs;

		// Probably only terminal!
		if(rhs.size() >= 2) {

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
					if(new_slg_nonterm_vec[rhs[j]].rhs.size() != 0){
						rhs_expansion.push_back(rhs[j]);
					}

					// RB --> not equal to 0.
					// If j is last element, we set it to RB in slg_nonterm but not in to the expansion
					if((j != rhs.size() - 1) and slg_nonterm_vec[rhs[j]].RB) {
						rhs_expansion.push_back(slg_nonterm_vec[rhs[j]].RB);
					}
				}

				// To handle last character/corner case
				rhs_expansion.push_back(rhs_expansion.back());

				vector<int> cap_rhs;

				for(int j=0; j<rhs_expansion.size()-1; j++) {
					if(rhs_expansion[j] < 0 && rhs_expansion[j+1] < 0) {
						if(left_set.find(rhs_expansion[j]) != left_set.end() && right_set.find(rhs_expansion[j+1]) != right_set.end()) {
							if(m.find({rhs_expansion[j], rhs_expansion[j+1]}) == m.end()) {
								m[{rhs_expansion[j], rhs_expansion[j+1]}] = recompression_rlslp->nonterm.size();
								recompression_rlslp->nonterm.push_back(RLSLPNonterm('1', rhs_expansion[j], rhs_expansion[j+1]));
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
		}
		else if(rhs.size() == 1) {

			// Terminal;

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
				cout << "Error: Not Found in Left and Right." << endl;
			}
		}
		// Size 0. ** doubt ** are these required? 
		else {
			new_slg_nonterm_vec.push_back(SLGNonterm());
		}


	}

	return make_unique<SLG>(new_slg_nonterm_vec);





}
unique_ptr<RecompressionRLSLP> recompression_on_slp(unique_ptr<InputSLP>& s) {

	unique_ptr<SLG> slg = make_unique<SLG>();

	unique_ptr<RecompressionRLSLP> recompression_rlslp = make_unique<RecompressionRLSLP>();

	map<pair<int, int>, int> m;

	// Initiate RLSLP
	for(int i=0; i<s->nonterm.size(); i++) {
		char type = s->nonterm[i].type;
		int first = s->nonterm[i].first;
		int second = s->nonterm[i].second;

		recompression_rlslp->nonterm.push_back(RLSLPNonterm(type, first, second));
	}

	// Compute S0 from S
	for(int i=0; i<s->nonterm.size(); i++) {
		char type = s->nonterm[i].type;
		int first = s->nonterm[i].first;
		int second = s->nonterm[i].second;

		vector<int> rhs;

		if(type == '0') {
			rhs = {first};
		}
		else {
			rhs = {first, second};
		}

		slg->nonterm.push_back(SLGNonterm(rhs));
	}

	printRecompressionRLSLP(recompression_rlslp);
	printSLG(slg);

	unique_ptr<SLG> new_slg = BComp(slg, recompression_rlslp, m);

	printRecompressionRLSLP(recompression_rlslp);
	printSLG(new_slg);

	unique_ptr<SLG> new_slg_pcomp = PComp(new_slg, recompression_rlslp, m);





	return recompression_rlslp;
}

 


int main() {

	vector<SLPNonterm> nonterm;

	nonterm.push_back(SLPNonterm('0', -1, 1));
	nonterm.push_back(SLPNonterm('0', -2, 1));
	nonterm.push_back(SLPNonterm('0', -3, 1));
	nonterm.push_back(SLPNonterm('1', 2, 1));
	nonterm.push_back(SLPNonterm('1', 3, 0));
	nonterm.push_back(SLPNonterm('1', 4, 2));
	nonterm.push_back(SLPNonterm('1', 5, 4));
	nonterm.push_back(SLPNonterm('1', 6, 5));


	unique_ptr<InputSLP> inputSLP = make_unique<InputSLP>(nonterm);

	recompression_on_slp(inputSLP);

}
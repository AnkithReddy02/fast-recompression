#include<bits/stdc++.h>
using namespace std;

class SLP {
public:
	struct NonTerminal {
		char type;
		int first; 
		int second;
		int lml;
		int rml;
		int vocc;
	};

	vector<NonTerminal> grammar;


	int computeVOccHelper(vector<NonTerminal> & reverseGrammar, vector<int> & dp, int & nonterm_node) {

		if(nonterm_node == reverseGrammar.size()-1) {
			return dp[nonterm_node] = 1;
		}

		if(dp[nonterm_node] != -1) {
			return dp[nonterm_node];
		}

		if(reverseGrammar[nonterm_node].first == reverseGrammar[nonterm_node].second) {
			return dp[nonterm_node] = 2 * computeVOccHelper(reverseGrammar, dp, nonterm_node.first);
		}
		else {
			return dp[nonterm_node] = 1 * computeVOccHelper(reverseGrammar, dp, nonterm_node.first) + 1 * computeVOccHelper(reverseGrammar, dp, nonterm_node.second);
		}


		return 0;
	}

	pair<int, int> compute_lml_rml(vector<NonTerminal> & grammar, int node) {

		// Terminal
		if(grammar[node].type == '0') {
			grammar[node].lml = grammar[node].rml = grammar[node].first;
			return {grammar[node].first, grammar[node].first};
		}

		pair<int, int> left = compute_lml_rml(grammar, grammar[node].first);
		pair<int, int> right = compute_lml_rml(grammar, grammar[node].second);

		grammar[node].lml = left.first;
		grammar[node].rml = right.second;

		return {grammar[node].lml, grammar[node].rml};

	}

	void computeVOcc() {
		vector<NonTerminal> reverseGrammar;

		vector<int> dp(reverseGrammar.size(), -1);

		for(int i=0; i<reverseGrammar.size(); i++) {
			grammar[i].vocc = computeVOccHelper(reverseGrammar, dp, i);
		}
	}

	void createAdjacencyList() {

	}





};
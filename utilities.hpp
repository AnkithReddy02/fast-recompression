#ifndef UTILITIES_H
#define UTILITIES_H

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

void expandRLSLP(int var, const vector<RLSLPNonterm> & rlslp_nonterm_vec, vector<int> & result) {
	if(rlslp_nonterm_vec[var].type == '0') {
		result.push_back(var);
		return;
	}

	char type = rlslp_nonterm_vec[var].type;
	int first = rlslp_nonterm_vec[var].first;
	int second = rlslp_nonterm_vec[var].second;

	if(type == '1') {
		expandRLSLP(first, rlslp_nonterm_vec, result);
		expandRLSLP(second, rlslp_nonterm_vec, result);
	}
	// Type == '2'
	else {
		while(second--) {
			expandRLSLP(first, rlslp_nonterm_vec, result);
		}
	}
	return;
}
void expandRLSLP(const unique_ptr<RecompressionRLSLP> & recompression_rlslp) {

	vector<int> result;
	vector<RLSLPNonterm> rlslp_nonterm_vec = recompression_rlslp->nonterm;
	expandRLSLP(rlslp_nonterm_vec.size()-1, rlslp_nonterm_vec, result);

	cout << result << endl;
}

void expandSLG(int var, vector<SLGNonterm> & slg_nonterm_vec, vector<int> & result) {

	if(var < 0) {
		result.push_back(abs(var));

		return;
	}

	vector<int> & rhs = slg_nonterm_vec[var].rhs;

	for(int rhs_var : rhs) {
		expandSLG(rhs_var, slg_nonterm_vec, result);
	}
}

#endif
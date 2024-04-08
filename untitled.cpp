#IFNDEF RECOMPRESSION_DEFINITIONS_HPP
#DEFINE RECOMPRESSION_DEFINITIONS_HPP

struct RLSLPNonterm {
	char type;
	int first;
	int second;
	int explen;

	RLSLPNonterm(char type, int first, int second) : type(type), first(first), second(second) {}
	RLSLPNonterm() : type('0'), first(0), second(0), explen(0) {
	}
};

class RecompressionRLSLP {
public:
	vector<RLSLPNonterm> nonterm;
};

struct SLGNonterm {
	// rhs empty represents empty variable.
	vector<int> rhs;
	int vOcc;
	int LMS;
	int RMS;
	// {-1, -1} is empty
	pair<int, int> LR;
	pair<int, int> RR;

	SLGNonterm(vector<int> rhs) : vOcc(0), LMS(0), RMS(0), LR({-1, -1}), RR({-1, -1}), rhs(rhs) {

	}

	SLGNonterm() : vOcc(0), LMS(0), RMS(0), LR({-1, -1}), RR({-1, -1}) {

	}
};

class SLG {

public:
	SLG() {

	}
	SLG(vector<SLGNonterm> & nonterm) : nonterm(nonterm) {

	}
	vector<SLGNonterm> nonterm;

};

struct SLPNonterm {
	char type;
	int first;
	int second;

	SLPNonterm(char type, int first, int second) : type(type), first(first), second(second) {

	}
};

class InputSLP {
public:
	InputSLP() {

	}
	InputSLP(SLPNonterm nonterm) : nonterm(nonterm) {

	}
	vector<SLPNonterm> nonterm;
};

#ENDIF
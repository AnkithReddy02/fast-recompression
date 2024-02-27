struct NonTerminal {
    char type;
    int first;
    int second;
    int explen;

    NonTerminal(char type, int first, int second, int explen) {
    	this->type = type;
    	this->first = first;
    	this->second = second;
    	this->explen = explen;
    }
};

struct Node {
    
    // Variable/Non-Terminal
    int var;
    // [l, r)
    int l;
    int r;

    Node() {
        
    }
    Node(int var, int l, int r) : var(var), l(l), r(r) {

    }
};
#include <iostream>
#include <vector>
#include <array>
#include <stack>
#include <unordered_set>
#include <algorithm>
#include <assert.h>
#include <ctime>
#include "io.h"
#include "radix_sort.h"
// #include "testlib.h"
#include "definitions.hpp"

using namespace std;

vector<int> BComp(vector<int> arr, vector<NonTerminal>& grammar) {
    int count = 1;
    size_t n = arr.size();
    vector<array<int, 3>> tuples;

    for (int i = n - 2; i >= 0; i--) {
        if (arr[i] == arr[i + 1]) {
            count++;
        } else {
            if (count >= 2) {
                tuples.push_back({ arr[i + 1], count, i + 1 });
            }
            count = 1;
        }
    }

    if (count >= 2) {
        tuples.push_back({ arr[0], count, 0 });
    }

    if (tuples.size() == 0) {
        return arr;
    }

    radixSort(tuples);

    int maxElement = *max_element(arr.begin(), arr.end());
    int freshCharacter = maxElement + 1;
    tuples.push_back(tuples.back());
    tuples.back()[0]++;

    count = 1;

    for (int i = 1; i < tuples.size(); i++) {
        if (tuples[i][0] == tuples[i - 1][0] && tuples[i][1] == tuples[i - 1][1]) {
            count++;
        } else {
            for (int j = i - count; j <= i - 1; j++) {
                for (int k = tuples[j][2]; k < tuples[j][2] + tuples[j][1]; k++) {
                    arr[k] = freshCharacter;
                }
            }

            if (grammar.size() <= freshCharacter) {
                char type = '2';
                int first = tuples[i - 1][0];
                int second = tuples[i - 1][1];
                int explen = second * grammar[first].explen;
                grammar.push_back(NonTerminal(type, first, second, explen));
            }

            count = 1;
            freshCharacter++;
        }
    }

    tuples.pop_back();

    vector<int> res;
    for (int& x : arr) {
        if (res.empty() || res.back() != x) {
            res.push_back(x);
        }
    }

    return res;
}

vector<array<int, 2>> createAdjList(vector<int>& arr) {
    vector<array<int, 2>> adjList;
    int n = arr.size();

    for (int i = 0; i < n - 1; i++) {
        adjList.push_back({ max(arr[i], arr[i + 1]), min(arr[i], arr[i + 1]) });
    }

    return adjList;
}

// only adjList as input. *
array<unordered_set<int>, 2> createPartition(vector<int>& arr, vector<array<int, 2>>& adjList) {
    
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

    for (int i = 0; i < arr.size() - 1; i++) {
        LRPairsCount += (leftSet.find(arr[i]) != leftSet.end()) && (rightSet.find(arr[i + 1]) != rightSet.end());
        RLPairsCount += (rightSet.find(arr[i]) != rightSet.end()) && (leftSet.find(arr[i + 1]) != leftSet.end());
    }

    if (RLPairsCount < LRPairsCount) {
        swap(leftSet, rightSet);
    }

    return { rightSet, leftSet };
}

vector<int> PComp(vector<int> arr, vector<NonTerminal>& grammar) {
    // assert(arr.size() > 1);
    vector<array<int, 2>> adjList = createAdjList(arr);
    radixSort(adjList);

    auto partition = createPartition(arr, adjList);
    unordered_set<int> partition1 = partition[0];
    unordered_set<int> partition2 = partition[1];
    int n = arr.size();
    vector<array<int, 3>> compressablePairs;

    for (int i = 0; i < n - 1; i++) {
        if (partition1.find(arr[i]) != partition1.end() && partition2.find(arr[i + 1]) != partition2.end()) {
            compressablePairs.push_back({ arr[i], arr[i + 1], i });
        }
    }

    radixSort(compressablePairs);
    compressablePairs.push_back(compressablePairs.back());
    compressablePairs.back()[0]++;
    int count = 1;
    int freshCharacter = *max_element(arr.begin(), arr.end()) + 1;

    for (int i = 1; i < compressablePairs.size(); i++) {
        if (compressablePairs[i - 1][0] == compressablePairs[i][0] && compressablePairs[i - 1][1] == compressablePairs[i][1]) {
            count++;
        } else {
            if (grammar.size() <= freshCharacter) {
                char type = '1';
                int first = compressablePairs[i - 1][0];
                int second = compressablePairs[i - 1][1];
                int explen = grammar[first].explen + grammar[second].explen;
                grammar.push_back(NonTerminal(type, first, second, explen));
            }

            for (int j = i - count; j < i; j++) {
                arr[compressablePairs[j][2]] = freshCharacter;
                arr[compressablePairs[j][2] + 1] = freshCharacter;
            }

            freshCharacter++;
            count = 1;
        }
    }

    compressablePairs.pop_back();
    arr.push_back(arr.back());
    arr.back()++;
    count = 1;

    vector<int> res;
    for (int i = 1; i < arr.size(); i++) {
        if (arr[i] == arr[i - 1]) {
            count++;
        } else {
            if (count == 1) {
                res.push_back(arr[i - 1]);
            } else {
                count /= 2;
                for (int j = 0; j < count; j++) {
                    res.push_back(arr[i - 1]);
                }
                count = 1;
            }
        }
    }

    return res;
}

int extract(int node, int i, int left, int right, vector<NonTerminal>& grammar) {
    if (left > right || i < left || i > right)
        return 0;

    const NonTerminal& nt = grammar[node];

    if (nt.type == '0') {
        return nt.first;
    } else if (nt.type == '1') {
        const NonTerminal& left_nt = grammar[nt.first];
        const NonTerminal& right_nt = grammar[nt.second];

        return extract(nt.first, i, left, left + left_nt.explen - 1, grammar) +
            extract(nt.second, i, left + left_nt.explen, right, grammar);

    } else {
        const NonTerminal& child_nt = grammar[nt.first];
        int child_explen = child_nt.explen;

        int newLeft = left + ((i - left) / child_explen) * child_explen;
        int newRight = left + ((i - left) / child_explen) * child_explen + child_explen - 1;

        return extract(nt.first, i, newLeft, newRight, grammar);
    }

    return -1;
}

void initialize_nodes(int node, const int& i, int left, int right, stack<Node>& ancestors, vector<NonTerminal>& grammar, Node& v) {
    if (left > right || i < left || i > right)
        return;

    const NonTerminal& nt = grammar[node];

    if (left == i) {
        v = Node(node, left, right + 1);
        return;
    }

    ancestors.push(Node(node, left, right + 1));

    if (nt.type == '0') {
        return;
    } else if (nt.type == '1') {
        const NonTerminal& left_nt = grammar[nt.first];
        const NonTerminal& right_nt = grammar[nt.second];

        initialize_nodes(nt.first, i, left, left + left_nt.explen - 1, ancestors, grammar, v);
        initialize_nodes(nt.second, i, left + left_nt.explen, right, ancestors, grammar, v);

    } else {
        const NonTerminal& child_nt = grammar[nt.first];
        int child_explen = child_nt.explen;

        int newLeft = left + ((i - left) / child_explen) * child_explen;
        int newRight = left + ((i - left) / child_explen) * child_explen + child_explen - 1;

        initialize_nodes(nt.first, i, newLeft, newRight, ancestors, grammar, v);
    }

    return;
}


Node getLeftMostChild(Node v, vector<NonTerminal> & grammar) {
	char type = grammar[v.var].type;

	const NonTerminal& nt = grammar[v.var];

    int left = v.l;
    int right = v.r - 1;

	if(type == '1') {
		const NonTerminal& left_nt = grammar[nt.first];
    	const NonTerminal& right_nt = grammar[nt.second];
		return Node(nt.first, left, left + left_nt.explen);
	}
	else if(type == '2') {
		const NonTerminal& child_nt = grammar[nt.first];
    	int child_explen = child_nt.explen;
		return Node(nt.first, left, left + child_explen);
	}

	return Node();
}


int getChildCount(Node parent, vector<NonTerminal> & grammar) {

	NonTerminal & nt = grammar[parent.var];

	if(nt.type == '0') {
		return 1;
	}
	else if(nt.type == '1') {
		return 2;
	}
	else {
		return nt.second;
	}
}

int getChildIndex(Node parent, Node v, vector<NonTerminal> & grammar) {
	NonTerminal & nt = grammar[parent.var];

	if(nt.type == '0') {
		return 1;
	}
	else if(nt.type == '1') {
		if(parent.l == v.l) {
			return 1;
		}
		else {
			return 2;
		}
	}
	else {
		return (v.l - parent.l)/(v.r - v.l) + 1;
	}

}

Node getKthSibling(Node parent, Node v, int k, vector<NonTerminal> &grammar ) {
	NonTerminal & nt = grammar[parent.var];

	int left = parent.l;
	int right = parent.r;

	// cout << "328: " << parent.var << ' ' << left << ' ' << right << ' ' << k << endl;
	

	return Node(v.var, left + (k-1)*(v.r - v.l), left + (k-1)*(v.r - v.l) + (v.r - v.l));
}

Node replaceWithHighestStartingAtPosition(Node v, stack<Node> &ancestors, vector<NonTerminal> & grammar) {

	Node child = v;
	while(ancestors.empty() == false and ancestors.top().r == v.r) {
		child = ancestors.top();
		ancestors.pop();
	}

	Node ancestor = ancestors.top();

	// cout << "sz: " << ' ' << ancestor.var << endl;

	NonTerminal nt = grammar[ancestor.var];

	// cout << nt.type << endl;


	if(nt.type == '1') {
		return Node(nt.second, child.r, ancestor.r);
	}
	else {
		int childIndex = getChildIndex(ancestor, child, grammar);

		// cout << "353: " << childIndex << endl;

		return getKthSibling(ancestor, child, childIndex + 1, grammar);
	}

	return Node();

}

int LCE(Node v1, Node v2, int i, stack<Node> & v1_ancestors, stack<Node> & v2_ancestors, vector<NonTerminal> &grammar) {

	int exp_len_v1 = v1.r - v1.l;

	int exp_len_v2 = v2.r - v2.l;

	// cout << exp_len_v1 << ' ' << exp_len_v2 << endl;

	// cout << "NEW" << endl;
	// cout << v1.var << ' ' << v1.l << ' ' << v1.r << endl;
	// cout << v2.var << ' ' << v2.l << ' ' << v2.r << endl;
	// cout << "END" << endl;

	if(exp_len_v1 == 1 and exp_len_v1 == exp_len_v2 and v1.var != v2.var) {
		return v1.l - i;
	}
	else if(exp_len_v1 > exp_len_v2) {
		v1_ancestors.push(v1);
		v1 = getLeftMostChild(v1, grammar);
		// cout << "374" << ' ' << v1.var << ' ' << v1.l << ' ' << v1.r << endl;
	}
	else if(exp_len_v1 < exp_len_v2) {
		v2_ancestors.push(v2);
		v2 = getLeftMostChild(v2, grammar);
	}
	else if(v1.var != v2.var) {
		v2_ancestors.push(v2);
		v1_ancestors.push(v1);
		v1 = getLeftMostChild(v1, grammar);
		v2 = getLeftMostChild(v2, grammar);
	}
	else {

		Node v1_parent = v1_ancestors.top();
		Node v2_parent = v2_ancestors.top();

		int j1 = getChildIndex(v1_parent, v1, grammar);
		int j2 = getChildIndex(v2_parent, v2, grammar);


		int d1 = getChildCount(v1_parent, grammar);
		int d2 = getChildCount(v2_parent, grammar);

		int lambda = min(d1 - j1, d2 - j2);

		// cout << "OK" << endl;

		// cout << j1 << ' ' << d1 << endl;
		// cout << j2 << ' ' << d2 << endl;

		// cout << lambda << endl;

		// cout << "OK END" << endl;


		if(lambda <= 1) {
			v1 = replaceWithHighestStartingAtPosition(v1, v1_ancestors, grammar);
			v2 = replaceWithHighestStartingAtPosition(v2, v2_ancestors, grammar);
		}
		else {
			v1 = getKthSibling(v1_parent, v1, j1 + lambda, grammar);
			v2 = getKthSibling(v2_parent, v2, j2 + lambda, grammar);
		}


	}

	return LCE(v1, v2, i, v1_ancestors, v2_ancestors, grammar);;
}

bool hasNumberInRange(const vector<int>& arr, int minRange, int maxRange) {
    for (int num = minRange; num <= maxRange; ++num) {
        if (find(arr.begin(), arr.end(), num) == arr.end()) {
            return false;
        }
    }
    return true;
}

vector<int> generateRandomNumbers(int n, int minRange, int maxRange) {
    // Seed the random number generator
    srand(time(nullptr));

    vector<int> numbers;

    for (int i = 0; i < n; ++i) {
        numbers.push_back(rand() % (maxRange - minRange + 1) + minRange);
    }

    if(!hasNumberInRange(numbers, minRange, maxRange)) {
    	return generateRandomNumbers(n, minRange, maxRange);
    }

    return numbers;
}

void solve() {
	vector<int> arr;
    vector<NonTerminal> grammar;
    grammar.push_back(NonTerminal('0', 0, 0, 0));

    int x = 61;
    for(int i=96; i<=96+x-1; i++) {
    	grammar.push_back(NonTerminal('0', i, 0, 1));
    }
    
    arr = { 3, 1, 1, 1, 2, 3, 4, 2, 2, 2, 1, 2, 1, 2, 3, 4, 1, 1, 2, 3, 4, 2, 2, 2, 1, 2, 1, 2 ,3, 4, 4 };

    arr = generateRandomNumbers(900, 1, 60);

    // for(int i=0;i<200;i++) {
    // 	int x;

    // 	cin >> x;



    // 	arr.push_back(x);
    // }


    vector<int> arrCopy = arr;
    int n = arr.size();

    arr.push_back(*max_element(arr.begin(), arr.end()) + 1);

    // cout << arr << endl;

    // cout << endl;



    int count = 0;

    while (arr.size() > 1) {
        if (count % 2 == 0) {
            arr = BComp(arr, grammar);
        } else {
            arr = PComp(arr, grammar);
        }
        count++;
    }


    int explen = grammar.back().explen;

    // grammar.push_back(NonTerminal('0', 100, 0, 1));
    // grammar.push_back(NonTerminal('1', grammar.size()-2, grammar.size()-1, explen + 1));

    vector<int> res;

    for (int i = 0; i < n; i++) {
        res.push_back(extract(grammar.size() - 1, i, 0, grammar.back().explen - 1, grammar) - 95);
    }

    // cout << res << endl;

    // assert(arrCopy == res);

    // cout << n << endl;

    // cout << arrCopy[43] << ' ' << arrCopy[116] << endl;

    // for(int i=43;i<=50;i++) {
    // 	cout << arrCopy[i] << ' ';
    // }

    // cout << endl;


    // for(int i=116; i<= 160; i++) {
    // 	cout << arrCopy[i] << ' ';
    // }

    // cout << endl;


    for(int i=0; i<n; i++) {
    	for(int j=i+1; j<n; j++) {

    		// if(i!=7 or j!=11) continue;

    		// cout << i << ' ' << j << endl;


    		Node v1, v2;
		    stack<Node> v1_ancestors, v2_ancestors;
		    // v1_ancestors.push(Node(grammar.size()-1, 0, 33));
		    // v2_ancestors.push(Node(grammar.size()-1, 0, 33));
		    initialize_nodes(grammar.size() - 1, i, 0, grammar.back().explen - 1, v1_ancestors, grammar, v1);
		    initialize_nodes(grammar.size() - 1, j, 0, grammar.back().explen - 1, v2_ancestors, grammar, v2);


		    // cout << v1.var << ' ' << v1.l << ' ' << v1.r << endl;
		    // cout << v2.var << ' ' << v2.l << ' ' << v2.r << endl;

		    int res1 = LCE(v1, v2, i, v1_ancestors, v2_ancestors, grammar);

		    int res2 = 0;

		    int ii = i;
		    int jj = j;

		    while(jj < n && arrCopy[ii] == arrCopy[jj]) {
		    	res2++;

		    	jj++;
		    	ii++;
		    }

		    // cout << i << ' ' << j << ' ' << res1 << ' ' << res2 << endl;
		    // assert(res1==res2);

		    if(res1 != res2) {
		    	cout << "ERROR" << endl;
		    	exit(0);
		    }


    	}
    }
}
int main(int argc, char* argv[]) {
    // registerGen(argc, argv, 1);

    for(int i=1; i<=10;i++) {
    	solve();
    }

    return 0;
}

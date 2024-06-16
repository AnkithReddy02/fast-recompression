#ifndef TEST_QUERIES_HPP
#define TEST_QUERIES_HPP

void get_random_queries(c_size_t text_size,
    space_efficient_vector<pair<c_size_t, c_size_t>> &pairs) {
    // pairs.reserve(1000000); 

    random_device rd; 
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> distrib(0, text_size - 1);

    // Generate 4096 pairs
    for(c_size_t i = 0; i < 1000000; ++i) {
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

        if(i == j) continue;
        if(i > j) swap(i, j);

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

void test_queries(RecompressionRLSLP *recompression_rlslp, const string &input_slp) {
    space_efficient_vector<RLSLPNonterm> & rlslp_nonterm_vec = recompression_rlslp->nonterm;
    space_efficient_vector<pair<c_size_t, c_size_t>> random_queries;

    typedef simple_slp<std::uint8_t, c_size_t> simple_slp_type;
    simple_slp_type *slp = new simple_slp_type(input_slp);

    get_random_queries(rlslp_nonterm_vec.back().explen, random_queries);

    auto start_time = std::chrono::high_resolution_clock::now();

    // for(auto x : random_queries) {
    for(c_size_t k = 0; k < random_queries.size(); ++k) {
        const auto& x = random_queries[k];

        c_size_t i = x.first;
        c_size_t j = x.second;

        if(i == j) continue;
        if(i > j) swap(i, j);

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
        c_size_t res2 = slp->lce(i, j);

        if(res1 != res2) {
            cerr << "Error: LCE Query didn't match for i: " << i << ", j: " << j << endl;
            exit(1);
        }


    }

    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

    cout << endl <<  "LCE Queries Passed!" << endl;
    // Output the duration
    cout << "Time taken for LCE Queries: " << duration_seconds << " seconds" << endl << endl;

    delete slp;
}

#endif
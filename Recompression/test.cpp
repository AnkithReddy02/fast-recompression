#include <bits/stdc++.h>
#include "hash_table.hpp"
#include "space_efficient_vector.hpp"
#include "utils.hpp"
#include "typedefs.hpp"
using namespace std;

template<>
std::uint64_t get_hash(const c_size_t &x) {
    return (std::uint64_t)x * (std::uint64_t)1232545666776865;
}

int main() {
    hash_table<c_size_t, c_size_t, c_size_t> m;

    m.insert(0, 1);
    m.insert(1, 2);
    m.insert(2, 3);
    m.insert(3, 4);
    m.insert(4, 5);


    for(size_t i=0; i<m.size(); i++) {
        cout << m.get(i).first << ' ' << m.get(i).second << endl;
    }

}
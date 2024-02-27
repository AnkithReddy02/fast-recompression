#ifndef IO_H
#define IO_H

#include <iostream>
#include <vector>
#include <array>
#include <unordered_set>
#include <unordered_map>

template <typename T, size_t N>
std::ostream & operator<<(std::ostream & os, const std::array<T, N> & arr) {
    for(const auto & x : arr) {
        os << x << ' ';
    }
    return os;
}

template <typename T>
std::ostream & operator<<(std::ostream & os, const std::unordered_set<T> & arr) {
    for(const auto & x : arr) {
        os << x << ' ';
    }
    return os;
}

template <typename T, typename S, size_t N>
std::ostream & operator<<(std::ostream & os, const std::unordered_map<T, std::array<S, N>> & arr) {
    for(const auto & x : arr) {
        os << x.first << " : [";
        for(const auto & y : x.second) {
            os << y << ' ';
        }
        os << ']';
        os << std::endl;
    }
    return os;
}

template <typename T, size_t N>
std::ostream & operator<<(std::ostream & os, const std::vector<std::array<T, N>> & arr) {
    for(const auto & x : arr) {
        os << x << std::endl;
    }
    return os;
}

template <typename T>
std::ostream & operator<<(std::ostream & os, const std::vector<T> & arr) {
    for(const auto & x : arr) {
        os << x << ' ';
    }
    return os;
}

#endif // IO_H

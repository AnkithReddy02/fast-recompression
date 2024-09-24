/**
 * @file    uint40.hpp
 * @section LICENCE
 *
 * This file is part of Lazy-AVLG v0.1.0
 * See: https://github.com/dominikkempa/lz77-to-slp
 *
 * Copyright (C) 2021
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#ifndef __INT40_HPP_INCLUDED
#define __INT40_HPP_INCLUDED

#include <cstdint>
#include <limits>

class int40 {
  private:
    std::uint32_t low;
    std::int8_t high;

  public:
    int40() {}
    int40(std::uint32_t l, std::int8_t h) : low(l), high(h) {}
    int40(const int40& a) : low(a.low), high(a.high) {}
    int40(const std::int32_t& a) : low(a), high((a < 0) ? -1 : 0) {}
    int40(const std::uint32_t& a) : low(a), high(0) {}
    int40(const std::int64_t& a) :
      low(a & 0xFFFFFFFFL), high((a >> 32) & 0xFF) {}
    int40(const std::uint64_t& a) :
      low(a & 0xFFFFFFFF), high((a >> 32) & 0xFF) {}

    inline operator int64_t() const {
        return (((std::int64_t)high) << 32) | (std::uint64_t)low;
    }

    // Relational
    template<typename T>
    inline bool operator == (const T &b) const {
        return int64_t(*this) == int64_t(b);
    }
    template<typename T>
    inline bool operator != (const T &b) const {
        return int64_t(*this) != int64_t(b);
    }
    template<typename T>
    inline bool operator < (const T &b) const {
        return int64_t(*this) < int64_t(b);
    }
    template<typename T>
    inline bool operator <= (const T &b) const {
        return int64_t(*this) <= int64_t(b);
    }
    template<typename T>
    inline bool operator > (const T &b) const {
        return int64_t(*this) > int64_t(b);
    }
    template<typename T>
    inline bool operator >= (const T &b) const {
        return int64_t(*this) >= int64_t(b);
    }

    // Assignment
    template<typename T>
    inline int40& operator += (const T &b) {
        int64_t temp(*this);
        temp += int64_t(b);
        *this = temp;
        return *this;
    }
    inline int40& operator -= (const int40 &b) {
        *this = *this - b;
        return *this;
    }
    inline int40& operator *= (const int40 &b) {
        *this = *this * b;
        return *this;
    }
    inline int40& operator /= (const int40 &b) {
        *this = *this / b;
        return *this;
    }

    // Arithmetic
    template<typename T>
    inline int40 operator + (const T &b) {
        int64_t result = int64_t(*this) + int64_t(b);
        return int40(result);
    }
    template<typename T>
    inline int40 operator - (const T &b) {
        int64_t result = int64_t(*this) - int64_t(b);
        return int40(result);
    }
    template<typename T>
    inline int40 operator * (const T &b) {
        int64_t result = int64_t(*this) * int64_t(b);
        return int40(result);
    }
    template<typename T>
    inline int40 operator / (const T &b) {
        int64_t result = int64_t(*this) / int64_t(b);
        return int40(result);
    }

    // inline long double operator + (const long doube &b) {
    //     long double result = int64_t(*this) + b;
    //     return result;
    // }
    
    // Unary Arithmetic
    // Prefix
    int40& operator ++ () {
        *this += 1;
        return *this;
    }
    // Postfix
    int40 operator ++ (int) {
        int40 temp = *this;
        ++(*this);
        return temp;
    }
    // Prefix
    int40& operator -- () {
        *this -= 1;
        return *this;
    }
    // Postfix
    int40 operator -- (int) {
        int40 temp = *this;
        --(*this);
        return temp;
    }


} __attribute__((packed));

namespace std {

template<>
struct is_signed<int40> {
  public:
    static const bool value = true;
};

template<>
class numeric_limits<int40> {
  public:
    
    /*
    ***********************************************
       HIGH  |             LOW
    -----------------------------------------------   
    10000000 | 00000000000000000000000000000000
    ***********************************************
    */
    static int40 min() {
      return int40(std::numeric_limits<std::uint32_t>::min(), std::numeric_limits<std::int8_t>::min());
    }

    /*
    ***********************************************
       HIGH  |             LOW
    -----------------------------------------------   
    01111111 | 11111111111111111111111111111111
    ***********************************************
    */
    static int40 max() {
      return int40(std::numeric_limits<std::uint32_t>::max(), std::numeric_limits<std::int8_t>::max());
    }
};

}  // namespace std

#endif  // __INT40_HPP_INCLUDED

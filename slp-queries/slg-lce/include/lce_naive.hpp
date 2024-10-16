/**
 * @file    lce_naive.hpp
 * @section LICENCE
 *
 * Copyright (C) 2024
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

#ifndef __LCE_NAIVE_HPP_INCLUDED
#define __LCE_NAIVE_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>

template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
struct lce_naive {
  static_assert(sizeof(char_type) <= sizeof(text_offset_type),
      "Error: sizeof(char_type) > sizeof(text_offset_type)!");

  private:

    char_type *m_text;
    std::uint64_t m_text_length;

  public:

    lce_naive(const std::string text_filename) {
      m_text_length = utils::file_size(text_filename) / sizeof(char_type);
      m_text = utils::allocate_array<char_type>(m_text_length);
      utils::read_from_file<char_type>(m_text, m_text_length, text_filename);
    }

    ~lce_naive() {
      utils::deallocate(m_text);
    }

    std::uint64_t lce(
        const std::uint64_t i,
        const std::uint64_t j) const {
      if (i == j)
        return m_text_length - i;

      std::uint64_t ans = 0;
      while (i + ans < m_text_length && j + ans < m_text_length &&
          m_text[i + ans] == m_text[j + ans]) ++ans;
      return ans;
    }
};

#endif  // __LCE_NAIVE_HPP_INCLUDED


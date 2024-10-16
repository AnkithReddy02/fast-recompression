/**
 * @file    simple_slp.hpp
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

#ifndef __SIMPLE_SLP_HPP_INCLUDED
#define __SIMPLE_SLP_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>

#include "utils/hash_table.hpp"
#include "utils/karp_rabin_hashing.hpp"
#include "utils/space_efficient_vector.hpp"
#include "utils/cache.hpp"
#include "utils/packed_pair.hpp"
#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"


//=============================================================================
// Class used to represent the lazy AVL grammar. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct simple_slp;

//=============================================================================
// Class used to represent nonterminal. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct nonterminal {
  private:

    //=========================================================================
    // Declare types.
    //=========================================================================
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;
    typedef text_offset_type ptr_type;
    typedef simple_slp<char_type, text_offset_type> slp_type;
    typedef packed_pair<text_offset_type, text_offset_type> pair_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    std::uint8_t m_height;
    std::uint8_t m_truncated_exp_len;
    ptr_type m_left_p;
    ptr_type m_right_p;

  public:

    //=========================================================================
    // Constructors.
    //=========================================================================
    nonterminal();
    nonterminal(const char_type);
    nonterminal(const std::uint8_t, const std::uint8_t,
      const ptr_type, const ptr_type);

    //=========================================================================
    // Access methods.
    //=========================================================================
    std::uint64_t get_height() const;
    std::uint64_t get_truncated_exp_len() const;
    char_type get_char() const;
    ptr_type get_left_p() const;
    ptr_type get_right_p() const;

    //=========================================================================
    // Key methods.
    //=========================================================================
    std::uint64_t write_expansion(const ptr_type, char_type * const,
        const slp_type * const) const;
    std::uint64_t get_prefix_kr_hash(const ptr_type,
        const std::uint64_t, const slp_type * const) const;

    //=========================================================================
    // Mostly unused.
    //=========================================================================
    char_type access_symbol(const ptr_type,
        std::uint64_t, const slp_type * const) const;
} __attribute__((packed));

//=============================================================================
// Hash functions of the appropriate type.
// Used in the hash table used to prune the grammar.
//=============================================================================
typedef const uint40 get_hash_ptr_type;

template<>
std::uint64_t get_hash(const get_hash_ptr_type &x) {
  return (std::uint64_t)x * (std::uint64_t)29996224275833;
}

template<>
std::uint64_t get_hash(const std::uint64_t &x) {
  return (std::uint64_t)x * (std::uint64_t)4972548694736365;
}

//=============================================================================
// Implementation of the slp class.
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
struct simple_slp {
  static_assert(sizeof(char_type) <= sizeof(text_offset_type),
      "Error: sizeof(char_type) > sizeof(text_offset_type)!");

  private:

    //=========================================================================
    // Declare typedefs.
    //=========================================================================
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;
    typedef text_offset_type ptr_type;
    typedef packed_pair<text_offset_type, text_offset_type> pair_type;
    typedef packed_pair<text_offset_type, std::uint64_t> hash_pair_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    nonterminal_type *m_nonterminals;
    pair_type *m_long_exp_len;
    hash_pair_type *m_long_exp_hashes;
    cache<ptr_type, std::uint64_t> *m_kr_hash_cache;
    char_type *m_snippet;

    std::uint64_t m_nonterminals_size;
    std::uint64_t m_long_exp_len_size;
    std::uint64_t m_long_exp_hashes_size;
    std::uint64_t m_text_length;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    simple_slp(const std::string slp_filename) {
      const std::uint64_t cache_size = (1 << 14);

      // Initialize Karp-Rabin fingerprints.
      karp_rabin_hashing::init();

      // Initialize KR hashing.
      m_kr_hash_cache = new cache<ptr_type, std::uint64_t>(cache_size);
      m_snippet = utils::allocate_array<char_type>(256);

      // Initialize the grammar reader.
      const std::uint64_t bufsize = (1 << 19);
      typedef async_stream_reader<text_offset_type> reader_type;
      reader_type *reader = new reader_type(slp_filename, bufsize, 4);

      // Compute the number of nonterminals.
      const std::uint64_t filesize = utils::file_size(slp_filename);
      std::uint64_t n_nonterminals = filesize / sizeof(pair_type);

      // Allocate array for nonterminals.
      m_nonterminals = utils::allocate_array<nonterminal_type>(n_nonterminals);
      m_nonterminals_size = 0;
      
      // Read nonterminals and compute number of long nonterminals.
      std::uint64_t n_long_exp = 0;
      for (std::uint64_t i = 0; i < n_nonterminals; ++i) {
        const std::uint64_t val1 = reader->read();
        const std::uint64_t val2 = reader->read();
        if (val1 == 0) {
          nonterminal_type nonterm((char_type)val2);
          add_nonterminal(nonterm);
        } else {
          const std::uint64_t left_p = val1 - 1;
          const std::uint64_t right_p = val2 - 1;
          const ptr_type nonterm_p = add_nonterminal(left_p, right_p);
          const nonterminal_type &nonterm = get_nonterminal(nonterm_p);
          if (nonterm.get_truncated_exp_len() == 255)
            ++n_long_exp;
        }
      }

      // Clean up.
      delete reader;

      // Allocate expansion length and KR hash for long nonterminals.
      m_long_exp_len_size = 0;
      m_long_exp_hashes_size = 0;
      m_long_exp_len = utils::allocate_array<pair_type>(n_long_exp);
      m_long_exp_hashes = utils::allocate_array<hash_pair_type>(n_long_exp);

      // Compute expansion length and KR hashes for long nonterminals.
      for (std::uint64_t i = 0; i < n_nonterminals; ++i) {
        const nonterminal_type &nonterm = get_nonterminal(i);
        const std::uint64_t truncated_exp_len = nonterm.get_truncated_exp_len();
        std::uint64_t exp_len = truncated_exp_len;
        if (truncated_exp_len == 255) {

          // Update list of expansions.
          const ptr_type left_p = nonterm.get_left_p();
          const ptr_type right_p = nonterm.get_right_p();
          const std::uint64_t left_exp_len = get_exp_len(left_p);
          const std::uint64_t right_exp_len = get_exp_len(right_p);
          exp_len = left_exp_len + right_exp_len;
          m_long_exp_len[m_long_exp_len_size++] = pair_type(i, exp_len);

          // Update list of hashes.
          const std::uint64_t left_hash = get_kr_hash(left_p);
          const std::uint64_t right_hash = get_kr_hash(right_p);
          const std::uint64_t kr_hash = karp_rabin_hashing::concat(
              left_hash, right_hash, right_exp_len);
          m_long_exp_hashes[m_long_exp_hashes_size++] =
            hash_pair_type(i, kr_hash);
        }
      }

      // Compute text length.
      {
        const std::uint64_t root_id = get_root_id();
        m_text_length = get_exp_len(root_id);
      }
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~simple_slp() {
      utils::deallocate(m_long_exp_hashes);
      utils::deallocate(m_long_exp_len);
      utils::deallocate(m_nonterminals);
      utils::deallocate(m_snippet);
      delete m_kr_hash_cache;
    }

    //=========================================================================
    // Return the number of nonterminals.
    //=========================================================================
    std::uint64_t number_of_nonterminals() const {
      return m_nonterminals_size;
    }

    //=========================================================================
    // Return the pointer to the root.
    //=========================================================================
    ptr_type get_root_id() const {
      return m_nonterminals_size - 1;
    }

    //=========================================================================
    // Returns true if the given nonterminal is a liteval (i.e., exp len = 1).
    //========================================================================
    bool is_literal(const ptr_type p) const {
      const nonterminal_type &nonterm = get_nonterminal(p);
      const std::uint64_t height = nonterm.get_height();
      return (height == 0);
    }

    //=========================================================================
    // Return the total length of right-hand sides of all expansions.
    //=========================================================================
    std::uint64_t total_rhs_length() const {
      std::uint64_t ret = 0;
      for (std::uint64_t i = 0; i < number_of_nonterminals(); ++i) {
        if (is_literal(i))
          ret += 1;
        else ret += 2;
      }
      return ret;
    }

  public:

    //=========================================================================
    // Gives access to a given nonterminal.
    //=========================================================================
    const nonterminal_type& get_nonterminal(const ptr_type id) const {
      return m_nonterminals[(std::uint64_t)id];
    }

    //=========================================================================
    // Return the expansion length of a given nonterminal.
    //=========================================================================
    std::uint64_t get_exp_len(const ptr_type id) const {
      const nonterminal_type &nonterm = get_nonterminal(id);
      std::uint64_t exp_len = 0;
      const std::uint64_t truncated_exp_len = nonterm.get_truncated_exp_len();
      if (truncated_exp_len == 255) {

        // Binary search in m_long_exp_len.
        std::uint64_t beg = 0;
        std::uint64_t end = m_long_exp_len_size;
        while (beg + 1 < end) {
          const std::uint64_t mid = (beg + end) / 2;
          const std::uint64_t mid_id = m_long_exp_len[mid].first;
          if (mid_id <= (std::uint64_t)id)
            beg = mid;
          else end = mid;
        }
        exp_len = m_long_exp_len[beg].second;
      } else exp_len = truncated_exp_len;
      return exp_len;
    }

    //=========================================================================
    // Return the Karp-Rabin hash of a given nonterminal.
    //=========================================================================
    std::uint64_t get_kr_hash(const ptr_type id) const {

      // Check, if the value is in cache.
      const std::uint64_t *cache_ret =
        m_kr_hash_cache->lookup(id);
      if (cache_ret != NULL)
        return *cache_ret;

      // Obtain/compute the hash.
      const nonterminal_type &nonterm = get_nonterminal(id);
      const std::uint64_t truncated_exp_len =
        nonterm.get_truncated_exp_len();
      std::uint64_t kr_hash = 0;
      if (truncated_exp_len < 255) {

        // Recompute the hash from scratch.
        (void) nonterm.write_expansion(id, m_snippet, this);
        kr_hash = karp_rabin_hashing::hash_string<char_type>(
            m_snippet, truncated_exp_len);
      } else {

        // Binary search in m_long_exp_hashes.
        std::uint64_t beg = 0;
        std::uint64_t end = m_long_exp_hashes_size;
        while (beg + 1 < end) {
          const std::uint64_t mid = (beg + end) / 2;
          const std::uint64_t mid_id = m_long_exp_hashes[mid].first;
          if (mid_id <= id)
            beg = mid;
          else end = mid;
        }
        kr_hash = m_long_exp_hashes[beg].second;
      }

      // Update cache.
      m_kr_hash_cache->insert(id, kr_hash);

      // Return the result.
      return kr_hash;
    }

    //=========================================================================
    // Add a nonterminal.
    //=========================================================================
    ptr_type add_nonterminal(const nonterminal_type &nonterm) {
      const ptr_type new_nonterm_p = number_of_nonterminals();
      m_nonterminals[m_nonterminals_size++] = nonterm;

      // Return the id of the nonterminal.
      return new_nonterm_p;
    }

    //=========================================================================
    // Add a new binary nonterminal.
    //=========================================================================
    ptr_type add_nonterminal(
        const ptr_type left_p,
        const ptr_type right_p) {
      const nonterminal_type &left = get_nonterminal(left_p);
      const nonterminal_type &right = get_nonterminal(right_p);

      // Compute values for the new nonterminal.
      const std::uint64_t truncated_left_exp_len =
        left.get_truncated_exp_len();
      const std::uint64_t truncated_right_exp_len =
        right.get_truncated_exp_len();
      const std::uint64_t truncated_exp_len =
        std::min((std::uint64_t)255,
            truncated_left_exp_len + truncated_right_exp_len);
      const std::uint8_t left_height = left.get_height();
      const std::uint8_t right_height = right.get_height();
      const std::uint8_t height = std::max(left_height, right_height) + 1;

      // Create and add new nonterminal.
      const ptr_type new_nonterm_p = number_of_nonterminals();
      nonterminal_type new_nonterm(height, truncated_exp_len, left_p, right_p);
      m_nonterminals[m_nonterminals_size++] = new_nonterm;

      // Return the ptr to the new nonterminal.
      return new_nonterm_p;
    }

  public:

    //=========================================================================
    // Check if given text = expansion of the grammar.
    //=========================================================================
    bool compare_expansion_to_text(
        const std::string text_filename,
        const std::uint64_t text_length) {
      
      static const std::uint64_t bufsize = (1 << 19);
      typedef async_stream_reader<char_type> reader_type;
      reader_type *reader = new reader_type(text_filename, bufsize, 4);

#ifdef NDEBUG
      static const std::uint64_t max_prefix_length = (10 << 20);
#else
      static const std::uint64_t max_prefix_length = (1 << 19);
#endif

      bool ok = true;
      for (std::uint64_t i = 0; i <
          std::min(max_prefix_length, text_length); ++i) {
        char_type text_i = reader->read();
        if (access_symbol(i) != text_i) {
          ok = false;
          break;
        }
      }

      // Clean up.
      delete reader;

      // Return the answer.
      return ok;
    }

    //=========================================================================
    // Return RAM use.
    //=========================================================================
    std::uint64_t ram_use() const {
      const std::uint64_t m_nonterminals_ram_use =
        number_of_nonterminals() * sizeof(nonterminal_type);
      const std::uint64_t m_long_exp_len_ram_use =
        m_long_exp_len_size * sizeof(pair_type);
      const std::uint64_t m_long_exp_hashes_ram_use =
        m_long_exp_hashes_size * sizeof(hash_pair_type);
      const std::uint64_t m_kr_hash_cache_ram_use =
        m_kr_hash_cache->ram_use();

      const std::uint64_t total =
        m_nonterminals_ram_use + 
        m_long_exp_hashes_ram_use +
        m_long_exp_len_ram_use +
        m_kr_hash_cache_ram_use;
      return total;
    }

    //=========================================================================
    // Print statistics.
    //=========================================================================
    void print_stats() const {

      // Print RAM usage breakdown.
      const std::uint64_t m_nonterminals_ram_use =
        number_of_nonterminals() * sizeof(nonterminal_type);
      const std::uint64_t m_long_exp_len_ram_use =
        m_long_exp_len_size * sizeof(pair_type);
      const std::uint64_t m_long_exp_hashes_ram_use =
        m_long_exp_hashes_size * sizeof(hash_pair_type);
      const std::uint64_t m_kr_hash_cache_ram_use =
        m_kr_hash_cache->ram_use();
      const std::uint64_t total =
        m_nonterminals_ram_use + 
        m_long_exp_hashes_ram_use +
        m_long_exp_len_ram_use +
        m_kr_hash_cache_ram_use;

      fprintf(stderr, "RAM use:\n");
      fprintf(stderr, "  m_nonterminals: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_nonterminals_ram_use) / (1 << 20),
          (100.L * m_nonterminals_ram_use) / total);
      fprintf(stderr, "  m_long_exp_len: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_long_exp_len_ram_use) / (1 << 20),
          (100.L * m_long_exp_len_ram_use) / total);
      fprintf(stderr, "  m_long_exp_hashes: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_long_exp_hashes_ram_use) / (1 << 20),
          (100.L * m_long_exp_hashes_ram_use) / total);
      fprintf(stderr, "  m_kr_hash_cache: %.2LfMiB (%.2Lf%%)\n",
          (1.L * m_kr_hash_cache_ram_use) / (1 << 20),
          (100.L * m_kr_hash_cache_ram_use) / total);
      fprintf(stderr, "Total: %.2LfMiB (%.2Lf%%)\n",
          (1.L * total) / (1 << 20),
          (100.L * total) / total);
    }

    //=========================================================================
    // Return symbol text[i].
    //=========================================================================
    char_type access_symbol(const std::uint64_t i) const {
      const ptr_type root_id = get_root_id();
      const nonterminal_type &nonterm = get_nonterminal(root_id);
      return nonterm.access_symbol(root_id, i, this);
    }

    //=========================================================================
    // Return KR hash of text[0..prefix_length).
    //=========================================================================
    std::uint64_t prefix_kr_hash(const std::uint64_t prefix_length) const {
      const ptr_type nonterm_p = get_root_id();
      const nonterminal_type &nonterm = get_nonterminal(nonterm_p);
      const std::uint64_t kr_hash =
        nonterm.get_prefix_kr_hash(nonterm_p, prefix_length, this);
      return kr_hash;
    }

    //=========================================================================
    // Return KR hash of text[i..j).
    //=========================================================================
    std::uint64_t substring_kr_hash(
        const std::uint64_t i,
        const std::uint64_t j) const {
      const std::uint64_t x = prefix_kr_hash(i);
      const std::uint64_t y = prefix_kr_hash(j);
      const std::uint64_t len = j - i;
      return karp_rabin_hashing::unconcat(x, y, len);
    }

    //=========================================================================
    // Return LCE(i, j). Complexity: O(n).
    //=========================================================================
    std::uint64_t lce_naive(
        const std::uint64_t i,
        const std::uint64_t j) const {

      // Handle the simple case.
      if (i == j)
        return m_text_length - i;

      // Compute LCE with a linear scan.
      const std::uint64_t max_lce = m_text_length - std::max(i, j);
      std::uint64_t ret = 0;
      while (ret < max_lce &&
          access_symbol(i + ret) == access_symbol(j + ret))
        ++ret;

      // Return the result.
      return ret;
    }

    //=========================================================================
    // Return LCE(i, j). Complexity: O(log LCE(i, j)).
    //=========================================================================
    std::uint64_t lce_using_kr_hash(
        const std::uint64_t i,
        const std::uint64_t j) const {

      // Handle the simple case.
      if (i == j)
        return m_text_length - i;

      // Compute a small upper bound on the answer.
      const std::uint64_t max_lce = m_text_length - std::max(i, j);
      std::uint64_t low = 0;
      std::uint64_t high = 1;
      while (high < max_lce &&
          substring_kr_hash(i, i + high) == substring_kr_hash(j, j + high)) {
        low = high;
        high = std::min(high << 1, max_lce);
      }

      // Finish using binary search.
      while (low != high) {

        // Invariant: answer is in the interval [log..high].
        std::uint64_t mid = (low + high + 1) / 2;
        if (substring_kr_hash(i, i + mid) == substring_kr_hash(j, j + mid))
          low = mid;
        else high = mid - 1;
      }

      // Return the result.
      return low;
    }

    //=========================================================================
    // Return LCE(i, j). Complexity: O(log LCE(i, j)). Uses hybrid approach,
    // where we first perform a linear scan, and then use binary search.
    //=========================================================================
    std::uint64_t lce(
        const std::uint64_t i,
        const std::uint64_t j) const {

      // Handle the simple case.
      if (i == j)
        return m_text_length - i;

      // Compute LCE with a linear scan.
      const std::uint64_t lce_limit = 32;
      const std::uint64_t lce_max = m_text_length - std::max(i, j);
      std::uint64_t lce_scan = 0;
      while (lce_scan < lce_limit && lce_scan < lce_max &&
          access_symbol(i + lce_scan) == access_symbol(j + lce_scan))
        ++lce_scan;

      // Return if linear scan computed the answer.
      if (lce_scan < lce_limit || lce_scan == lce_max)
        return lce_scan;

      // Compute a small upper bound on the answer.
      // Note that at this point: 0 < lce_scan < lce_max.
      std::uint64_t low = lce_scan;
      std::uint64_t high = lce_scan + 1;
      while (high < lce_max &&
          substring_kr_hash(i, i + high) == substring_kr_hash(j, j + high)) {
        low = high;
        high = std::min(high << 1, lce_max);
      }

      // Finish using binary search.
      while (low != high) {

        // Invariant: answer is in the interval [log..high].
        std::uint64_t mid = (low + high + 1) / 2;
        if (substring_kr_hash(i, i + mid) == substring_kr_hash(j, j + mid))
          low = mid;
        else high = mid - 1;
      }

      // Return the result.
      return low;
    }
};

//=============================================================================
// Default constructor.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal()
  : m_height(0),
    m_truncated_exp_len(1),
    m_left_p(std::numeric_limits<text_offset_type>::max()),
    m_right_p(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for a nonterminal expanding to a single symbol.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(const char_type c)
  : m_height(0),
    m_truncated_exp_len(1),
    m_left_p((text_offset_type)((std::uint64_t)c)),
    m_right_p(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for non-single-symbol nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(
      const std::uint8_t height,
      const std::uint8_t truncated_exp_len,
      const ptr_type left_p,
      const ptr_type right_p)
  : m_height(height),
    m_truncated_exp_len(truncated_exp_len),
    m_left_p(left_p),
    m_right_p(right_p) {}

//=============================================================================
// Get nonterminal height.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>::get_height() const {
  return (std::uint64_t)m_height;
}

//=============================================================================
// Get nonterminal expansion length.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>
::get_truncated_exp_len() const {
  return (std::uint64_t)m_truncated_exp_len;
}

//=============================================================================
// Return the char stored in a nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
char_type nonterminal<char_type, text_offset_type>::get_char() const {
  return (char_type)((std::uint64_t)m_left_p);
}

//=============================================================================
// Get nonterminal left ptr.
//=============================================================================
template<typename char_type, typename text_offset_type>
text_offset_type nonterminal<char_type, text_offset_type>::get_left_p() const {
  return m_left_p;
}

//=============================================================================
// Get nonterminal left ptr.
//=============================================================================
template<typename char_type, typename text_offset_type>
text_offset_type nonterminal<char_type, text_offset_type>::get_right_p() const {
  return m_right_p;
}

//=============================================================================
// Write the expansion into the given array.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>::write_expansion(
    const ptr_type x_p,
    char_type * const text,
    const slp_type * const g) const {
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t x_height = x.get_height();
  if (x_height == 0) {
    const char_type my_char = x.get_char();
    text[0] = my_char;
    return 1;
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    const std::uint64_t x_left_exp_len =
      x_left.write_expansion(x_left_p, text, g);
    const std::uint64_t x_right_exp_len =
      x_right.write_expansion(x_right_p, text + x_left_exp_len, g);
    return x_left_exp_len + x_right_exp_len;
  }
}

//=============================================================================
// Get the KR hash of an expansion prefix.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>::get_prefix_kr_hash(
    const ptr_type x_p,
    const std::uint64_t prefix_length,
    const slp_type * const g) const {

  // Handle special case.
  if (prefix_length == 0)
    return 0;

  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint64_t x_exp_length = g->get_exp_len(x_p);
  if (prefix_length == x_exp_length)
    return g->get_kr_hash(x_p);
  else {

    // We are guaranteed that height of x is > 0.
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);

    const std::uint64_t x_left_exp_len = g->get_exp_len(x_left_p);
    if (prefix_length <= x_left_exp_len)
      return x_left.get_prefix_kr_hash(x_left_p, prefix_length, g);
    else {
      const std::uint64_t x_left_kr_hash =
        x_left.get_prefix_kr_hash(x_left_p, x_left_exp_len, g);
      const std::uint64_t right_prefix = prefix_length - x_left_exp_len;
      const std::uint64_t x_right_kr_hash =
        x_right.get_prefix_kr_hash(x_right_p, right_prefix, g);
      return karp_rabin_hashing::concat(x_left_kr_hash,
          x_right_kr_hash, right_prefix);
    }
  }
}

//=============================================================================
// Access i-th leftmost symbol in the expansion (i = 0,1,..).
//=============================================================================
template<typename char_type, typename text_offset_type>
char_type nonterminal<char_type, text_offset_type>::access_symbol(
    ptr_type x_p,
    const std::uint64_t i,
    const slp_type * const g) const {
  const nonterminal_type x = g->get_nonterminal(x_p);
  const std::uint64_t height = x.get_height();
  if (height == 0)
    return get_char();
  ptr_type x_left_p = x.get_left_p();
  ptr_type x_right_p = x.get_right_p();
  const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
  const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
  const std::uint64_t left_exp_len = g->get_exp_len(x_left_p);
  if (i < left_exp_len)
    return x_left.access_symbol(x_left_p, i, g);
  else
    return x_right.access_symbol(x_right_p, i - left_exp_len, g);
}

#endif  // __SIMPLE_SLP_LCE_HPP_INCLUDED

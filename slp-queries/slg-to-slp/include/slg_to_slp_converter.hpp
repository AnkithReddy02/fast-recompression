/**
 * @file    slg_to_slp_converter.hpp
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

#ifndef __SLG_TO_SLP_CONVERTER_HPP_INCLUDED
#define __SLG_TO_SLP_CONVERTER_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>

#include "utils/space_efficient_vector.hpp"
#include "utils/packed_pair.hpp"
#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"


//=============================================================================
// Class used to represent the lazy AVL grammar. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct slg_to_slp_converter;

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
    typedef slg_to_slp_converter<char_type, text_offset_type> slg_to_slp_converter_type;
    typedef packed_pair<text_offset_type, text_offset_type> pair_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    std::uint8_t m_type;
    ptr_type m_left_p;
    ptr_type m_right_p;

  public:

    //=========================================================================
    // Constructors.
    //=========================================================================
    nonterminal();
    nonterminal(const char_type);
    nonterminal(const ptr_type, const ptr_type);

    //=========================================================================
    // Access methods.
    //=========================================================================
    char_type get_char() const;
    std::uint8_t get_type() const;
    ptr_type get_left_p() const;
    ptr_type get_right_p() const;

    //=========================================================================
    // Mostly unused.
    //=========================================================================
    char_type access_symbol(const ptr_type,
        std::uint64_t, const slg_to_slp_converter_type * const) const;
} __attribute__((packed));

//=============================================================================
// Implementation of the slg_to_slp_converter class.
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
struct slg_to_slp_converter {
  static_assert(sizeof(char_type) <= sizeof(text_offset_type),
      "Error: sizeof(char_type) > sizeof(text_offset_type)!");

  private:

    //=========================================================================
    // Declare typedefs.
    //=========================================================================
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;
    typedef text_offset_type ptr_type;
    typedef packed_pair<text_offset_type, text_offset_type> pair_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    nonterminal_type *m_nonterminals;
    std::uint64_t m_nonterminals_size;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    slg_to_slp_converter(const std::string slg_filename) {

      // Initialize the grammar reader.
      const std::uint64_t bufsize = (1 << 19);
      typedef async_stream_reader<text_offset_type> reader_type;
      reader_type *reader = new reader_type(slg_filename, bufsize, 4);

      // Read the number of roots.
      const std::uint64_t n_roots = reader->read();

      // Allocate array for roots, accouting for the sentinel root.
      ptr_type *roots_vec = utils::allocate_array<ptr_type>(n_roots);
      std::uint64_t roots_vec_size = 0;

      // Read roots.
      for (std::uint64_t i = 0; i < n_roots; ++i) {
        const std::uint64_t root_id = ((std::uint64_t)reader->read()) - 1;
        roots_vec[roots_vec_size++] = root_id;
      }

      // Compute the number of nonterminals.
      const std::uint64_t filesize = utils::file_size(slg_filename);
      const std::uint64_t skip_length =
        (n_roots + 1) * sizeof(text_offset_type);
      std::uint64_t n_nonterminals =
        (filesize - skip_length) / sizeof(pair_type);

      // Compute the number of nonterminals that will be added as a resuls
      // or merging the roots.
      std::uint64_t extra_nonterminals = 0;
      {
        std::uint64_t cur_roots_vec_size = roots_vec_size;
        while (cur_roots_vec_size > 1) {
          extra_nonterminals += cur_roots_vec_size / 2;
          cur_roots_vec_size = (cur_roots_vec_size + 1) / 2;
        }
      }

      // Allocate array for nonterminals.
      m_nonterminals = utils::allocate_array<nonterminal_type>(n_nonterminals + extra_nonterminals);
      m_nonterminals_size = 0;
      
      // Read nonterminals and compute number of long nonterminals.
      for (std::uint64_t i = 0; i < n_nonterminals; ++i) {
        const std::uint64_t val1 = reader->read();
        const std::uint64_t val2 = reader->read();
        if (val1 == 0) {
          nonterminal_type nonterm((char_type)val2);
          add_nonterminal(nonterm);
        } else {
          const std::uint64_t left_p = val1 - 1;
          const std::uint64_t right_p = val2 - 1;
          add_nonterminal(left_p, right_p);
        }
      }

      // Clean up.
      delete reader;

      // Add additional nonterminals so that SLG only has a single root.
      while (roots_vec_size > 1) {

        // compute the number of roots in the new round.
        std::uint64_t new_roots_vec_size = (roots_vec_size + 1) / 2;
        ptr_type *new_roots_vec = utils::allocate_array<ptr_type>(new_roots_vec_size);

        // Merge roots pairwise.
        for (std::uint64_t i = 0; i < roots_vec_size / 2; ++i) {
          const std::uint64_t left_p = roots_vec[2*i];
          const std::uint64_t right_p = roots_vec[2*i+1];
          const ptr_type nonterm_p = add_nonterminal(left_p, right_p);
          ++n_nonterminals;
          new_roots_vec[i] = nonterm_p;
        }

        // Copy the last root if roots_vec_size was odd.
        if (roots_vec_size & 1)
          new_roots_vec[new_roots_vec_size - 1] = roots_vec[roots_vec_size - 1];

        // Replace old roots with new ones.
        utils::deallocate(roots_vec);
        roots_vec = new_roots_vec;
        roots_vec_size = new_roots_vec_size;
      }
    }

    void write_slp_to_file(const std::string slp_filename) const {

      // Create the output writer.
      const std::uint64_t bufsize = (1 << 19);
      typedef async_stream_writer<text_offset_type> writer_type;
      writer_type *writer = new writer_type(slp_filename, bufsize, 4);

      // Write expansions of non-root nonterminals.
      // Note: some may be unused in the grammar!
      for (std::uint64_t i = 0; i < m_nonterminals_size; ++i) {
        const nonterminal_type &nonterm = get_nonterminal(i);
        const std::uint8_t type = nonterm.get_type();
        if (type == 0) {
          const char_type c = nonterm.get_char();
          writer->write((text_offset_type)0);
          writer->write((text_offset_type)c);
        } else {
          const std::uint64_t left = (std::uint64_t)nonterm.get_left_p() + 1;
          const std::uint64_t right = (std::uint64_t)nonterm.get_right_p() + 1;
          writer->write((text_offset_type)left);
          writer->write((text_offset_type)right);
        }
      }

      // Clean up.
      delete writer;
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~slg_to_slp_converter() {
      utils::deallocate(m_nonterminals);
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

  public:

    //=========================================================================
    // Gives access to a given nonterminal.
    //=========================================================================
    const nonterminal_type& get_nonterminal(const ptr_type id) const {
      return m_nonterminals[(std::uint64_t)id];
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

      // Create and add new nonterminal.
      const ptr_type new_nonterm_p = number_of_nonterminals();
      nonterminal_type new_nonterm(left_p, right_p);
      m_nonterminals[m_nonterminals_size++] = new_nonterm;

      // Return the ptr to the new nonterminal.
      return new_nonterm_p;
    }
};

//=============================================================================
// Default constructor.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal()
  : m_type(0),
    m_left_p(std::numeric_limits<text_offset_type>::max()),
    m_right_p(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for a nonterminal expanding to a single symbol.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(const char_type c)
  : m_type(0),
    m_left_p((text_offset_type)((std::uint64_t)c)),
    m_right_p(std::numeric_limits<text_offset_type>::max()) {}

//=============================================================================
// Constructor for non-single-symbol nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
nonterminal<char_type, text_offset_type>::nonterminal(
      const ptr_type left_p,
      const ptr_type right_p)
  : m_type(1),
    m_left_p(left_p),
    m_right_p(right_p) {}

//=============================================================================
// Return the char stored in a nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
char_type nonterminal<char_type, text_offset_type>::get_char() const {
  return (char_type)((std::uint64_t)m_left_p);
}

//=============================================================================
// Return the type of a nonterminal.
//=============================================================================
template<typename char_type, typename text_offset_type>
std::uint8_t nonterminal<char_type, text_offset_type>::get_type() const {
  return m_type;
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

#endif  // __SLG_TO_SLP_CONVERTER_HPP_INCLUDED


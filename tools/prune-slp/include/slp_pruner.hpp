/**
 * @file    slp_pruner.hpp
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

#ifndef __SLP_PRUNER_HPP_INCLUDED
#define __SLP_PRUNER_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>

#include "utils/space_efficient_vector.hpp"
#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"


//=============================================================================
// Class used to represent the lazy AVL grammar. Forward declaration.
//=============================================================================
template<typename char_type, typename text_offset_type>
struct slp_pruner;

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
    typedef slp_pruner<char_type, text_offset_type> slp_pruner_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    std::uint64_t m_name;
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
    std::uint64_t get_name() const;
    void set_name(const std::uint64_t);
    std::uint8_t get_type() const;
    ptr_type get_left_p() const;
    ptr_type get_right_p() const;

    //=========================================================================
    // Mostly unused.
    //=========================================================================
    char_type access_symbol(
        const ptr_type,
        std::uint64_t,
        const slp_pruner_type * const) const;
    std::uint64_t visit(
        const ptr_type,
        uint64_t *names,
        const slp_pruner_type * const,
        async_stream_writer<text_offset_type> *,
        std::uint64_t &) const;
} __attribute__((packed));

//=============================================================================
// Implementation of the slp_pruner class.
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
struct slp_pruner {
  static_assert(sizeof(char_type) <= sizeof(text_offset_type),
      "Error: sizeof(char_type) > sizeof(text_offset_type)!");

  private:

    //=========================================================================
    // Declare typedefs.
    //=========================================================================
    typedef nonterminal<char_type, text_offset_type> nonterminal_type;
    typedef text_offset_type ptr_type;

    //=========================================================================
    // Class members.
    //=========================================================================
    nonterminal_type *m_nonterminals;
    std::uint64_t m_nonterminals_size;

  public:

    //=========================================================================
    // Constructor.
    //=========================================================================
    slp_pruner(const std::string slp_filename) {

      // Initialize the grammar reader.
      const std::uint64_t bufsize = (1 << 19);
      typedef async_stream_reader<text_offset_type> reader_type;
      reader_type *reader = new reader_type(slp_filename, bufsize, 4);

      // Compute the number of nonterminals.
      const std::uint64_t filesize = utils::file_size(slp_filename);
      std::uint64_t n_nonterminals = filesize / (2 * sizeof(text_offset_type));

      // Allocate array for nonterminals.
      m_nonterminals = utils::allocate_array<nonterminal_type>(n_nonterminals);
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
    }

    void write_reachable_from_root(const std::string pruned_slp_filename) const {

      // Allocate array for "reduced" nonterminal identifiers.
      std::uint64_t *names = utils::allocate_array<std::uint64_t>(m_nonterminals_size);
      std::fill(names, names + m_nonterminals_size, std::numeric_limits<std::uint64_t>::max());

      // Initialize writer
      const std::uint64_t bufsize = (1 << 19);
      typedef async_stream_writer<text_offset_type> writer_type;
      writer_type *writer = new writer_type(pruned_slp_filename, bufsize, 4);

      // Write pruned grammar to file
      std::uint64_t visit_counter = 0;
      ptr_type root_id = get_root_id();
      const nonterminal_type &nonterm = get_nonterminal(root_id);
      nonterm.visit(root_id, names, this, writer, visit_counter);

      // Print stats.
      fprintf(stderr, "Size of original SLP = %lu\n", m_nonterminals_size);
      fprintf(stderr, "Size of pruned SLP = %lu (%.2Lf%%)\n",
          visit_counter, (100.L * visit_counter) / m_nonterminals_size);

      // Clean up.
      delete writer;
      utils::deallocate(names);
    }

    //=========================================================================
    // Destructor.
    //=========================================================================
    ~slp_pruner() {
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

template<typename char_type, typename text_offset_type>
std::uint64_t nonterminal<char_type, text_offset_type>::visit(
    ptr_type x_p,
    uint64_t *names,
    const slp_pruner_type * const g,
    async_stream_writer<text_offset_type> *writer,
    std::uint64_t &visited_count) const {
  if (names[(std::uint64_t)x_p] != std::numeric_limits<std::uint64_t>::max())
    return names[(std::uint64_t)x_p];
  const nonterminal_type &x = g->get_nonterminal(x_p);
  const std::uint8_t type = x.get_type();
  if (type == 0) {
    const char_type c = x.get_char();
    writer->write((text_offset_type)0);
    writer->write((text_offset_type)c);
  } else {
    const ptr_type x_left_p = x.get_left_p();
    const ptr_type x_right_p = x.get_right_p();
    const nonterminal_type &x_left = g->get_nonterminal(x_left_p);
    const nonterminal_type &x_right = g->get_nonterminal(x_right_p);
    const std::uint64_t left_name = x_left.visit(x_left_p, names, g, writer, visited_count);
    const std::uint64_t right_name = x_right.visit(x_right_p, names, g, writer, visited_count);
    writer->write((text_offset_type)left_name);
    writer->write((text_offset_type)right_name);
  }
  const std::uint64_t name = ++visited_count;
  names[(std::uint64_t)x_p] = name;
  return name;
}

#endif  // __SLP_PRUNER_HPP_INCLUDED


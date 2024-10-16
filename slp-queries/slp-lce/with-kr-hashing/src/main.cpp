/**
 * @file    main.cpp
 * @section LICENCE
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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <string>
#include <ctime>
#include <unistd.h>
#include <getopt.h>

#include "../include/types/uint40.hpp"
#include "../include/utils/utils.hpp"
#include "../include/simple_slp.hpp"
#include "../include/lce_naive.hpp"


//=============================================================================
// Read the text and compare to the text encoded by the grammar.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
void test_queries(
    const std::string text_filename,
    const std::string slp_filename) {

  // Get get length.
  std::uint64_t text_length =
    utils::file_size(text_filename) / sizeof(char_type);

  // Print initial messages.
  fprintf(stderr, "Test queries on SLP\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Text filename = %s\n", text_filename.c_str());
  fprintf(stderr, "SLP filename = %s\n", slp_filename.c_str());
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n",
      text_length, (1.L * text_length * sizeof(char_type) / (1 << 20)));
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Read grammar from file.
  typedef simple_slp<char_type, text_offset_type> slp_type;
  slp_type *slp = nullptr;
  {
    fprintf(stderr, "Read grammar from file... ");
    long double start = utils::wclock();
    slp = new slp_type(slp_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "DONE (%.2Lfs)\n", elapsed);
  }

  fprintf(stderr, "\n");

  // Obtain statistics.
  const std::uint64_t n_nonterminals = slp->number_of_nonterminals();
  const std::uint64_t total_rhs_length = slp->total_rhs_length();
  const std::uint64_t grammar_size = total_rhs_length;
  const std::uint64_t ram_use = slp->ram_use();

  // Print info. Note that the grammar may
  // still contain unused nonterminals.
  fprintf(stderr, "Statistics:\n");
  fprintf(stderr, "  Number of nonterminals = %lu\n", n_nonterminals);
  fprintf(stderr, "  Grammar size = %lu\n", grammar_size);
  fprintf(stderr, "  Grammar RAM use = %.2LfMiB\n",
      (1.L * ram_use) / (1UL << 20));
  fprintf(stderr, "\n");

  // Print RAM breakdown.
  slp->print_stats();
  fprintf(stderr, "\n");

  // Compare grammar output to text.
  {

    // Start the timer.
    fprintf(stderr, "Compare grammar output to text... ");
    long double start = utils::wclock();

    // Run the comparison.
    bool eq = slp->compare_expansion_to_text(text_filename, text_length);

    // Print summary.
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lf ", elapsed);
    fprintf(stderr, "(%s)\n", eq ? "OK" : "FAILED");
  }

  // Test correctness of random access queries.
  {

    // Print info.
    fprintf(stderr, "Test correctness of random access queries... ");

    // Prepare queries.
    static const std::uint64_t n_queries = (1 << 20);
    std::uint64_t *queries = utils::allocate_array<std::uint64_t>(n_queries);
    for (std::uint64_t i = 0; i < n_queries; ++i)
      queries[i] = utils::random_int<std::uint64_t>(0, text_length - 1);

    char_type *text = utils::allocate_array<char_type>(text_length);
    utils::read_from_file(text, text_length, text_filename);
    
    // Start the queries.
    for (std::uint64_t i = 0; i < n_queries; ++i) {
      std::uint64_t x = queries[i];
      std::uint64_t correct_access = text[x];
      std::uint64_t computed_access = slp->access_symbol(x);
      if (correct_access != computed_access) {
        fprintf(stderr, "\nError:\n");
        fprintf(stderr, "\tx = %lu\n", x);
        fprintf(stderr, "\tcorrect RA = %lu\n", correct_access);
        fprintf(stderr, "\tcomputed RA = %lu\n", computed_access);
        std::exit(EXIT_FAILURE);
      }
    }

    // Clean up.
    utils::deallocate(text);
    utils::deallocate(queries);
    fprintf(stderr, "OK\n");
  }

  // Test correctness of LCE queries.
  {

    // Print info.
    fprintf(stderr, "Test correctness of LCE queries... ");

    // Prepare queries.
    static const std::uint64_t n_queries = (1 << 23);
    typedef std::pair<std::uint64_t, std::uint64_t> pair_type;
    pair_type *queries = utils::allocate_array<pair_type>(n_queries);
    for (std::uint64_t i = 0; i < n_queries; ++i) {
      const std::uint64_t pos_i = utils::random_int<std::uint64_t>(0, text_length - 1);
      const std::uint64_t pos_j = utils::random_int<std::uint64_t>(0, text_length - 1);
      queries[i] = std::make_pair(pos_i, pos_j);
    }

    typedef lce_naive<char_type, text_offset_type> lce_naive_type;
    lce_naive_type *naive = new lce_naive_type(text_filename);
    
    // Start the queries.
    std::uint64_t lce_sum = 0;
    std::uint64_t lce_max = 0;
    for (std::uint64_t i = 0; i < n_queries; ++i) {
      std::uint64_t x = queries[i].first;
      std::uint64_t y = queries[i].second;
      std::uint64_t correct_lce = naive->lce(x, y);
      std::uint64_t computed_lce = slp->lce(x, y);
      if (correct_lce != computed_lce) {
        fprintf(stderr, "\nError:\n");
        fprintf(stderr, "\tx = %lu\n", x);
        fprintf(stderr, "\ty = %lu\n", y);
        fprintf(stderr, "\tcorrect LCE = %lu\n", correct_lce);
        fprintf(stderr, "\tcomputed LCE = %lu\n", computed_lce);
        std::exit(EXIT_FAILURE);
      }

      lce_sum += computed_lce;
      lce_max = std::max(lce_max, computed_lce);
    }

    // Clean up.
    delete naive;
    utils::deallocate(queries);
    fprintf(stderr, "OK\n");
    fprintf(stderr, "  Max LCE = %lu\n", lce_max);
    fprintf(stderr, "  Avg LCE = %.2Lf\n",
        (long double)lce_sum / (long double)n_queries);
  }

  // Measure time for random access.
  {

    // Print info.
    fprintf(stderr, "Benchmark random access... ");

    // Prepare queries.
#ifndef NDEBUG
    static const std::uint64_t n_queries = (1 << 16);
#else
    static const std::uint64_t n_queries = (1 << 20);
#endif
    std::uint64_t *queries = utils::allocate_array<std::uint64_t>(n_queries);
    for (std::uint64_t i = 0; i < n_queries; ++i)
      queries[i] = utils::random_int<std::uint64_t>(0, text_length - 1);
    const static std::uint64_t reps = 10;

    // Initialize stats.
    std::uint64_t checksum = 0;
    std::uint64_t ptr = 0;

    // Start the queries.
    long double start = utils::wclock();
    for (std::uint64_t i = 0; i < n_queries * reps; ++i) {
      checksum += slp->access_symbol(queries[ptr++]);
      if (ptr == n_queries)
        ptr = 0;
    }
    if (checksum == 849302)
      std::exit(EXIT_FAILURE);

    // Measure the time.
    long double elapsed = utils::wclock() - start;
    elapsed *= 1000000.0L;
    elapsed /= (n_queries * reps);
    fprintf(stderr, "%.2Lfus/query\n", elapsed);

    // Clean up.
    utils::deallocate(queries);
  }

  // Measure time for LCE queries.
  {

    // Print info.
    fprintf(stderr, "Benchmark LCE queries... ");

    // Prepare queries.
#ifndef NDEBUG
    static const std::uint64_t n_queries = (1 << 16);
#else
    static const std::uint64_t n_queries = (1 << 20);
#endif
    typedef std::pair<std::uint64_t, std::uint64_t> pair_type;
    pair_type *queries = utils::allocate_array<pair_type>(n_queries);
    for (std::uint64_t i = 0; i < n_queries; ++i) {
      const std::uint64_t pos_i =
        utils::random_int<std::uint64_t>(0, text_length - 1);
      const std::uint64_t pos_j =
        utils::random_int<std::uint64_t>(0, text_length - 1);
      queries[i] = std::make_pair(pos_i, pos_j);
    }
    const static std::uint64_t reps = 10;

    // Initialize stats.
    std::uint64_t checksum = 0;
    std::uint64_t ptr = 0;

    // Start the queries.
    long double start = utils::wclock();
    for (std::uint64_t i = 0; i < n_queries * reps; ++i) {
      checksum += slp->lce(queries[ptr].first, queries[ptr].second);
      ++ptr;
      if (ptr == n_queries)
        ptr = 0;
    }
    if (checksum == 849302)
      std::exit(EXIT_FAILURE);

    // Measure the time.
    long double elapsed = utils::wclock() - start;
    elapsed *= 1000000.0L;
    elapsed /= (n_queries * reps);
    fprintf(stderr, "%.2Lfus/query\n", elapsed);

    // Clean up.
    utils::deallocate(queries);
  }

  // Clean up.
  delete slp;

  // Print summary.
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
}

int main(int argc, char **argv) {
  if (argc != 3)
    std::exit(EXIT_FAILURE);

  // Init rand.
  srand(time(0) + getpid());

  // Initialize runtime statistics.
  utils::initialize_stats();

  // Declare types.
  typedef std::uint8_t char_type;
  typedef uint40 text_offset_type;

  // Obtain filenames.
  std::string text_filename = argv[1];
  std::string slp_filename = argv[2];

  // Turn paths absolute.
  text_filename = utils::absolute_path(text_filename);
  slp_filename = utils::absolute_path(slp_filename);

  // Run the algorithm.
  test_queries<char_type, text_offset_type>(
      text_filename, slp_filename);
}


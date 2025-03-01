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
#include "../include/slp_pruner.hpp"
#include "../include/typedefs.hpp"

//=============================================================================
// Read the text and compare to the text encoded by the grammar.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
void prune_slp(
    const std::string slp_filename,
    const std::string pruned_slp_filename) {

  // Print initial messages.
  fprintf(stderr, "Prune SLP\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "SLP filename = %s\n", slp_filename.c_str());
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Read grammar from file.
  typedef slp_pruner<char_type, text_offset_type> slp_pruner_type;
  slp_pruner_type *slp_pruner = nullptr;
  {
    fprintf(stderr, "Read SLP from file... ");
    long double start = utils::wclock();
    slp_pruner = new slp_pruner_type(slp_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "DONE (%.2Lfs)\n", elapsed);
  }

  // Prune SLP.
  slp_pruner->write_reachable_from_root(pruned_slp_filename);

  // Clean up.
  delete slp_pruner;

  // Print summary.
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
}

int main(int argc, char **argv) {
  if (argc != 3)
    std::exit(EXIT_FAILURE);

  // Initialize runtime statistics.
  utils::initialize_stats();

  // Declare types.
  typedef std::uint8_t char_type;
  typedef c_size_t text_offset_type;

  // Obtain filenames.
  std::string slp_filename = argv[1];
  std::string pruned_slp_filename = argv[2];

  // Turn paths absolute.
  slp_filename = utils::absolute_path(slp_filename);
  pruned_slp_filename = utils::absolute_path(pruned_slp_filename);

  // Run the algorithm.
  prune_slp<char_type, text_offset_type>(
      slp_filename, pruned_slp_filename);
}


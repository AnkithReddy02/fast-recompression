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
#include "../include/slg_to_slp_converter.hpp"
#include "../include/typedefs.hpp"

//=============================================================================
// Convert SLG to SLP.
//=============================================================================
template<
  typename char_type,
  typename text_offset_type>
void convert_slg_to_slp(
    const std::string slg_filename,
    const std::string slp_filename) {

  // Print initial messages.
  fprintf(stderr, "Convert SLG to SLP\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "SLG filename = %s\n", slg_filename.c_str());
  fprintf(stderr, "SLP filename = %s\n", slp_filename.c_str());
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Read SLG from file.
  typedef slg_to_slp_converter<char_type, text_offset_type> converter_type;
  converter_type *converter = nullptr;
  {
    fprintf(stderr, "Read SLG from file and convert to SLP... ");
    long double start = utils::wclock();
    converter = new converter_type(slg_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "DONE (%.2Lfs)\n", elapsed);
  }

  // Obtain statistics.
  const std::uint64_t n_nonterminals = converter->number_of_nonterminals();

  // Print info. Note that the grammar may
  // still contain unused nonterminals.
  fprintf(stderr, "Statistics:\n");
  fprintf(stderr, "  Number of nonterminals = %lu\n", n_nonterminals);
  fprintf(stderr, "\n");

  // Write SLP to file
  {
    fprintf(stderr, "Write SLP to file... ");
    long double start = utils::wclock();
    converter->write_slp_to_file(slp_filename);
    long double elapsed = utils::wclock() - start;
    fprintf(stderr, "DONE (%.2Lfs)\n", elapsed);
  }

  // Clean up.
  delete converter;

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
  std::string slg_filename = argv[1];
  std::string slp_filename = argv[2];

  // Turn paths absolute.
  slg_filename = utils::absolute_path(slg_filename);
  slp_filename = utils::absolute_path(slp_filename);

  // Run the conversion
  convert_slg_to_slp<char_type, text_offset_type>(slg_filename, slp_filename);
}


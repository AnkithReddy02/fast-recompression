/**
 * @file    main.cpp
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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include <ctime>
#include <unistd.h>
#include <getopt.h>

#include "../include/types/uint40.hpp"
#include "../include/utils/utils.hpp"
#include "../include/io/async_stream_reader.hpp"


//=============================================================================
// Extend LZ parsing by a given number of phrases.
//=============================================================================
template<
  typename char_type = std::uint8_t,
  typename text_offset_type = std::uint64_t>
void extend_lz(
    std::string input_lz_filename,
    std::string output_lz_filename,
    const std::uint64_t max_phrase_length,
    const std::uint64_t n_new_phrases) {

  // Declare typedefs.
  typedef std::pair<text_offset_type, text_offset_type> pair_type;

  // Turn paths absolute.                                                       
  input_lz_filename = utils::absolute_path(input_lz_filename);
  output_lz_filename = utils::absolute_path(output_lz_filename);

  // Print parameters.
  fprintf(stderr, "Extend LZ77 parsing\n");
  fprintf(stderr, "Timestamp = %s", utils::get_timestamp().c_str());
  fprintf(stderr, "Input LZ filename = %s\n", input_lz_filename.c_str());
  fprintf(stderr, "Output LZ filename = %s\n", output_lz_filename.c_str());
  fprintf(stderr, "sizeof(char_type) = %lu\n", sizeof(char_type));
  fprintf(stderr, "sizeof(text_offset_type) = %lu\n", sizeof(text_offset_type));
  fprintf(stderr, "\n\n");

  // Compute number of phrases.
  const std::uint64_t n_phrases =
    utils::file_size(input_lz_filename) / (2 * sizeof(text_offset_type));

  // Read the parsing.
  std::vector<pair_type> parsing;
  parsing.reserve(n_phrases);
  {
    // Initialize parsing reader.
    typedef async_stream_reader<pair_type> reader_type;
    reader_type *reader = new reader_type(input_lz_filename);

    // Stream parsing.
    for (std::uint64_t phrase_id = 0; phrase_id < n_phrases; ++phrase_id)
      parsing.push_back(reader->read());

    // Clean up.
    delete reader;
  }

  // Compute text length.
  std::uint64_t text_length = 0;
  for (std::uint64_t phrase_id = 0; phrase_id < n_phrases; ++phrase_id) {
    const std::uint64_t len = parsing[phrase_id].second;
    text_length += std::max(len, (std::uint64_t)1);
  }

  // Print text length and number of phrases.
  fprintf(stderr, "Parsing size = %lu (%.2LfMiB)\n",
      n_phrases,
      (long double)(n_phrases * sizeof(pair_type)) / (1 << 20));
  fprintf(stderr, "Text length = %lu (%.2LfMiB)\n",
      text_length,
      (long double)text_length / (1 << 20));

  // Extend the parsing.
  std::uint64_t cur_text_length = text_length;
  for (std::uint64_t i = 0; i < n_new_phrases; ++i) {
    std::uint64_t pos = utils::random_int<std::uint64_t>(0, cur_text_length - 1);
    std::uint64_t len = utils::random_int<std::uint64_t>(1, max_phrase_length);
    parsing.push_back(
        std::make_pair(
          (text_offset_type)pos,
          (text_offset_type)len));
    cur_text_length += len;
  }

  // Print new text length and new number of phrases.
  fprintf(stderr, "New parsing size = %lu (%.2LfMiB)\n",
      (std::uint64_t)parsing.size(),
      (long double)(parsing.size() * sizeof(pair_type)) / (1 << 20));
  fprintf(stderr, "New text length = %lu (%.2LfMiB)\n",
      cur_text_length,
      (long double)cur_text_length / (1 << 20));

  // Write extended parsing to file.
  {
    fprintf(stderr, "Write parsing to file... ");
    const long double start = utils::wclock();
    const pair_type * const parsing_data = parsing.data();
    utils::write_to_file<pair_type>(parsing_data,
        parsing.size(), output_lz_filename);
    const long double elapsed = utils::wclock() - start;
    fprintf(stderr, "%.2Lfs\n", elapsed);
  }

  // Print summary.
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  RAM allocation: cur = %lu bytes, peak = %.2LfMiB\n",
      utils::get_current_ram_allocation(),
      (1.L * utils::get_peak_ram_allocation()) / (1UL << 20));
}

int main(int argc, char **argv) {
  if (argc != 5)
    std::exit(EXIT_FAILURE);

  // Initialize runtime statistics.
  utils::initialize_stats();

  // Declare types.
  typedef std::uint8_t char_type;
  typedef uint40 text_offset_type;

  // Obtain filenames.
  std::string input_lz_filename = argv[1];
  std::uint64_t max_phrase_length = std::atoi(argv[2]);
  std::uint64_t n_new_phrases = std::atoi(argv[3]);
  std::string output_lz_filename = argv[4];

  // Run the algorithm.
  extend_lz<char_type, text_offset_type>(
      input_lz_filename,
      output_lz_filename,
      max_phrase_length,
      n_new_phrases);
}


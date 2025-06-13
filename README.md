# Fast-Recompression

## Description

This is a C++ implementation of a fast and space-efficient grammar recompression algorithm that converts a straight-line program (SLP) into a run-length straight-line program (RLSLP). This implementation accompanies the manuscript in [1]. Recompression is a powerful technique for constructing locally consistent grammars introduced in [2]. In our implementation, we follow closely the details from [3].

The system includes a full compression pipeline that builds an RLSLP from plain text using approximate LZ77 parsing (via Bentley–McIlroy factorization), construction of grammars (SLG), SLP generation, and pruning.

The repository also includes a validation tool that checks correctness of the recompressed grammar using randomized LCE queries compared against a trusted SLP implementation.

---

## Compilation and Usage

Fast-Recompression is compiled by simply typing `make` in the directory containing this README. This builds the executable:

- `recomp` — Builds an RLSLP from SLP

To build the binary:

```bash
$ make
```

---

### Example: Recompressing an SLP

```bash
$ ./recomp [-v or -vv] input.slp -o output.rlslp [-t input.txt]
```

This reads a straight-line program from `input.slp` and recompresses it into a run-length SLP written to `output.rlslp`. The optional `-t input.txt` flag can be provided to validate that the recompressed grammar expands to exactly the original string `input.txt`, which corresponds to the input that was originally compressed into `input.slp`.

Use `-v` or `-vv` for verbose output.

---

### Full Pipeline: Text -> Approx. LZ77 -> SLG -> SLP -> PRUNED_SLP -> RLSLP

A script is provided to run the full compression pipeline starting from a plain text file:

```bash
$ ./text_rlslp_script.sh input.txt [block_size]
```

- `input.txt` — required text file  
- `block_size` — optional positive integer parameter used for Bentley–McIlroy factorization LZ77 factorization (default: 50)

All output files generated at each stage are written to the same directory as the input, using the prefix `input.txt.` (e.g., `input.txt.rlslp`, etc.).

All intermediate steps are logged in `text_rlslp_script.log`.

To run individual steps manually, use the corresponding binaries in the `tools/` directory.

---

### RLSLP Correctness Verification via Random LCE Queries

To build the test binary:

```bash
$ cd test
$ make
```

Once built, you can verify the correctness of a recompressed RLSLP by comparing it against the original SLP using randomly generated Longest Common Extension (LCE) queries:

```bash
$ ./test input.rlslp input.slp
```

Here, `input.slp` is the grammar input to `recomp`, and `input.rlslp` is its recompressed output.

- Computes `LCE(i, j)` on both grammars for random positions `i` and `j`,
- Uses a trusted implementation [`test/include/simple_slp.hpp`](test/include/simple_slp.hpp), based on Karp–Rabin fingerprinting, to evaluate the result on the original SLP,
- Compares the result against the corresponding RLSLP implementation defined in [`src/recompression_definitions.cpp`](src/recompression_definitions.cpp),
- Stops and reports the first mismatch (if any).

---

## Dependencies

- C++17 compiler  
- GNU Make  
- Standard UNIX environment (`bash`, `awk`, etc.)

---

## References

[1] Ankith Reddy Adudodla, Dominik Kempa. *Engineering Fast and Space-Efficient Recompression from SLP-Compressed Text*. 2025.<br>
[2] Artur Jez. *Recompression: A Simple and Powerful Technique for Word Equations*. J. ACM 63(1): 4:1–4:51 (2016).<br>
[3] Dominik Kempa, Tomasz Kociumaka. *Collapsing the Hierarchy of Compressed Data Structures: Suffix Arrays in Optimal Compressed Space*. FOCS 2023: 1877–1886.

---

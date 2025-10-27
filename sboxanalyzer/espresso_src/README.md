# Espresso Logic Minimizer - Modernized Edition

A modernized version of the Espresso heuristic logic minimizer with ANSI C conformance and enhanced reliability.

## About

This repository contains a maintained and modernized version of the Espresso logic minimizer, originally developed at the University of California, Berkeley. Espresso is a heuristic multi-valued PLA (Programmable Logic Array) minimizer that has been widely used in electronic design automation and logic synthesis.

### Original Source

The original Espresso implementation was developed by the University of California, Berkeley and is available from:
- Official Distribution: https://ptolemy.berkeley.edu/projects/embedded/pubs/downloads/espresso/index.htm
- Copyright (c) 1988, 1989, Regents of the University of California

The original software was distributed under a permissive BSD-style license allowing use, copying, and preparation of derivative works.

## What's New in This Version

This modernized version addresses several improvements while maintaining 100% functional compatibility with the original implementation:

### Code Modernization

- **Cross-Platform Support**: Builds cleanly on Linux, macOS, and Windows (MSVC)

- **ANSI C Compliance**: Converted all Kernighan & Ritchie (K&R) C function declarations to ANSI C standard
  - Original: 2,173 compiler warnings on modern compilers (AppleClang 17.0.0)
  - Modernized: 0 warnings with strict compilation flags
  
- **Enhanced Type Safety**: Modern function prototypes with explicit parameter types

- **Improved Maintainability**: Clearer function signatures for easier understanding and maintenance

### Verification

All changes have been rigorously verified to ensure zero functional changes:
- Tested against 183 example files across three test suites
- Binary-identical output verified using SHA-256 hashing
- All modifications preserve original algorithms and logic flow

## Building

### Prerequisites

- CMake 3.10 or higher
- C compiler (GCC, Clang, MSVC, or compatible)
- Make or Ninja build system (Windows: Visual Studio or MSBuild)

### Build Instructions

**Unix/Linux/macOS:**

```bash
# Configure and build
cmake -B build
cmake --build build

# Optional: Install system-wide (requires permissions)
sudo cmake --install build

# The espresso binary will be in build/espresso
```

**Windows (Visual Studio):**

```bash
# Configure with Visual Studio generator
cmake -G "Visual Studio 17 2022" -A x64 -B build

# Build
cmake --build build --config Release

# The espresso.exe will be in build\Release\
```

## Usage

After building, you can run Espresso on PLA files:

```bash
cd build

# Basic minimization (outputs to stdout)
./espresso ../examples/examples/b2

# Save output to a file
./espresso ../examples/examples/alu2 > minimized.pla

# Get help and see all available options
./espresso

# Use exact minimization algorithm for optimal results
./espresso -Dexact ../examples/examples/b3

# Enable execution trace for detailed progress information
./espresso -t input.pla

# Show execution summary with timing and statistics
./espresso -s input.pla

# Advanced: Exact algorithm with strong minimization and fd output format
./espresso -Dexact -estrong -ofd input.pla
```

The output is a minimized PLA representation in standard format. For complete command-line options, run `./espresso` without arguments. For PLA file format specification, see `man/espresso.5`.

### Command-Line Options

**Main Algorithm Modes (-D):**
- `-Dexact` - Exact minimization algorithm (guarantees minimum number of product terms, and heuristically minimizes number of literals). Potentially expensive.
- `-Dmany` - Reads and minimizes PLA's until end-of-file is detected. PLA's in the same file are separated by `.e`
- `-Dsimplify` - Quick simplification of the cover (complementation followed by single containment)
- `-Dso` - Minimize each function one at a time as a single-output function. Terms will not be shared among the functions
- `-Dso_both` - Minimize each function one at a time as a single-output function, but choose the function or its complement based on which has fewer terms
- `-Dopo` - Perform output phase optimization (i.e., determine which functions to complement to reduce the number of terms needed to implement the function)
- `-Dopoall` - Minimize the function with all possible phase assignments (exponential number of minimizations)

**Verification and Analysis (-D):**
- `-Dcheck` - Checks that the function is a partition of the entire space (i.e., that the ON-set, OFF-set and DC-set are pairwise disjoint, and that their union is the Universe)
- `-Dstats` - Provide simple statistics on the size of the function
- `-Dverify` - Checks for Boolean equivalence of two PLA's. Reads two filenames from the command line, each containing a single PLA

**Espresso Options (-e):**
- `-efast` - Stop after the first EXPAND and IRREDUNDANT operations (i.e., do not iterate over the solution)
- `-estrong` - Uses the alternate strategy SUPER_GASP (as a replacement for LAST_GASP) which is more expensive, but occasionally provides better results
- `-eeat` - Normally comments are echoed from the input file to the output file. This option discards any comments in the input file
- `-enirr` - The result will not necessarily be made irredundant in the final step which removes redundant literals
- `-eness` - Essential primes will not be detected
- `-epos` - Swaps the ON-set and OFF-set of the function after reading the function. This can be used to minimize the OFF-set of a function
- `-eonset` - Recompute the ON-set before the minimization. Useful when the PLA has a large number of product terms (e.g., an exhaustive list of minterms)
- `-enunwrap` - The ON-set will not be unwrapped before beginning the minimization

**Output Options (-o):**
- `-of`, `-ofd`, `-ofr`, `-ofdr` - Select any combination of the ON-set (f), the OFF-set (r), or the DC-set (d)
- `-oeqntott` - Output algebraic equations
- `-opleasure` - Output an unmerged PLA format

**Other Options:**
- `-s` - Print a short summary of the execution including initial cost, final cost, and computer resources used
- `-t` - Produce a trace showing the execution of the program. After each main step, a line is printed reporting processor time used and current cost
- `-x` - Suppress printing of the solution
- `-d` - Enable debugging (useful only for those familiar with the algorithms used)
- `-v[type]` - Verbose debugging detail (not generally useful)
- `-Sn` - Select strategy for subcommands (opo, opoall, pair, pairall, so, so_both)
- `-rn-m` - Select range for subcommands (first and last outputs or variables)

### Example Files

The repository includes comprehensive example files in three categories:

- `examples/examples/` - Standard test cases (123 files)
- `examples/hard_examples/` - Computationally intensive cases (19 files)
- `examples/tlex/` - Additional test cases (41 files)

## Testing

A comprehensive test script is provided to verify the correctness of the implementation:

```bash
# Run all tests (183 example files with multiple modes)
./test.sh
```

### Test Suite Features

- **Multi-Mode Testing**: Each example is tested with 2-3 different Espresso modes (default, fast, strong, etc.)
- **Comprehensive Coverage**: Tests all 183 example files across 8 different algorithm modes
- **Smart Timeout Handling**: 59-second timeout per test to handle computationally intensive cases
- **Deterministic Execution**: Files processed in sorted order for consistent results across systems
- **SHA-256 Verification**: Cryptographically secure hashing validates output correctness
- **Mode-Specific Optimizations**: 
  - Hard examples use only fast modes (no exact minimization)
  - Specific problematic files excluded from certain modes to prevent timeouts
  - Special handling for large examples (e.g., o64.pla uses only `-Dsimplify`)

### Expected Test Results

```
Total files: 183
Total tests run: 369 (multiple modes per file)
Unique modes tested: 7
Successful tests: 369
Failed: 0
Timed out: 0

Final combined hash (sha256):
53f911764ba4d1799b25b43c20b23f08abe0df036fa8c76cccaf3854b8d7dd56

✓ All tests passed successfully!
  Hash matches expected value
```

The test suite verifies that the modernized code produces **identical output** to the original implementation, ensuring that all ANSI C conversions and code modernizations preserve the exact functionality of the original Espresso algorithm.

## Documentation

- Manual pages are available in `man/espresso.1` and `man/espresso.5`
- Original documentation is preserved in the `man/` directory
- Format specification: `man/espresso.5`

## Project Structure

```
espresso/
├── espresso/           # Core source files
├── utility/            # Utility functions and headers
├── examples/           # Test cases and examples
│   ├── examples/       # Standard examples
│   ├── hard_examples/  # Computationally intensive examples
│   └── tlex/          # Additional test cases
├── man/               # Manual pages
├── CMakeLists.txt     # Build configuration
└── test.sh            # Comprehensive test suite
```

## License

This modernized version is released under the GNU General Public License v3.0 (GPLv3).

The original Espresso software is copyright (c) 1988, 1989, Regents of the University of California, and was distributed under a permissive BSD-style license. The original copyright notice and license terms have been preserved in `utility/copyright.h`.

As required by the original license:
```
Copyright (c) 1988, 1989, Regents of the University of California.
All rights reserved.

Use and copying of this software and preparation of derivative works
based upon this software are permitted. However, any distribution of
this software or derivative works must include the above copyright notice.
```

The derivative work in this repository includes the required copyright notice and is distributed under GPLv3, which is compatible with the original BSD-style license.

## Contributing

Contributions are welcome. When contributing, please:

- Maintain ANSI C compliance
- Ensure all changes pass the test suite
- Preserve compatibility with the original Espresso behavior
- Follow the existing code style
- Include appropriate tests for new features

## Acknowledgments

- Original authors and contributors at UC Berkeley
- The OCT Tools Distribution 3.0 team
- CHIPS Alliance for their Espresso repository (https://github.com/chipsalliance/espresso)
- All contributors to the Espresso project over the decades

## Contact

For issues, questions, or contributions related to this modernized version, please use the GitHub issue tracker.

## References

- Original Espresso Distribution: https://ptolemy.berkeley.edu/projects/embedded/pubs/downloads/espresso/index.htm
- UC Berkeley Embedded Systems: https://ptolemy.berkeley.edu/projects/embedded/

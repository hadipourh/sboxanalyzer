# S-box Analyzer

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**S-box Analyzer**

S-box Analyzer is a tool for analyzing (vectorial) Boolean functions such as S-boxes against differential, linear, differential-linear, boomerang, and integral attacks.

Particularly, it derives the minimized CP/MILP and SMT/SAT constraints to encode a set of binary vectors, (vectorial) Boolean functions, the Differential Distribution Table (DDT), Linear Approximation Table (LAT), [Differential-Linear Connectivity Table (DLCT)](https://ia.cr/2024/255), Boomerang Connectivity Table (BCT), and [Monomial Prediction Table (MPT)](https://tosc.iacr.org/index.php/ToSC/article/view/9715) of S-boxes.

---

![logo](./images/sboxanalyzer.svg)

- [S-box Analyzer](#s-box-analyzer)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Convert a Set of Binary Vectors to CNF/MILP Constraints](#convert-a-set-of-binary-vectors-to-cnfmilp-constraints)
  - [Convert a Truth Table to CNF/MILP Constraints](#convert-a-truth-table-to-cnfmilp-constraints)
  - [DDT Encoding](#ddt-encoding)
    - [Encoding the DDT of SKINNY-64](#encoding-the-ddt-of-skinny-64)
    - [Encoding the DDT of Ascon](#encoding-the-ddt-of-ascon)
    - [Encoding the DDT of PRESENT](#encoding-the-ddt-of-present)
    - [Encoding the DDT of SKINNY-128](#encoding-the-ddt-of-skinny-128)
    - [Encoding the DDT of AES](#encoding-the-ddt-of-aes)
    - [Encoding the \*-DDT and \*-LAT](#encoding-the--ddt-and--lat)
    - [Encoding the DDT for CryptoSMT](#encoding-the-ddt-for-cryptosmt)
  - [Deriving and Encoding Deterministic Differential Propagation](#deriving-and-encoding-deterministic-differential-propagation)
  - [LAT Encoding](#lat-encoding)
    - [Encoding the LAT of SKINNY-64](#encoding-the-lat-of-skinny-64)
    - [Encoding the LAT of Ascon](#encoding-the-lat-of-ascon)
  - [Deriving and Encoding Deterministic Linear Propagation](#deriving-and-encoding-deterministic-linear-propagation)
  - [MPT Encoding](#mpt-encoding)
    - [Encoding the MPT of PRESENT](#encoding-the-mpt-of-present)
    - [Encoding the MPT of Ascon](#encoding-the-mpt-of-ascon)
  - [Convex Hull Constraints for MILP](#convex-hull-constraints-for-milp)
  - [DLCT Encoding](#dlct-encoding)
    - [Encoding the \*-DLCT of KNOT](#encoding-the--dlct-of-knot)
    - [Encoding the \*-DLCT of Midori](#encoding-the--dlct-of-midori)
    - [Verification of Hadipour et al.'s Theorem](#verification-of-hadipour-et-als-theorem)
  - [BCT Encoding](#bct-encoding)
  - [Exact Pair-Relation Encoding](#exact-pair-relation-encoding)
    - [Quasidifferential Trail Constraints](#quasidifferential-trail-constraints)
    - [Mix-Basis Trail Constraints](#mix-basis-trail-constraints)
    - [Valuation-Weighted Constraints](#valuation-weighted-constraints)
    - [Exact Boomerang Pair-Relation Encoding](#exact-boomerang-pair-relation-encoding)
  - [Citation](#citation)
  - [Papers](#papers)
  - [Road Map](#road-map)
  - [License](#license)

---

## Dependencies

S-box Analyzer requires Python 3.10 or higher and uses [PassageMath](https://github.com/passagemath/passagemath) (modular SageMath) for cryptographic primitives and ESPRESSO for logic minimization.

**All dependencies are automatically installed via pip:**

- **PassageMath** (≥10.6.30) - The distributions `passagemath-modules`, `passagemath-brial`, `passagemath-repl`, and `passagemath-symbolics`, which provide SageMath's S-box utilities and Boolean polynomial rings; PassageMath Standard contains them as well
- **ESPRESSO** - Pre-built binaries included for all major platforms (macOS, Linux, Windows)

**Optional (only if pre-built binary doesn't work):**
- **CMake** - For building ESPRESSO from source as fallback
- **C compiler** (gcc/clang/MSVC) - For building ESPRESSO from source

**Optional (only for the [convex hull constraints](#convex-hull-constraints-for-milp)):**
- **passagemath-polyhedra** - For `milp_method="convex_hull"`; it is part of PassageMath Standard

## Installation

Install directly from PyPI using pip:

```bash
pip install sboxanalyzer
```

To also install the optional dependency of the [convex hull constraints](#convex-hull-constraints-for-milp):

```bash
pip install "sboxanalyzer[hull]"
```

Or install from source:

```bash
git clone https://github.com/hadipourh/sboxanalyzer.git
cd sboxanalyzer
pip install -e .
```

**That's it!** The package includes pre-built ESPRESSO binaries for:
- macOS (Intel & Apple Silicon)
- Linux (x86_64 & ARM64)
- Windows (x86_64)

If the pre-built binary doesn't work on your system, ESPRESSO will be automatically built from source on first use (requires CMake and a C compiler).

## Usage

After installation, you can use S-box Analyzer directly in Python:

```python
from pprint import pprint
from sboxanalyzer import SboxAnalyzer
from sage.crypto.sboxes import GIFT as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_diff_constraints()
pprint(milp)
```

Example output:

```text
Number of constraints: 49
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
Weight: 3.0000 p0 + 2.0000 p1 + 1.4150 p2
['- p0 - p1 >= -1',
 '- p0 - p2 >= -1',
 '- p1 - p2 >= -1',
 'a1 - a3 + p0 >= 0',
 '- a0 + p0 + p1 >= 0',
 'a2 + b1 - p2 >= 0',
 'b0 + b2 + b3 - p0 >= 0',
 'a0 + a3 - b0 + p0 >= 0',
 '- a3 - b2 + b3 + p0 >= -1',
 '- a0 - a1 - b0 - p1 >= -3',
 '- a2 + a3 - b0 - p1 >= -2',
 'a1 + a2 + b2 - p1 >= 0',
 'a1 + a2 + b3 - p1 >= 0',
 '- a0 - a1 + a2 + a3 - b3 >= -2',
 'a1 - a2 + a3 - b2 - b3 >= -2',
 'a0 + a1 + b0 - b2 - b3 >= -1',
 '- a1 - a3 + b0 - b2 - b3 >= -3',
 '- a0 + b0 - b1 + b2 - b3 >= -2',
 'a0 + b0 + b1 + b2 - b3 >= 0',
 '- a1 + b0 - b1 - b2 + b3 >= -2',
 'a1 + b0 + b1 - b2 + b3 >= 0',
 'a1 - a2 + a3 + b2 + b3 >= 0',
 'a0 + a3 + b2 + b3 - p0 >= 0',
 'a0 - a1 - a2 - b1 + p0 >= -2',
 'a0 + a2 + b1 - b2 + p0 >= 0',
 '- a1 - a3 + b2 - b3 + p0 >= -2',
 'a0 - b1 + b2 + b3 + p0 >= 0',
 '- a1 + b1 + b2 + b3 + p0 >= 0',
 'a3 + b0 + b2 - b3 + p1 >= 0',
 'a3 + b0 - b2 + b3 + p1 >= 0',
 'a0 - a1 - a2 - b0 - b1 - b3 >= -4',
 '- a0 - a2 - a3 - b1 + b2 - b3 >= -4',
 '- a0 + a2 - b0 + b1 + b2 - b3 >= -2',
 'a0 - a1 - b0 - b2 - b3 + p1 >= -3',
 '- a0 + a1 - b0 - b2 - b3 + p1 >= -3',
 '- a0 - a1 - a3 + b2 + b3 + p1 >= -2',
 'a0 + a2 + a3 - b1 - b2 + p2 >= -1',
 'a1 + a2 + a3 - p0 + p1 + p2 >= 0',
 '- a1 + a2 - a3 - b0 + b1 + b2 - b3 >= -3',
 '- a0 - a2 - a3 - b0 + b1 - b2 + b3 >= -4',
 '- a0 - a1 + a2 - a3 - b1 + b3 + p1 >= -3',
 'a0 - a1 - a2 + a3 + b1 - b3 + p2 >= -2',
 'a0 + a1 + a2 - b0 - b1 + b2 - b3 >= -2',
 'a1 + a2 - b0 - b1 - b2 + b3 >= -2',
 'a0 + a1 - a2 - a3 + b1 + b2 - b3 >= -2',
 'a1 - a2 - a3 + b1 - b2 + b3 >= -2',
 'a0 - a1 + a2 - b0 + b1 - b2 >= -2',
 'a0 - a1 - a2 - a3 - b1 - b2 >= -4',
 '- a0 - a1 - a2 + a3 - b2 >= -3']
```

The call also returns CNF and CP constraints; only the MILP inequalities are shown here.

Interpretation of the outputs:

- `Input: a0||a1||a2||a3; msb: a0`: The binary vector $`a = a_{0}||a_{1}||a_{2}||a_{3}`$ encodes the input difference where $`a_{0}`$ is the most significant bit of $`a`$.
- `Output: b0||b1||b2||b3; msb: b0`: The binary vector $`b = b_{0}||b_{1}||b_{2}||b_{3}`$ encodes the output difference where $`b_{0}`$ is the most significant bit of $`b`$.
- `Weight: 3.0000 p0 + 2.0000 p1 + 1.4150 p2`: The linear function $`3 \cdot p_0 + 2 \cdot p_1 + 1.4150 \cdot p_2`$ encodes the weight of differential transition $`a \rightarrow b`$, where $`\Pr (a \rightarrow b) = 2^{-(3 \cdot p_0 + 2 \cdot p_1 + 1.4150 \cdot p_2)}`$. The additional variables $`p_{0}, p_{1}`$, and $`p_{2}`$ are binary decision variables to encode the probability of differential transitions.

**Note:** Python code blocks can be copied directly into scripts or notebooks. Example output is shown separately. The `sage.crypto.sboxes` imports are provided by PassageMath, which is installed with the package; no SageMath shell is required.

## Convert a Set of Binary Vectors to CNF/MILP Constraints

Here, we show how to describe a set of binary vectors over $`\mathbb{F}_{2}^{n}`$ to (minimized) CNF/MILP constraints.

Assume that we have a set of binary vectors over $`\mathbb{F}_{2}^{17}`$ as follows:

```
S = {
  (1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0),
  (1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0),
  (1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0),
  (1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0),
  (0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0),
  (0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0),
  (0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1)
}
```

We can use S-box Analyzer to encode $`S`$ to (minimized) CNF/MILP constraints.
The following code shows how to do this:

```python
from sboxanalyzer import SboxAnalyzer as SA
S = [
    (1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0),
    (1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0),
    (1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0),
    (1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0),
    (0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0),
    (0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0),
    (0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1)
]
cnf, milp, cp = SA.encode_set_of_binary_vectors(S, mode=6)
print(cnf)
```

Example output:

```text
Number of constraints: 27
Variables: x0||x1||x2||x3||x4||x5||x6||x7||x8||x9||x10||x11||x12||x13||x14||x15||x16; msb: x0
(~x1 | ~x3) & (~x3 | x5) & (x5 | ~x6) & (x6 | ~x7) & (~x2 | ~x9) & (x1 | x9) & (~x6 | x13) & (~x11 | x14) & (x12 | x14) & (x13 | ~x15) & (x8 | x15) & (x3 | ~x16) & (x7 | ~x16) & (x14 | ~x16) & (x3 | x4 | ~x7) & (x3 | x6 | x8) & (~x0 | x4 | ~x9) & (x4 | ~x8 | ~x10) & (~x4 | ~x8 | x11) & (~x4 | ~x10 | ~x13) & (x0 | x7 | ~x14) & (~x6 | x11 | ~x15) & (~x1 | x12 | ~x15) & (~x0 | x3 | x15) & (x10 | ~x12 | x16) & (~x5 | ~x8 | ~x9 | x12) & (x2 | ~x5 | ~x11 | ~x13)
```

## Convert a Truth Table to CNF/MILP Constraints

Here, we show how to convert a truth table of a Boolean function to CNF/MILP constraints.
Let's generate a random truth table of a Boolean function $`f:\mathbb{F}_{2}^{4} \rightarrow \mathbb{F}_{2}`$:

```python
from random import randint
BF = [randint(0, 1) for _ in range(16)]
print(BF)
```

Example output:

```text
[0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0]
```

The randomly generated truth table and its constraints vary between runs.
We can use S-box Analyzer to encode the truth table of $`f`$ to (minimized) CNF/MILP constraints as follows:

```python
from pprint import pprint
from sboxanalyzer import SboxAnalyzer as SA
cnf, milp, cp = SA.encode_boolean_function(BF, mode=6)
pprint(milp)
print(cnf)
```

Example output:

```text
Number of constraints: 3
Variables: x0||x1||x2||x3; msb: x0
['x0 + x3 >= 1', '- x1 - x2 >= -1', '- x1 + x3 >= 0']
(x0 | x3) & (~x1 | ~x2) & (~x1 | x3)
```

## DDT Encoding

### Encoding the DDT of SKINNY-64

```python
from sboxanalyzer import *
from sage.crypto.sboxes import SKINNY_4 as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_diff_constraints()
```

Example output:

```text
Number of constraints: 39
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
Weight: 3.0000 p0 + 2.0000 p1
```

To make a trade-off between the time of simplification and the solution's optimality, S-Box Analyzer supports seven different modes, i.e., `[mode=1,...,mode=7]`. The default mode is 6, which is the best choice for both simplification time and optimality. For example, using the following command, we can minimize the number of constraints a little more:

```python
cnf, milp, cp = sa.minimized_diff_constraints(mode=5)
```

Example output:

```text
Number of constraints: 37
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
Weight: 3.0000 p0 + 2.0000 p1
```

### Encoding the DDT of Ascon

```python
from sage.crypto.sboxes import Ascon as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_diff_constraints()
```

Example output:

```text
Number of constraints: 77
Input: a0||a1||a2||a3||a4; msb: a0
Output: b0||b1||b2||b3||b4; msb: b0
Weight: 4.0000 p0 + 3.0000 p1 + 2.0000 p2
```

### Encoding the DDT of PRESENT

```python
from sage.crypto.sboxes import PRESENT as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_diff_constraints()
```

Example output:

```text
Number of constraints: 54
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
Weight: 3.0000 p0 + 2.0000 p1
```

### Encoding the DDT of SKINNY-128

To encode the DDT of large S-boxes (8-bit S-boxes), we usually divide the DDT into several sub-DDTs and encode the sub-DDTs seperately. Lastly, we put the constraints for all sub-DDTs together to encode the whole DDT. The following code shows how to encode the DDT of SKINNY-128.

***Encode 2-DDT***

```python
from sage.crypto.sboxes import SKINNY_8 as sb
sa = SboxAnalyzer(sb)
print(sa.get_differential_spectrum())
cnf, milp, cp = sa.minimized_diff_constraints(subtable=2, mode=2)
```

Example output:

```text
[2, 4, 6, 8, 12, 16, 20, 24, 28, 32, 40, 48, 64]
Number of constraints: 206
Input: a0||a1||a2||a3||a4||a5||a6||a7; msb: a0
Output: b0||b1||b2||b3||b4||b5||b6||b7; msb: b0
```

***Encode 4-DDT***

```python
cnf, milp, cp = sa.minimized_diff_constraints(subtable=4)
```

Example output:

```text
Number of constraints: 279
Input: a0||a1||a2||a3||a4||a5||a6||a7; msb: a0
Output: b0||b1||b2||b3||b4||b5||b6||b7; msb: b0
```

***Encode 6-DDT***

```python
cnf, milp, cp = sa.minimized_diff_constraints(subtable=6)
```

Example output:

```text
Number of constraints: 33
Input: a0||a1||a2||a3||a4||a5||a6||a7; msb: a0
Output: b0||b1||b2||b3||b4||b5||b6||b7; msb: b0
```

***Encode 8-DDT***

```python
cnf, milp, cp = sa.minimized_diff_constraints(subtable=8)
```

Example output:

```text
Number of constraints: 235
Input: a0||a1||a2||a3||a4||a5||a6||a7; msb: a0
Output: b0||b1||b2||b3||b4||b5||b6||b7; msb: b0
```

You can encode other sub-DDTs of SKINNY-128's S-box in a similar way. Moreover, you may achieve a more optimal solution using the different modes such as `mode=2` or `mode=5`. However, the simplification time might be higher. The default mode is `mode=6` since it generates the nearly optimum result and is sufficient to get a remarkable speed up in automatic differential analysis based on MILP or SAT/SMT.

### Encoding the DDT of AES

```python
from sage.crypto.sboxes import AES as sb
sa = SboxAnalyzer(sb)
print(sa.get_differential_spectrum())
cnf, milp, cp = sa.minimized_diff_constraints(subtable=2)
cnf, milp, cp = sa.minimized_diff_constraints(subtable=4)
```

Example output:

```text
[2, 4]
Number of constraints: 7967
Input: a0||a1||a2||a3||a4||a5||a6||a7; msb: a0
Output: b0||b1||b2||b3||b4||b5||b6||b7; msb: b0
Number of constraints: 321
Input: a0||a1||a2||a3||a4||a5||a6||a7; msb: a0
Output: b0||b1||b2||b3||b4||b5||b6||b7; msb: b0
```

As can be seen our results concerning encoding the DDT of AES's S-box is much better than the results reported in [this paper](https://tosc.iacr.org/index.php/ToSC/article/view/805).

### Encoding the *-DDT and *-LAT

In impossible differential attack (or zero correlation linear attacks) we only encode the possibility of the differential transitions (or a linear transitions), i.e., the *-DDT (or *-LAT). As illustrated in the following example, by setting the `subtable` argument to `star` we can simply encode the *-DDT.

```python
from sboxanalyzer import *
from sage.crypto.sboxes import Midori_Sb0 as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_diff_constraints(subtable="star", mode=5)
```

Example output:

```text
Number of constraints: 47
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
```

### Encoding the DDT for CryptoSMT

By setting the `cryptosmt_compatible` argument to `True`, you can generate an SMT formula compatible with CryptoSMT. For example, to encode the DDT of [CRAFT](https://tosc.iacr.org/index.php/ToSC/article/view/8466) in a format compatible with CryptoSMT, you can use the following commands:

```python
from sboxanalyzer import *
from sage.crypto.sboxes import CRAFT as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_diff_constraints(cryptosmt_compatible=True)
print(cnf)
```

Example output:

```text
Number of constraints: 54
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
Weight: p0 + p1 + p2
(~a2 | p1) & (~b2 | p1) & (~p0 | p1) & (~p1 | p2) & (a1 | ~a2 | a3 | ~p0) & (a1 | ~a3 | b2 | p0) & (a2 | ~a3 | b2 | p0) & (~a1 | a3 | b2 | p0) & (a2 | b1 | ~b3 | p0) & (a2 | ~b1 | b3 | p0) & (~a0 | b2 | b3 | p0) & (a1 | a2 | a3 | ~b1 | ~b3) & (a0 | ~a2 | a3 | ~b2 | b3) & (~a1 | ~a2 | ~a3 | b2 | b3) & (~a1 | ~a3 | b0 | b1 | p0) & (a0 | a1 | ~b1 | ~b3 | p0) & (a1 | ~a3 | ~b1 | ~b3 | p0) & (~a1 | a3 | ~b1 | ~b3 | p0) & (~a0 | a3 | b1 | ~b3 | p0) & (~a1 | ~b0 | b1 | ~b3 | p0) & (a1 | ~a3 | ~b0 | b3 | p0) & (~a3 | ~b0 | ~b1 | b3 | p0) & (a0 | a1 | a2 | a3 | ~p2) & (b0 | b1 | b2 | b3 | ~p2) & (~a1 | ~a2 | ~b0 | b1 | b2 | ~b3) & (~a2 | a3 | b0 | ~b1 | b2 | ~p0) & (~a2 | ~a3 | b0 | b2 | ~b3 | ~p0) & (~a0 | a1 | ~a3 | ~b2 | b3 | ~p0) & (a0 | a3 | ~b0 | b1 | b3 | p0) & (a1 | a2 | a3 | b1 | b3 | ~p2) & (a0 | a3 | ~b0 | b1 | b2 | ~b3 | ~p0) & (~a0 | a1 | a2 | ~a3 | b0 | b3 | ~p0) & (a0 | a1 | ~b0 | ~b1 | b2 | b3 | ~p0) & (~a1 | ~a3 | b1 | b3 | ~p0 | ~p1 | ~p2) & (~a0 | ~a2 | b0 | ~b2 | p0 | ~p1 | ~p2) & (~a0 | ~a1 | a3 | ~b0 | b1 | ~b2 | ~p1 | ~p2) & (~a0 | a1 | ~a2 | ~b0 | ~b1 | b3 | ~p1 | ~p2) & (a0 | a2 | ~a3 | b1 | ~b2 | ~p0 | ~p1 | ~p2) & (a2 | b0 | ~b1 | ~b2 | ~b3 | ~p0 | ~p1 | ~p2) & (a0 | ~a1 | a2 | ~b2 | b3 | ~p0 | ~p1 | ~p2) & (a0 | ~a1 | ~a2 | ~a3 | ~b2 | p0 | ~p1 | ~p2) & (a0 | a1 | ~a2 | b1 | ~b2 | p0 | ~p1 | ~p2) & (~a0 | ~a1 | a2 | a3 | b0 | b1 | ~p0 | ~p1 | ~p2) & (~a2 | a3 | b1 | ~b2 | ~p0) & (a1 | ~a2 | ~b2 | b3 | ~p0) & (~a1 | ~a3 | ~b1 | ~b3 | ~p0) & (a2 | ~b0 | ~b1 | ~b3) & (~a0 | ~a3 | b0 | b1 | b2) & (a1 | ~a2 | b1 | ~b2 | ~p0) & (~a0 | ~a1 | b0 | b2 | b3) & (~a2 | a3 | ~b2 | b3 | ~p0) & (a0 | a1 | a2 | ~b0 | ~b3) & (a0 | a2 | a3 | ~b0 | ~b1) & (~a0 | ~a1 | ~a3 | b2)
```

## Deriving and Encoding Deterministic Differential Propagation

The following example shows how to derive and encode the deterministic differential propagation of the Ascon S-box.

```python
from sboxanalyzer import *
from sage.crypto.sboxes import Ascon as sb
sa = SboxAnalyzer(sb)
ddp = sa.encode_deterministic_differential_behavior()
cp = sa.generate_cp_constraints(ddp)
print(cp)
```

Example output:

```text
Input: a0||a1||a2||a3||a4; msb: a0
Output: b0||b1||b2||b3||b4; msb: b0
if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = 0 /\ b1 = 0 /\ b2 = 0 /\ b3 = 0 /\ b4 = 0)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 1) then (b0 = -1 /\ b1 = 1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 1 /\ a4 == 0) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 1 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = 0 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = 1 /\ b4 = 0)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == 1) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 1 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 1 /\ a4 == 1) then (b0 = 0 /\ b1 = -1 /\ b2 = -1 /\ b3 = 1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == -1 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 0)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == -1 /\ a3 == 1 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 1)
elseif (a0 == 0 /\ a1 == 1 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = 1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 1 /\ a2 == 0 /\ a3 == 1 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = 1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 1 /\ a2 == 1 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 0 /\ b3 = 0 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 1 /\ a2 == 1 /\ a3 == 1 /\ a4 == 0) then (b0 = -1 /\ b1 = 0 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 1 /\ a2 == 1 /\ a3 == 1 /\ a4 == 1) then (b0 = -1 /\ b1 = 1 /\ b2 = -1 /\ b3 = 0 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = 1 /\ b2 = 0 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 1) then (b0 = 1 /\ b1 = 0 /\ b2 = -1 /\ b3 = -1 /\ b4 = 1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 0 /\ a3 == 1 /\ a4 == 1) then (b0 = 0 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 0)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == 0) then (b0 = 0 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 1 /\ a3 == 1 /\ a4 == 0) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 1 /\ a3 == 1 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 0)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == -1 /\ a3 == 0 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == -1 /\ a3 == 1 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 0)
elseif (a0 == 1 /\ a1 == 1 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 1 /\ a2 == 1 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 0 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 1 /\ a2 == 1 /\ a3 == 1 /\ a4 == 0) then (b0 = -1 /\ b1 = 1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 1 /\ a2 == 1 /\ a3 == 1 /\ a4 == 1) then (b0 = -1 /\ b1 = 0 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == -1 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 0 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == -1 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == -1 /\ a1 == 1 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == -1 /\ a1 == 1 /\ a2 == 1 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = 0 /\ b3 = -1 /\ b4 = -1)
else (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
endif
```

## LAT Encoding

Here, we show how to encode the (squared) LAT of S-boxes.

### Encoding the LAT of SKINNY-64

```python
from sboxanalyzer import *
from sage.crypto.sboxes import SKINNY_4 as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_linear_constraints()
```

Example output:

```text
Number of constraints: 29
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
Weight: 4.0000 p0 + 2.0000 p1
```

Note that `sa.minimized_linear_constraints()` encodes squared correlations, rather than signed LAT entries.

### Encoding the LAT of Ascon

```python
from sage.crypto.sboxes import Ascon as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_linear_constraints()
```

Example output:

```text
Number of constraints: 96
Input: a0||a1||a2||a3||a4; msb: a0
Output: b0||b1||b2||b3||b4; msb: b0
Weight: 4.0000 p0 + 2.0000 p1
```

## Deriving and Encoding Deterministic Linear Propagation

The following example shows how to derive and encode the deterministic linear propagation through the inverse of the Ascon S-box.

```python
from sboxanalyzer import *
from sage.crypto.sboxes import Ascon as sb
sa = SboxAnalyzer(sb.inverse())
dlp = sa.encode_deterministic_linear_behavior()
cp = sa.generate_cp_constraints(dlp)
print(cp)
```

Example output:

```text
Input: a0||a1||a2||a3||a4; msb: a0
Output: b0||b1||b2||b3||b4; msb: b0
if (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = 0 /\ b1 = 0 /\ b2 = 0 /\ b3 = 0 /\ b4 = 0)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = 0 /\ b3 = 1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == -1) then (b0 = -1 /\ b1 = -1 /\ b2 = 0 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 1 /\ a4 == 0) then (b0 = -1 /\ b1 = 1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 1 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 0 /\ a3 == 1 /\ a4 == -1) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == 0) then (b0 = 0 /\ b1 = 1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == -1) then (b0 = -1 /\ b1 = -1 /\ b2 = 1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 1 /\ a4 == 0) then (b0 = -1 /\ b1 = 0 /\ b2 = 0 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 1 /\ a4 == 1) then (b0 = -1 /\ b1 = -1 /\ b2 = 0 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == 1 /\ a3 == 1 /\ a4 == -1) then (b0 = -1 /\ b1 = -1 /\ b2 = 0 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 0 /\ a2 == -1 /\ a3 == 0 /\ a4 == 0) then (b0 = 0 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 1 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 1)
elseif (a0 == 0 /\ a1 == 1 /\ a2 == 1 /\ a3 == 0 /\ a4 == 0) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 0 /\ a1 == 1 /\ a2 == -1 /\ a3 == 0 /\ a4 == 0) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 0) then (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = 1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 0 /\ a3 == 0 /\ a4 == 1) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = 0 /\ b4 = 1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == 1 /\ a3 == 0 /\ a4 == 1) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 0 /\ a2 == -1 /\ a3 == 0 /\ a4 == 1) then (b0 = 1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 1 /\ a2 == 0 /\ a3 == 0 /\ a4 == 1) then (b0 = 0 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = 0)
elseif (a0 == 1 /\ a1 == 1 /\ a2 == 1 /\ a3 == 0 /\ a4 == 1) then (b0 = 0 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
elseif (a0 == 1 /\ a1 == 1 /\ a2 == -1 /\ a3 == 0 /\ a4 == 1) then (b0 = 0 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
else (b0 = -1 /\ b1 = -1 /\ b2 = -1 /\ b3 = -1 /\ b4 = -1)
endif
```

## MPT Encoding

Here, we show how to encode the propagation of monomial trails through S-boxes.

### Encoding the MPT of PRESENT

```python
from pprint import pprint
from sboxanalyzer import *
from sage.crypto.sboxes import PRESENT as sb
sa = SboxAnalyzer(sb)
mpt = sa.monomial_prediction_table()
pprint(mpt)
cnf, milp, cp = sa.minimized_integral_constraints()
```

Example output:

```text
[[1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
 [0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
 [0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0],
 [0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0],
 [0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
 [0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
 [0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0],
 [0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0],
 [0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
 [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
 [0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1],
 [0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1],
 [0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
 [0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0],
 [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]
Number of constraints: 41
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
```

### Encoding the MPT of Ascon

```python
from sboxanalyzer import *
from sage.crypto.sboxes import Ascon as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_integral_constraints()
```

Example output:

```text
Number of constraints: 97
Input: a0||a1||a2||a3||a4; msb: a0
Output: b0||b1||b2||b3||b4; msb: b0
```

## Convex Hull Constraints for MILP

By default, the MILP constraints are the clauses minimized by ESPRESSO, each written as a linear inequality. With `milp_method="convex_hull"`, the MILP constraints are instead the facets of the convex hull of the valid points:

```python
from pprint import pprint
from sboxanalyzer import *
from sage.crypto.sboxes import PRESENT as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_integral_constraints(milp_method="convex_hull")
pprint(milp[:2])
```

Example output:

```text
Number of constraints: 41
Number of convex hull constraints: 209
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
['- 2 a0 - a1 - a2 - 2 a3 + 5 b0 + 5 b1 + 5 b2 + 2 b3 >= 0',
 '- 2 a0 - a1 - a2 + 3 b0 + 3 b1 + 3 b2 + 2 b3 >= 0']
```

Both sets of constraints accept exactly the same binary points. They differ in their linear relaxation: the convex hull gives the tightest relaxation that an exact set of linear inequalities can have, at the cost of more constraints. A tighter relaxation can help a MILP solver, in particular when it has to prove that a model is infeasible, whereas more constraints make each solver iteration more expensive. Which one is faster depends on the model, so it is worth comparing both.

- Only the MILP constraints (the second output) change. The SAT clauses and CP constraints are still derived from ESPRESSO, and `mode` only affects those.
- The coefficients are integers that are not necessarily in $`\{-1, 1\}`$. The bounds $`0 \leq x \leq 1`$ are not included, since the variables are supposed to be binary.
- It is supported by `minimized_diff_constraints`, `minimized_linear_constraints`, `minimized_integral_constraints`, and `encode_set_of_binary_vectors`.
- It requires `passagemath-polyhedra`, which is part of PassageMath Standard and can also be installed using `pip install "sboxanalyzer[hull]"`.
- It is limited to 12 binary variables (`SboxAnalyzer.convex_hull_max_variables`), which is enough for the DDT, LAT, and MPT of the common 4-bit S-boxes. The number of facets grows quickly with the number of variables: the convex hull of the MPT of PRESENT has 209 facets, whereas that of its DDT (with the probability variables) has 1645.

## DLCT Encoding

Here, we show how to encode the DLCT of S-boxes.
Use `differential_linear_connectivity_table()` to compute the DLCT and `minimized_differential_linear_constraints(subtable="star")` to encode its nonzero support. The UDLCT and LDLCT are computed by `upper_differential_linear_connectivity_table()` and `lower_differential_linear_connectivity_table()`, respectively.

### Encoding the *-DLCT of KNOT

```python
from sboxanalyzer import *
from sage.crypto.sboxes import KNOT as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_differential_linear_constraints(subtable='star')
```

Example output:

```text
Number of constraints: 34
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
```

### Encoding the *-DLCT of Midori

```python
from sboxanalyzer import *
from sage.crypto.sboxes import Midori_Sb0 as sb
sa = SboxAnalyzer(sb)
cnf, milp, cp = sa.minimized_differential_linear_constraints(subtable='star')
```

Example output:

```text
Number of constraints: 43
Input: a0||a1||a2||a3; msb: a0
Output: b0||b1||b2||b3; msb: b0
```

### Verification of Hadipour et al.'s Theorem

The following code verifies the Hadipour et al.'s theorem, (Proposition 2 in [2024/255](https://ia.cr/2024/255)) for the S-box of Ascon.

```python
from sboxanalyzer import *
from sage.crypto.sboxes import Ascon as sb
sa = SboxAnalyzer(sb)
check = sa.check_hadipour_theorem()
print(check)
```

Example output:

```text
The Hadipour et al.'s theorem is satisfied.
True
```

## BCT Encoding

For a permutation S-box, `boomerang_connectivity_table()` computes the BCT, and `minimized_boomerang_constraints(subtable="star")` encodes its nonzero support:

```python
from sboxanalyzer import SboxAnalyzer
from sage.crypto.sboxes import PRESENT

sa = SboxAnalyzer(PRESENT)
cnf, milp, cp = sa.minimized_boomerang_constraints(subtable="star")
```

This BCT encoding uses input and output differences. The [exact boomerang pair-relation interface](#exact-boomerang-pair-relation-encoding) additionally represents value masks in the chosen bases.

## Exact Pair-Relation Encoding

The following examples generate local S-box constraints for quasidifferential and mix-basis trails. These relations describe both a value mask and a difference (or a label in the selected difference basis), rather than a difference alone.

Use `minimized_exact_pair_constraints()` and choose the input and output bases explicitly. The supported names are `quasidifferential` (a character on the value and one fixed difference), `linear` (characters on both the value and the difference), `subset-difference`, and `ultrametric-difference`. The last two combine a value character with subset or ultrametric coordinates on the difference. Call `sa.supported_exact_pair_bases()` to list them.

### Quasidifferential Trail Constraints

We use the 3-bit PRINTcipher S-box to keep the examples small. Selecting `quasidifferential` at both ends gives the local support constraints for a quasidifferential trail:

```python
from sboxanalyzer import SboxAnalyzer
from sage.crypto.sboxes import PRINTcipher

sa = SboxAnalyzer(PRINTcipher)
m = sa.input_size()
q0 = [f"q0_{i}" for i in range(2 * m)]
q1 = [f"q1_{i}" for i in range(2 * m)]
cnf_q, milp_q, cp_q = sa.minimized_exact_pair_constraints(
    input_basis="quasidifferential",
    output_basis="quasidifferential",
    subtable="star",
    input_variables=q0,
    output_variables=q1,
)
```

The return values are SAT clauses, MILP inequalities, and CP constraints. For an $`m`$-bit input and $`n`$-bit output, the method uses $`2m`$ input and $`2n`$ output variables. In quasidifferential coordinates, the difference bits come first, followed by the value-mask bits, each most-significant bit first. Thus `q0[:m]` describes the input difference and `q0[m:]` its value mask. Unlike a DDT-only encoding, this relation also keeps the value masks.

### Mix-Basis Trail Constraints

Continuing with the same S-box, we switch from quasidifferential input coordinates to linear output coordinates:

```python
l2 = [f"l2_{i}" for i in range(2 * m)]
cnf_mix, milp_mix, cp_mix = sa.minimized_exact_pair_constraints(
    input_basis="quasidifferential",
    output_basis="linear",
    subtable="star",
    input_variables=q1,
    output_variables=l2,
)
```

Here the first half of `l2` is a mask on the output difference, and the second half is a mask on the output value. To switch to ultrametric coordinates on the difference instead, use `output_basis="ultrametric-difference"`.

For two S-boxes connected directly, the shared variables `q1` connect the first example to this switch. Their combined support constraints are:

```python
cnf_trail = f"({cnf_q}) & ({cnf_mix})"
milp_trail = milp_q + milp_mix
```

This models trails through two successive S-box applications, with quasidifferential coordinates at their connection and linear coordinates at the output. For a complete cipher, add its linear-layer, key-addition, and endpoint constraints. Shared labels must use the same basis, and all MILP label variables are binary.

### Valuation-Weighted Constraints

To encode local 2-adic valuations as well as support, use `subtable=None` and `weight_model="v2"`:

```python
cnf_v2, milp_v2, cp_v2 = sa.minimized_exact_pair_constraints(
    input_basis="quasidifferential",
    output_basis="linear",
    subtable=None,
    weight_model="v2",
)
print(sa.exact_pair_objective)
```

Set `output_basis="quasidifferential"` in this call for quasidifferential-to-quasidifferential valuation constraints. Each nonzero integer table entry $`R`$ receives cost $`v_2(R)`$, encoded by one-hot cost variables. For these quasidifferential-input examples, the normalized transition coefficient is $`R/2^m`$, so its valuation is $`v_2(R)-m`$. These costs are not negative logarithms of probabilities. When combining weighted components, give each component distinct `cost_variables` and account for the normalization in the trail objective.

Use `exact_pair_transition_table()` to obtain the integer table, `get_exact_pair_transition_spectrum()` for its nonzero entries, and `get_exact_pair_transition_v2_spectrum()` for their valuation classes. These methods accept the same `input_basis` and `output_basis` arguments.

### Exact Boomerang Pair-Relation Encoding

The same interface is available for the exact boomerang relation of a permutation S-box. For every input pair state it counts, with the signs of the chosen bases, the output pair states whose boomerang quartet returns through the S-box. Restricted to zero value masks, the `quasidifferential`-to-`quasidifferential` table coincides with the BCT.

```python
from sboxanalyzer import SboxAnalyzer
from sage.crypto.sboxes import PRESENT
sa = SboxAnalyzer(PRESENT)
print(sa.supported_exact_boomerang_bases())
cnf, milp, cp = sa.minimized_exact_boomerang_constraints(
    input_basis="quasidifferential",
    output_basis="quasidifferential",
    subtable="star",
    mode=6,
)
```

Example output:

```text
('quasidifferential', 'linear', 'subset-difference', 'ultrametric-difference')
Number of constraints: 1570
Input boomerang pair-state: a0||a1||a2||a3||a4||a5||a6||a7; msb: a0
Output boomerang pair-state: b0||b1||b2||b3||b4||b5||b6||b7; msb: b0
```

The options `subtable` and `weight_model` have the same meaning as for `minimized_exact_pair_constraints`.
The corresponding table and spectrum methods are `exact_boomerang_transition_table()`, `get_exact_boomerang_transition_spectrum()`, and `get_exact_boomerang_transition_v2_spectrum()`.

## Citation

If you use our tools in a project resulting in an academic publication, please acknowledge it by citing our papers:

```
@article{DBLP:journals/tosc/HadipourNE22,
  author    = {Hosein Hadipour and
               Marcel Nageler and
               Maria Eichlseder},
  title     = {Throwing Boomerangs into Feistel Structures Application to CLEFIA,
               WARP, LBlock, LBlock-s and {TWINE}},
  journal   = {IACR Trans. Symmetric Cryptol.},
  volume    = {2022},
  number    = {3},
  pages     = {271--302},
  year      = {2022},
  doi       = {10.46586/tosc.v2022.i3.271-302}
}

@misc{cryptoeprint:2026/2471,
  author       = {Hosein Hadipour},
  title        = {Connecting Cryptanalytic Distinguishers through Divisibility},
  howpublished = {Cryptology {ePrint} Archive, Paper 2026/2471},
  year         = {2026},
  url          = {https://eprint.iacr.org/2026/2471}
}
```

## Papers

Sbox Analyzer has been used in the following papers:

- [Throwing Boomerangs into Feistel Structures: Application to CLEFIA, WARP, LBlock, LBlock-s and TWINE](https://ia.cr/2022/745)
- [Integral Cryptanalysis of WARP based on Monomial Prediction](https://ia.cr/2022/729)
- [Improved Search for Integral, Impossible-Differential and Zero-Correlation Attacks: Application to Ascon, ForkSKINNY, SKINNY, MANTIS, PRESENT and QARMAv2](https://ia.cr/2023/1701)
- [Revisiting Differential-Linear Attacks via a Boomerang Perspective With Application to AES, Ascon, CLEFIA, SKINNY, PRESENT, KNOT, TWINE, WARP, LBlock, Simeck, and SERPENT](https://ia.cr/2024/255)
- [Gleeok: A Family of Low-Latency PRFs andits Applications to Authenticated Encryption](https://tches.iacr.org/index.php/TCHES/article/view/11439/10944)
- [Breaking the Twinkle Authenticated Encryption Scheme and Analyzing Its Underlying Permutation](https://ia.cr/2025/1339)
- [On Extending Integral Distinguishers](https://ia.cr/2026/1402)
- [Connecting Cryptanalytic Distinguishers through Divisibility](https://ia.cr/2026/2471)

## Road Map

The function names below are public methods of `SboxAnalyzer`.

- [X] Encoding DDT — `minimized_diff_constraints()`
- [X] Encoding LAT support and squared correlations — `minimized_linear_constraints()`
- [X] Encoding [MPT](https://tosc.iacr.org/index.php/ToSC/article/view/9715) — `minimized_integral_constraints()`
- [X] Encoding DLCT support — `minimized_differential_linear_constraints(subtable="star")`
- [X] Computing UDLCT and LDLCT — `upper_differential_linear_connectivity_table()` and `lower_differential_linear_connectivity_table()`
- [X] Checking [Hadipour et al.'s theorem](https://eprint.iacr.org/2024/255) — `check_hadipour_theorem()`
- [X] Encoding BCT support — `minimized_boomerang_constraints(subtable="star")`
- [X] Encoding exact pair-state support and local 2-adic valuation constraints — `minimized_exact_pair_constraints()`
- [X] Encoding exact boomerang pair-state support and local 2-adic valuation constraints — `minimized_exact_boomerang_constraints()`
- [X] Deriving MILP constraints from the convex hull — `milp_method="convex_hull"` in the [supported methods](#convex-hull-constraints-for-milp)
- [X] Delivering a pip-installable Python package
- [ ] Full integration into [CryptoSMT](https://github.com/kste/cryptosmt) (DDT encoding with `cryptosmt_compatible=True` is already available)

Any contributions or comments regarding the development of the tool are warmly welcome.

## License

S-box Analyzer is released under the [MIT license](LICENSE).

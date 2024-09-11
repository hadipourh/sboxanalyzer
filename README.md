# S-box Analyzer
[![license](./images/license-MIT-green.svg)](./LICENSE.txt)


**S-box Analyzer**

S-box Analyzer is a tool for analyzing S-boxes and (vectorial) Boolean functions against differential, linear, differential-linear, boomerang, and integral attacks. 

Particularly, it derives the minimized CP/MILP and SMT/SAT constraints to encode a set of binary vectors, (vectorial) Boolean functions, the Differential Distribution Table (DDT), Linear Approximation Table (LAT), [Differential-Linear Connectivity Table (DLCT)](https://ia.cr/2024/255) and [Monomial Prediction Table (MPT)](https://tosc.iacr.org/index.php/ToSC/article/view/9715) of S-boxes.

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
  - [DLCT Encoding](#dlct-encoding)
    - [Encoding the \*-DLCT of KNOT](#encoding-the--dlct-of-knot)
    - [Encoding the \*-DLCT of Midori](#encoding-the--dlct-of-midori)
    - [Verification of Hadipour et al.'s Theorem](#verification-of-hadipour-et-als-theorem)
  - [Citation](#citation)
  - [Papers](#papers)
  - [Road Map](#road-map)
  - [License](#license)


---

## Dependencies

S-box Analyzer has been implemented as a SageMath module and employs ESPRESSO for logic minimization. Thus, it requires the following software:

- [ESPRESSO](https://ptolemy.berkeley.edu/projects/embedded/pubs/downloads/espresso/index.htm)
- [SageMath](https://www.sagemath.org/)

## Installation

- ESPRESSO

  To build ESPRESSO, navigate into the espresso folder and run the following commands:

  ```bash
  cd espresso
  mkdir build && cd build
  cmake ..
  make -j8
  ```

  To see more details about the ESPRESSO logic minimizer, see [here](https://ptolemy.berkeley.edu/projects/embedded/pubs/downloads/espresso/index.htm). The updated versions of ESPRESSO, compatible with new compilers, are available [here](https://github.com/classabbyamp/espresso-logic) and [here](https://github.com/chipsalliance/espresso/releases/tag/v2.4).

- SageMath

  To see the installation recipe of SageMath, see [here](https://doc.sagemath.org/html/en/installation/index.html).

## Usage

Run the SageMath in the same directory as the S-box Analyzer. Next, import `sboxanalyzer` into the SageMath and then simply use it. The following example shows how to use the S-box Analyzer to encode the DDT of [GIFT](https://giftcipher.github.io/gift/)'s S-box:

```python
sage: from sboxanalyzer import *  
sage: from sage.crypto.sboxes import GIFT as sb
sage: sa = SboxAnalyzer(sb)                               
sage: cnf, milp = sa.minimized_diff_constraints()
                                                       
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 49
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1 + 1.4150 p2


sage: print(milp)
['- a0 - a1 - a2 + a3 - b2 >= -3'
'a0 - a1 - a2 - a3 - b1 - b2 >= -4'
'a0 - a1 + a2 - b0 + b1 - b2 >= -2'
'a1 + a2 - b0 - b1 - b2 + b3 >= -2'
'a1 - a2 - a3 + b1 - b2 + b3 >= -2'
'a0 + a1 + a2 - b0 - b1 + b2 - b3 >= -2'
'a0 + a1 - a2 - a3 + b1 + b2 - b3 >= -2'
'- p0 - p1 >= -1'
'- p0 - p2 >= -1'
'- p1 - p2 >= -1'
'a1 - a3 + p0 >= 0'
'- b0 + b2 + p0 >= 0'
'- b0 + b3 + p0 >= 0'
'a2 + b1 - p2 >= 0'
'b0 + b2 + b3 - p0 >= 0'
'a0 - a3 + b0 + p0 >= 0'
'- a0 - a3 + b1 + p0 >= -1'
'a1 + a2 + b2 - p1 >= 0'
'- a0 + a3 + b0 + p1 >= 0'
'- a0 - a1 + a2 + a3 - b3 >= -2'
'a1 - a2 + a3 - b2 - b3 >= -2'
'a0 + a1 + b0 - b2 - b3 >= -1'
'- a1 - a3 + b0 - b2 - b3 >= -3'
'- a0 + b0 - b1 + b2 - b3 >= -2'
'a0 + b0 + b1 + b2 - b3 >= 0'
'a1 + a2 + a3 - b1 + b3 >= 0'
'a1 + b0 + b1 - b2 + b3 >= 0'
'a0 - a1 + a3 + b2 + b3 >= 0'
'a1 - a2 + a3 + b2 + b3 >= 0'
'- a0 + b0 + b1 + b2 + b3 >= 0'
'a0 - a1 - a2 - b1 + p0 >= -2'
'- a0 - a1 - b1 - b2 + p0 >= -3'
'a1 + a2 + a3 - b0 + p1 >= 0'
'a3 + b0 + b2 - b3 + p1 >= 0'
'- a1 + b0 - b1 + b3 + p1 >= -1'
'a3 + b0 - b2 + b3 + p1 >= 0'
'a0 - a1 - a2 - b0 - b1 - b3 >= -4'
'a0 - a1 + a2 - a3 + b1 - b3 >= -2'
'- a0 - a2 - a3 - b1 + b2 - b3 >= -4'
'- a0 + a2 - b0 + b1 + b2 - b3 >= -2'
'a0 - a1 - b0 - b2 - b3 + p1 >= -3'
'- a0 + a1 - b0 - b2 - b3 + p1 >= -3'
'- a0 - a1 - a3 + b2 + b3 + p1 >= -2'
'a0 + a2 + a3 - b1 - b2 + p2 >= -1'
'a0 + a2 + a3 - b2 + p0 + p2 >= 0'
'- a0 - a1 - a2 - a3 - b0 + b1 + b3 >= -4'
'- a0 - a1 + a2 - a3 - b1 - b2 + b3 >= -4'
'a0 - a1 - a2 + a3 + b1 - b3 + p2 >= -2']

sage: print(cnf)
(~a0 | ~a1 | ~a2 | a3 | ~b2) & (a0 | ~a1 | ~a2 | ~a3 | ~b1 | ~b2) & (a0 | ~a1 | a2 | ~b0 | b1 | ~b2) & (a1 | a2 | ~b0 | ~b1 | ~b2 | b3) & (a1 | ~a2 | ~a3 | b1 | ~b2 | b3) & (a0 | a1 | a2 | ~b0 | ~b1 | b2 | ~b3) & (a0 | a1 | ~a2 | ~a3 | b1 | b2 | ~b3) & (~p0 | ~p1) & (~p0 | ~p2) & (~p1 | ~p2) & (a1 | ~a3 | p0) & (~b0 | b2 | p0) & (~b0 | b3 | p0) & (a2 | b1 | ~p2) & (b0 | b2 | b3 | ~p0) & (a0 | ~a3 | b0 | p0) & (~a0 | ~a3 | b1 | p0) & (a1 | a2 | b2 | ~p1) & (~a0 | a3 | b0 | p1) & (~a0 | ~a1 | a2 | a3 | ~b3) & (a1 | ~a2 | a3 | ~b2 | ~b3) & (a0 | a1 | b0 | ~b2 | ~b3) & (~a1 | ~a3 | b0 | ~b2 | ~b3) & (~a0 | b0 | ~b1 | b2 | ~b3) & (a0 | b0 | b1 | b2 | ~b3) & (a1 | a2 | a3 | ~b1 | b3) & (a1 | b0 | b1 | ~b2 | b3) & (a0 | ~a1 | a3 | b2 | b3) & (a1 | ~a2 | a3 | b2 | b3) & (~a0 | b0 | b1 | b2 | b3) & (a0 | ~a1 | ~a2 | ~b1 | p0) & (~a0 | ~a1 | ~b1 | ~b2 | p0) & (a1 | a2 | a3 | ~b0 | p1) & (a3 | b0 | b2 | ~b3 | p1) & (~a1 | b0 | ~b1 | b3 | p1) & (a3 | b0 | ~b2 | b3 | p1) & (a0 | ~a1 | ~a2 | ~b0 | ~b1 | ~b3) & (a0 | ~a1 | a2 | ~a3 | b1 | ~b3) & (~a0 | ~a2 | ~a3 | ~b1 | b2 | ~b3) & (~a0 | a2 | ~b0 | b1 | b2 | ~b3) & (a0 | ~a1 | ~b0 | ~b2 | ~b3 | p1) & (~a0 | a1 | ~b0 | ~b2 | ~b3 | p1) & (~a0 | ~a1 | ~a3 | b2 | b3 | p1) & (a0 | a2 | a3 | ~b1 | ~b2 | p2) & (a0 | a2 | a3 | ~b2 | p0 | p2) & (~a0 | ~a1 | ~a2 | ~a3 | ~b0 | b1 | b3) & (~a0 | ~a1 | a2 | ~a3 | ~b1 | ~b2 | b3) & (a0 | ~a1 | ~a2 | a3 | b1 | ~b3 | p2)
```

Interpretation of the outputs:

- `Input:	a0||a1||a2||a3; a0: msb`: The binary vector $a = a_{0}||a_{1}||a_{2}||a_{3}$ encodes the input difference where $a_{0}$ is the most significant bit of $a$.
- `Output:	b0||b1||b2||b3; b0: msb`: The binary vector $b = b_{0}||b_{1}||b_{2}||b_{3}$ encodes the output difference where $b_{0}$ is the most significant bit of $b$.
- `Weight: 3.0000 p0 + 2.0000 p1 + 1.4150 p2`: The linear function $3 \cdot p_0 + 2 \cdot p_1 + 1.4150 \cdot p_2$ encodes the weight of differential transition $a \rightarrow b$, where $\Pr (a \rightarrow b) = 2^{-(3 \cdot p_0 + 2 \cdot p_1 + 1.4150 \cdot p_2)}$. The additional variables $p_{0}, p_{1}$, and $p_{2}$ are binary decision variables to encode the probability of differential transitions.

## Convert a Set of Binary Vectors to CNF/MILP Constraints

Here, we show how to describe a set of binary vectors over $\mathbb{F}_{2}^{n}$ to (minimized) CNF/MILP constraints.

Assume that we have a set of binary vectors over $\mathbb{F}_{2}^{17}$ as follows:

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

We can use S-box Analyzer to encode $S$ to (minimized) CNF/MILP constraints. 
The following code shows how to do this:

```python
sage: from sboxanalyzer import SboxAnalyzer as SA
sage: S = [
....:     (1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0),
....:     (1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0),
....:     (1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0),
....:     (1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0),
....:     (0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0),
....:     (0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0),
....:     (0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1)
....: ]
sage: cnf, milp = SA.encode_set_of_binary_vectors(S, mode=6)
Generateing and simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.12 seconds
Number of constraints: 27
Variables: x0||x1||x2||x3||x4||x5||x6||x7||x8||x9||x10||x11||x12||x13||x14||x15||x16; msb: x0
sage: pretty_print(milp)
['- x1 - x3 >= -1',
 '- x3 + x5 >= 0',
 'x5 - x6 >= 0',
 'x6 - x7 >= 0',
 '- x2 - x9 >= -1',
 'x1 + x9 >= 1',
 '- x6 + x13 >= 0',
 '- x11 + x14 >= 0',
 'x12 + x14 >= 1',
 'x13 - x15 >= 0',
 'x8 + x15 >= 1',
 'x3 - x16 >= 0',
 'x7 - x16 >= 0',
 'x14 - x16 >= 0',
 'x3 + x4 - x7 >= 0',
 'x3 + x6 + x8 >= 1',
 '- x0 + x4 - x9 >= -1',
 'x4 - x8 - x10 >= -1',
 '- x4 - x8 + x11 >= -1',
 '- x4 - x10 - x13 >= -2',
 'x0 + x7 - x14 >= 0',
 '- x6 + x11 - x15 >= -1',
 '- x1 + x12 - x15 >= -1',
 '- x0 + x3 + x15 >= 0',
 'x10 - x12 + x16 >= 0',
 '- x5 - x8 - x9 + x12 >= -2',
 'x2 - x5 - x11 - x13 >= -2']
sage: print(cnf)
(~x1 | ~x3) & (~x3 | x5) & (x5 | ~x6) & (x6 | ~x7) & (~x2 | ~x9) & (x1 | x9) & (~x6 | x13) & (~x11 | x14) & (x12 | x14) & (x13 | ~x15) & (x8 | x15) & (x3 | ~x16) & (x7 | ~x16) & (x14 | ~x16) & (x3 | x4 | ~x7) & (x3 | x6 | x8) & (~x0 | x4 | ~x9) & (x4 | ~x8 | ~x10) & (~x4 | ~x8 | x11) & (~x4 | ~x10 | ~x13) & (x0 | x7 | ~x14) & (~x6 | x11 | ~x15) & (~x1 | x12 | ~x15) & (~x0 | x3 | x15) & (x10 | ~x12 | x16) & (~x5 | ~x8 | ~x9 | x12) & (x2 | ~x5 | ~x11 | ~x13)
```

## Convert a Truth Table to CNF/MILP Constraints

Here, we show how to convert a truth table of a Boolean function to CNF/MILP constraints.
Let's generate a random truth table of a Boolean function $f:\mathbb{F}_{2}^{4} \rightarrow \mathbb{F}_{2}$:

```python
sage: BF = [randint(0, 1) for _ in range(16)]; BF
[1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0]
```

We can use S-box Analyzer to encode the truth table of $f$ to (minimized) CNF/MILP constraints as follows:

```python
sage: from sboxanalyzer import SboxAnalyzer as SA
sage: cnf, milp = SA.encode_boolean_function(BF, mode=6)
Generateing and simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 5
Variables: x0||x1||x2||x3; msb: x0
sage: pretty_print(milp)
['- x0 - x2 + x3 >= -1',
 'x0 - x1 + x2 + x3 >= 0',
 '- x0 + x1 + x2 - x3 >= -1',
 'x1 - x2 + x3 >= 0',
 '- x1 - x2 - x3 >= -2']
sage: pretty_print(cnf)
'(~x0 | ~x2 | x3) & (x0 | ~x1 | x2 | x3) & (~x0 | x1 | x2 | ~x3) & (x1 | ~x2 | x3) & (~x1 | ~x2 | ~x3)'
```

## DDT Encoding

### Encoding the DDT of SKINNY-64

```python
sage: from sboxanalyzer import *                                               
sage: from sage.crypto.sboxes import SKINNY_4 as sb                            
sage: sa = SboxAnalyzer(sb)                                                    
sage: cnf, milp = sa.minimized_diff_constraints()

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 39
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
```

To make a trade-off between the time of simplification and the solution's optimality, S-Box Analyzer supports seven different modes, i.e., `[mode=1,...,mode=7]`. The default mode is 6, which is the best choice for both simplification time and optimality. For example, using the following command, we can minimize the number of constraints a little more:

```python
sage: cnf, milp = sa.minimized_diff_constraints(mode=5)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.22 seconds
Number of constraints: 37
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
```

### Encoding the DDT of Ascon

```python
sage: from sage.crypto.sboxes import Ascon as sb                               
sage: sa = SboxAnalyzer(sb)                                                    
sage: cnf, milp = sa.minimized_diff_constraints()

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.04 seconds
Number of constraints: 77
Input:	a0||a1||a2||a3||a4; a0: msb
Output:	b0||b1||b2||b3||b4; b0: msb
Weight: 4.0000 p0 + 3.0000 p1 + 2.0000 p2
```

### Encoding the DDT of PRESENT

```python
sage: from sage.crypto.sboxes import PRESENT as sb                             
sage: sa = SboxAnalyzer(sb)                                                    
sage: cnf, milp = sa.minimized_diff_constraints()

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.03 seconds
Number of constraints: 54
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
```

### Encoding the DDT of SKINNY-128

To encode the DDT of large S-boxes (8-bit S-boxes), we usually divide the DDT into several sub-DDTs and encode the sub-DDTs seperately. Lastly, we put the constraints for all sub-DDTs together to encode the whole DDT. The following code shows how to encode the DDT of SKINNY-128.

***Encode 2-DDT***

```python
sage: from sage.crypto.sboxes import SKINNY_8 as sb                            
sage: sa = SboxAnalyzer(sb)                                                    
sage: sa.get_differential_spectrum()                                                     
[2, 4, 6, 8, 12, 16, 20, 24, 28, 32, 40, 48, 64]

sage: cnf, milp = sa.minimized_diff_constraints(subtable=2, mode=2)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.50 seconds
Number of constraints: 206
Input:	a0||a1||a2||a3||a4||a5||a6||a7; a0: msb
Output:	b0||b1||b2||b3||b4||b5||b6||b7; b0: msb
```

***Encode 4-DDT***

```python
sage: cnf, milp = sa.minimized_diff_constraints(subtable=4)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.67 seconds
Number of constraints: 279
Input:	a0||a1||a2||a3||a4||a5||a6||a7; a0: msb
Output:	b0||b1||b2||b3||b4||b5||b6||b7; b0: msb
```

***Encode 6-DDT***

```python
sage: cnf, milp = sa.minimized_diff_constraints(subtable=6)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.03 seconds
Number of constraints: 33
Input:	a0||a1||a2||a3||a4||a5||a6||a7; a0: msb
Output:	b0||b1||b2||b3||b4||b5||b6||b7; b0: msb
```

***Encode 8-DDT***

```python
sage: cnf, milp = sa.minimized_diff_constraints(subtable=8)                    
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.75 seconds
Number of constraints: 235
Input:	a0||a1||a2||a3||a4||a5||a6||a7; a0: msb
Output:	b0||b1||b2||b3||b4||b5||b6||b7; b0: msb
```
You can encode other sub-DDTs of SKINNY-128's S-box in a similar way. Moreover, you may achieve a more optimal solution using the different modes such as `mode=2` or `mode=5`. However, the simplification time might be higher. The default mode is `mode=6` since it generates the nearly optimum result and is sufficient to get a remarkable speed up in automatic differential analysis based on MILP or SAT/SMT.

### Encoding the DDT of AES

```python
sage: from sage.crypto.sboxes import AES as sb                                 
sage: sa = SboxAnalyzer(sb)                                                    
sage: sa.get_differential_spectrum()                                                         
[2, 4]

sage: cnf, milp = sa.minimized_diff_constraints(subtable=2)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 65.11 seconds
Number of constraints: 7967
Input:	a0||a1||a2||a3||a4||a5||a6||a7; a0: msb
Output:	b0||b1||b2||b3||b4||b5||b6||b7; b0: msb

sage: cnf, milp = sa.minimized_diff_constraints(subtable=4)                    
Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 869.08 seconds
Number of constraints: 321
Input:	a0||a1||a2||a3||a4||a5||a6||a7; a0: msb
Output:	b0||b1||b2||b3||b4||b5||b6||b7; b0: msb
```

As can be seen our results concerning encoding the DDT of AES's S-box is much better than the results reported in [this paper](https://tosc.iacr.org/index.php/ToSC/article/view/805).

### Encoding the *-DDT and *-LAT

In impossible differential attack (or zero correlation linear attacks) we only encode the possibility of the differential transitions (or a linear transitions), i.e., the *-DDT (or *-LAT). As illustrated in the following example, by setting the `subtable` argument to `star` we can simply encode the *-DDT.

```python
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import Midori_Sb0 as sb
sage: sa = SboxAnalyzer(sb)
sage: cnf, milp = sa.minimized_diff_constraints(subtable="star", mode=5)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 47
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
```

### Encoding the DDT for CryptoSMT

By setting the `cryptosmt_compatible` argument to `True`, you can generate an SMT formula compatible with CryptoSMT. For example, to encode the DDT of [CRAFT](https://tosc.iacr.org/index.php/ToSC/article/view/8466) in a format compatible with CryptoSMT, you can use the following commands:

```python
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import CRAFT as sb
sage: sa = SboxAnalyzer(sb)                                         
sage: cnf, milp = sa.minimized_diff_constraints(cryptosmt_compatible=True)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.02 seconds
Number of constraints: 54
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: p0 + p1 + p2

sage: print(cnf)

'(~a2 | p1) & (~b2 | p1) & (~p0 | p1) & (~p1 | p2) & (a1 | ~a2 | a3 | ~p0) & (a1 | ~a3 | b2 | p0) & (a2 | ~a3 | b2 | p0) & (~a1 | a3 | b2 | p0) & (a2 | b1 | ~b3 | p0) & (a2 | ~b1 | b3 | p0) & (~a0 | b2 | b3 | p0) & (a1 | a2 | a3 | ~b1 | ~b3) & (a0 | ~a2 | a3 | ~b2 | b3) & (~a1 | ~a2 | ~a3 | b2 | b3) & (~a1 | ~a3 | b0 | b1 | p0) & (a0 | a1 | ~b1 | ~b3 | p0) & (a1 | ~a3 | ~b1 | ~b3 | p0) & (~a1 | a3 | ~b1 | ~b3 | p0) & (~a0 | a3 | b1 | ~b3 | p0) & (~a1 | ~b0 | b1 | ~b3 | p0) & (a1 | ~a3 | ~b0 | b3 | p0) & (~a3 | ~b0 | ~b1 | b3 | p0) & (a0 | a1 | a2 | a3 | ~p2) & (b0 | b1 | b2 | b3 | ~p2) & (~a1 | ~a2 | ~b0 | b1 | b2 | ~b3) & (~a2 | a3 | b0 | ~b1 | b2 | ~p0) & (~a2 | ~a3 | b0 | b2 | ~b3 | ~p0) & (~a0 | a1 | ~a3 | ~b2 | b3 | ~p0) & (a0 | a3 | ~b0 | b1 | b3 | p0) & (a1 | a2 | a3 | b1 | b3 | ~p2) & (a0 | a3 | ~b0 | b1 | b2 | ~b3 | ~p0) & (~a0 | a1 | a2 | ~a3 | b0 | b3 | ~p0) & (a0 | a1 | ~b0 | ~b1 | b2 | b3 | ~p0) & (~a1 | ~a3 | b1 | b3 | ~p0 | ~p1 | ~p2) & (~a0 | ~a2 | b0 | ~b2 | p0 | ~p1 | ~p2) & (~a0 | ~a1 | a3 | ~b0 | b1 | ~b2 | ~p1 | ~p2) & (~a0 | a1 | ~a2 | ~b0 | ~b1 | b3 | ~p1 | ~p2) & (a0 | a2 | ~a3 | b1 | ~b2 | ~p0 | ~p1 | ~p2) & (a2 | b0 | ~b1 | ~b2 | ~b3 | ~p0 | ~p1 | ~p2) & (a0 | ~a1 | a2 | ~b2 | b3 | ~p0 | ~p1 | ~p2) & (a0 | ~a1 | ~a2 | ~a3 | ~b2 | p0 | ~p1 | ~p2) & (a0 | a1 | ~a2 | b1 | ~b2 | p0 | ~p1 | ~p2) & (~a0 | ~a1 | a2 | a3 | b0 | b1 | ~p0 | ~p1 | ~p2) & (~a2 | a3 | b1 | ~b2 | ~p0) & (a1 | ~a2 | ~b2 | b3 | ~p0) & (~a1 | ~a3 | ~b1 | ~b3 | ~p0) & (a2 | ~b0 | ~b1 | ~b3) & (~a0 | ~a3 | b0 | b1 | b2) & (a1 | ~a2 | b1 | ~b2 | ~p0) & (~a0 | ~a1 | b0 | b2 | b3) & (~a2 | a3 | ~b2 | b3 | ~p0) & (a0 | a1 | a2 | ~b0 | ~b3) & (a0 | a2 | a3 | ~b0 | ~b1) & (~a0 | ~a1 | ~a3 | b2)'
```

## Deriving and Encoding Deterministic Differential Propagation

The following example shows how to derive and encode the deterministic differential propagation of the Ascon S-box.

```python
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import Ascon as sb
sage: sa = SboxAnalyzer(sb)
sage: ddp = sa.encode_deterministic_differential_behavior()
sage: cp = sa.generate_cp_constraints(ddp)
Input:	a0||a1||a2||a3||a4; a0: msb
Output:	b0||b1||b2||b3||b4; b0: msb
sage: print(cp)

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
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import SKINNY_4 as sb
sage: sa = SboxAnalyzer(sb)
sage: cnf, milp = sa.minimized_linear_constraints()

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 29
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
```
Note that `sa.minimized_linear_constraints()` encode the squared LAT scaled by correlation. 

### Encoding the LAT of Ascon

```python
sage: from sage.crypto.sboxes import Ascon as sb
sage: sa = SboxAnalyzer(sb)
sage: cnf, milp = sa.minimized_linear_constraints()

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.04 seconds
Number of constraints: 96
Input:	a0||a1||a2||a3||a4; a0: msb
Output:	b0||b1||b2||b3||b4; b0: msb
Weight: 4.0000 p0 + 2.0000 p1
```

## Deriving and Encoding Deterministic Linear Propagation

The following example shows how to derive and encode the deterministic linear propagation through the inverse of the Ascon S-box.

```python
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import Ascon as sb
sage: sa = SboxAnalyzer(sb.inverse())
sage: dlp = sa.encode_deterministic_linear_behavior()
sage: cp = sa.generate_cp_constraints(dlp)
Input:	a0||a1||a2||a3||a4; a0: msb
Output:	b0||b1||b2||b3||b4; b0: msb
sage: print(cp)
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
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import PRESENT as sb
sage: sa = SboxAnalyzer(sb)
sage: mpt = sa.monomial_prediction_table(); mpt
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
sage: cnf, milp = sa.minimized_integral_constraints()

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 41
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
```

### Encoding the MPT of Ascon

```python
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import Ascon as sb
sage: sa = SboxAnalyzer(sb)
sage: cnf, milp = sa.minimized_integral_constraints()

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 97
Input:	a0||a1||a2||a3||a4; a0: msb
Output:	b0||b1||b2||b3||b4; b0: msb
```

## DLCT Encoding

Here, we show how to encode the DLCT of S-boxes. 

### Encoding the *-DLCT of KNOT

```python 
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import KNOT as sb
sage: sa = SboxAnalyzer(sb)
sage: cnf, milp = sa.minimized_differential_linear_constraints(subtable='star')

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 34
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
```

### Encoding the *-DLCT of Midori

```python
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import Midori_Sb0 as sb
sage: sa = SboxAnalyzer(sb)
sage: cnf, milp = sa.minimized_differential_linear_constraints(subtable='star')

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.00 seconds
Number of constraints: 43
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
```

### Verification of Hadipour et al.'s Theorem

The following code verifies the Hadipour et al.'s theorem, (Proposition 2 in [2024/255](https://ia.cr/2024/255)) for the S-box of Ascon.

```python
sage: from sboxanalyzer import *
sage: from sage.crypto.sboxes import Ascon as sb
sage: sa = SboxAnalyzer(sb)
sage: check = sa.check_hadipour_theorem(); check
The Hadipour et al.'s theorem is satisfied.
True
```

## Citation

If you use our tools in a project resulting in an academic publication, please acknowledge it by citing our paper:

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
```

## Papers

Sbox Analyzer has been used in the following papers:

- [Throwing Boomerangs into Feistel Structures: Application to CLEFIA, WARP, LBlock, LBlock-s and TWINE](https://ia.cr/2022/745)
- [Integral Cryptanalysis of WARP based on Monomial Prediction](https://ia.cr/2022/729)
- [Improved Search for Integral, Impossible-Differential and Zero-Correlation Attacks: Application to Ascon, ForkSKINNY, SKINNY, MANTIS, PRESENT and QARMAv2](https://ia.cr/2023/1701)
- [Revisiting Differential-Linear Attacks via a Boomerang Perspective With Application to AES, Ascon, CLEFIA, SKINNY, PRESENT, KNOT, TWINE, WARP, LBlock, Simeck, and SERPENT](https://ia.cr/2024/255)
- [Gleeok: A Family of Low-Latency PRFs andits Applications to Authenticated Encryption](https://tches.iacr.org/index.php/TCHES/article/view/11439/10944)
  
## Road Map

 - [x] Encoding DDT
 - [x] Encoding LAT
 - [x] Encoding [MPT](https://tosc.iacr.org/index.php/ToSC/article/view/9715)
 - [x] Adding [DLCT, UDLCT, LDLCT and Hadipour et al's theorem](https://eprint.iacr.org/2024/255)
 - [ ] Integrating the tool into the [SageMath](https://www.sagemath.org/)
 - [ ] Integrating the tool into the [CryptoSMT](https://github.com/kste/cryptosmt)

Any contributions or comments regarding the development of the tool are warmly welcome. 

## License

S-box Analyszer is released under the [MIT license](LICENSE).

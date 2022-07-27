# S-box Analyzer

S-box Analyzer is a tool for deriving the minimized MILP and SMT/SAT constraints to encode the Differential Distribution Table (DDT), Linear Approximation Table (LAT), or Monomial Prediction Table (MPT) of S-boxes.

---
![logo](./images/sboxanalyzer.svg)

- [S-box Analyzer](#s-box-analyzer)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
  - [Usage](#usage)
  - [More Examples](#more-examples)
    - [Encoding the DDT of SKINNY-64](#encoding-the-ddt-of-skinny-64)
    - [Encoding the DDT of Ascon](#encoding-the-ddt-of-ascon)
    - [Encoding the DDT of PRESENT](#encoding-the-ddt-of-present)
    - [Encoding the DDT of SKINNY-128](#encoding-the-ddt-of-skinny-128)
    - [Encoding the DDT of AES](#encoding-the-ddt-of-aes)
    - [Encoding the DDT for CryptoSMT](#encoding-the-ddt-for-cryptosmt)
  - [Paper](#paper)
  - [Road Map](#road-map)
  - [License](#license)


---

## Dependencies

S-box Analyzer has been implemented as a SageMath module and employs ESPRESSO for logic minimization. Thus, it requires the following software:

- [ESPRESSO](https://ptolemy.berkeley.edu/projects/embedded/pubs/downloads/espresso/index.htm)
- [SageMath](https://www.sagemath.org/)

## Installation

- ESPRESSO

  To build ESPRESSO, use the following commands:

  ```bash
  cd espresso
  make
  ```

  To see more details about the ESPRESSO logic minimizer, see [here](https://ptolemy.berkeley.edu/projects/embedded/pubs/downloads/espresso/index.htm). An updated version of ESPRESSO which is compatible with new compilers is available [here](https://github.com/classabbyamp/espresso-logic).

- SageMath

  To see the installation recipe of SageMath, see [here](https://doc.sagemath.org/html/en/installation/index.html).

## Usage

Open SageMath in the same directory as the S-box analyzer. Import `sboxanalyzer` into the SageMath and then simply use it. The following example shows you how to use the S-box analyzer to encode the DDT of an S-box:

```python
sage: from sboxanalyzer import *  
sage: from sage.crypto.sboxes import PRINTcipher as sb
sage: sa = SboxAnalyzer(sb)                               
sage: cnf, milp = sa.minimized_diff_constraints()

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 17
Input:	a0||a1||a2; a0: msb
Output:	b0||b1||b2; b0: msb
Weight: 2.0000 p0

sage: print(milp)
['- a1 + p0 >= 0',
'- b0 + p0 >= 0', 
'- b2 + p0 >= 0', 
'a0 + a1 - a2 + b2 >= 0', 
'a1 + b0 - b1 + b2 >= 0', 
'- a0 + b0 + b1 + b2 >= 0', 
'a0 + a1 + a2 - p0 >= 0', 
'a1 + a2 + b0 - p0 >= 0', 
'a0 + a2 + b1 - p0 >= 0', 
'a2 + b0 + b1 - p0 >= 0', 
'a0 + b1 + b2 - p0 >= 0', 
'- a1 - a2 + b0 - b1 - b2 >= -3', 
'a0 - a1 - a2 - b1 - b2 >= -3', 
'- a0 - a2 - b0 + b1 - b2 >= -3', 
'- a0 - a1 - b0 - b1 + b2 >= -3', 
'- a0 + a1 - a2 - b0 - b2 >= -3', 
'- a0 - a1 + a2 - b0 - b1 >= -3']

sage: print(cnf)
(~a1 | p0) & (~b0 | p0) & (~b2 | p0) & (a0 | a1 | ~a2 | b2) & 
(a1 | b0 | ~b1 | b2) & (~a0 | b0 | b1 | b2) & (a0 | a1 | a2 | ~p0) & 
(a1 | a2 | b0 | ~p0) & (a0 | a2 | b1 | ~p0) & (a2 | b0 | b1 | ~p0) & (a0 | b1 | b2 | ~p0) & 
(~a1 | ~a2 | b0 | ~b1 | ~b2) & (a0 | ~a1 | ~a2 | ~b1 | ~b2) & (~a0 | ~a2 | ~b0 | b1 | ~b2) & 
(~a0 | ~a1 | ~b0 | ~b1 | b2) & (~a0 | a1 | ~a2 | ~b0 | ~b2) & (~a0 | ~a1 | a2 | ~b0 | ~b1)
```

Interpretation of the output:

- `Input:	a0||a1||a2; a0: msb`: $a_{0}||a_{1}||a_{2}$ encode the input difference vector $a$ where $a_{0}$ is the most significant bit of $a$.
- `Output:	b0||b1||b2; b0: msb`: $b_{0}||b_{1}||b_{2}$ encode the output difference vector $b$ where $b_{0}$ is the most significant bit of $b$.
- `Weight: 2.0000 p0`: $2\cdot p0$ is the weight of differential transition $a \rightarrow b$, i.e.,  $\Pr \{a \rightarrow b\} = 2^{-2\cdot p_0}$, where $p_{0}$ is a binary decision variable encoding the probability of valid transitions.

## More Examples

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

**Our tool supports 7 different modes, i.e., `[mode=1,...,mode=7]`, for minimization to make a trade off between the time of simplification and the optimality of the solution.** For example using the following command we can minimize the number of constraints further:

```python
sage: cnf, milp = sa.minimized_diff_constraints(mode=5)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.22 seconds
Number of constraints: 37
Input:	a0||a1||a2||a3; a0: msb
Output:	b0||b1||b2||b3; b0: msb
Weight: 3.0000 p0 + 2.0000 p1
```

**The default mode is `6` since it is the best choice in terms of both time and optimality.**

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

To encode the DDT of larges S-boxes (8-bit S-boxes), we usually divide the DDT into several sub-DDTs and encode the sub-DDTs seperately. Lastly, we put the constraints for all sub-DDTs together to encode the whole DDT. The following code shows how to encode the DDT of SKINNY-128.

***Encode 2-DDT***

```python
sage: from sage.crypto.sboxes import SKINNY_8 as sb                            
sage: sa = SboxAnalyzer(sb)                                                    
sage: sa.diff_spectrum                                                         
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

You can achieve more optimal solutions by `mode=2` or `mode=5`. However, the simplification time will be higher. I do not recommend it since the solution derived by `mode=2` is nearly optimum and sufficient to get a remarkable speed up in automatic differential analysis based on MILP or SAT/SMT. You can encode other sub-DDTs of SKINNY-8 in a similar way.

### Encoding the DDT of AES

```python
sage: from sage.crypto.sboxes import AES as sb                                 
sage: sa = SboxAnalyzer(sb)                                                    
sage: sa.diff_spectrum                                                         
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

As you can see our results concerning the encoding of AES DDT is much better than the results reported in [this paper](https://tosc.iacr.org/index.php/ToSC/article/view/805).

### Encoding the DDT for CryptoSMT

By setting the `cryptosmt_compatible` argument to `True`, you can generate an SMT formula compatible with CryptoSMT. For example, to encode the DDT of [CRAFT](https://tosc.iacr.org/index.php/ToSC/article/view/8466) in a format compatible with CryptoSMT, you can use the following commands:

```python
sage: from sboxanalyzer import *                                                
sage: from sage.crypto.sboxes import PRINTcipher as sb                          
sage: sa = SboxAnalyzer(sb)                                                     
sage: cnf, milp = sa.minimized_diff_constraints(cryptosmt_compatible=True)

Simplifying the MILP/SAT constraints ...
Time used to simplify the constraints: 0.01 seconds
Number of constraints: 18
Input:	a0||a1||a2; a0: msb
Output:	b0||b1||b2; b0: msb
Weight: p0 + p1

sage: print(cnf)

'(p0 | ~p1) & (~a1 | p1) & (~b0 | p1) & (~b2 | p1) & (a1 | a2 | b0 | ~b2) & (a0 | a2 | b1 | ~b2) & (a2 | b0 | b1 | ~b2) & (a0 | a1 | ~a2 | b2) & (a1 | b0 | ~b1 | b2) & (~a0 | b0 | b1 | b2) & (a0 | b1 | b2 | ~p0) & (a0 | a1 | a2 | ~p1) & (~a1 | ~a2 | b0 | ~b1 | ~b2) & (a0 | ~a1 | ~a2 | ~b1 | ~b2) & (~a0 | ~a2 | ~b0 | b1 | ~b2) & (~a0 | ~a1 | ~b0 | ~b1 | b2) & (~a0 | a1 | ~a2 | ~b0 | ~b2) & (~a0 | ~a1 | a2 | ~b0 | ~b1)'

```

## Paper

If you use our tools in a project resulting in an academic publication, please acknowledge it by citing our paper:

```bib
@misc{cryptoeprint:2022/745,
      author = {Hosein Hadipour and Marcel Nageler and Maria Eichlseder},
      title = {Throwing Boomerangs into Feistel Structures: Application to CLEFIA, WARP, LBlock, LBlock-s and TWINE},
      howpublished = {Cryptology ePrint Archive, Paper 2022/745},
      year = {2022},
      note = {\url{https://eprint.iacr.org/2022/745}},
      url = {https://eprint.iacr.org/2022/745}
}
```

## Road Map

 - [x] Encoding DDT
 - [ ] Encoding LAT
 - [ ] Encoding MPT
 - [ ] Integrating the tool into the [SageMath](https://www.sagemath.org/)
 - [ ] Integrating the tool into the [CryptoSMT](https://github.com/kste/cryptosmt)

The LAT and the MPT encoders are not implemented yet. However, they will follow exactly the same template as the DDT encoder and can be easily implemented. I will include them as soon as I get a chance. Any contributions or comments regarding the development of the tool are warmly welcome. 

## License

S-box Analyszer is released under the [GPlv3 license](LICENSE.txt).

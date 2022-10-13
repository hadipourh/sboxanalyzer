# S-box Analyzer

S-box Analyzer is a tool for deriving the minimized MILP and SMT/SAT constraints to encode the Differential Distribution Table (DDT), Linear Approximation Table (LAT), and [Monomial Prediction Table (MPT)](https://tosc.iacr.org/index.php/ToSC/article/view/9715) of S-boxes.

---
![logo](./images/sboxanalyzer.svg)

- [S-box Analyzer](#s-box-analyzer)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Examples](#examples)
    - [Encoding the DDT of SKINNY-64](#encoding-the-ddt-of-skinny-64)
    - [Encoding the DDT of Ascon](#encoding-the-ddt-of-ascon)
    - [Encoding the DDT of PRESENT](#encoding-the-ddt-of-present)
    - [Encoding the DDT of SKINNY-128](#encoding-the-ddt-of-skinny-128)
    - [Encoding the DDT of AES](#encoding-the-ddt-of-aes)
    - [Encoding the *-DDT and *-LAT](#encoding-the--ddt-and--lat)
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

## Examples

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
You can encode other sub-DDTs of SKINNY-128's S-box in a similar way. Moreover, you may achieve a more optimal solution using the different modes such as `mode=2` or `mode=5`. However, the simplification time might be higher. The default mode is `mode=6` since it generates the nearly optimum result and is sufficient to get a remarkable speed up in automatic differential analysis based on MILP or SAT/SMT.

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
sage: from sage.crypto.sboxes import CRAFT as sb                                                            
sage: sa = SboxAnalyze(sb)                                                                                                                                   
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

## Paper

If you use our tools in a project resulting in an academic publication, please acknowledge it by citing our paper:

```bib
@article{DBLP:journals/tosc/HadipourNE22,
  author    = {Hosein Hadipour and
               Marcel Nageler and
               Maria Eichlseder},
  title     = {Throwing Boomerangs into Feistel Structures Application to CLEFIA,
               WARP, LBlock, LBlock-s and {TWINE}},
  journal   = {{IACR} Trans. Symmetric Cryptol.},
  volume    = {2022},
  number    = {3},
  pages     = {271--302},
  year      = {2022},
  doi       = {10.46586/tosc.v2022.i3.271-302}
}
```

## Road Map

 - [x] Encoding DDT
 - [ ] Encoding LAT
 - [ ] Encoding [MPT](https://tosc.iacr.org/index.php/ToSC/article/view/9715)
 - [ ] Integrating the tool into the [SageMath](https://www.sagemath.org/)
 - [ ] Integrating the tool into the [CryptoSMT](https://github.com/kste/cryptosmt)

The LAT and the MPT encoders are not implemented yet. However, they will follow the same template as the DDT encoder and can be easily implemented. I will include them as soon as I get a chance. Any contributions or comments regarding the development of the tool are warmly welcome. 

## License

S-box Analyszer is released under the [GPlv3 license](LICENSE.txt).

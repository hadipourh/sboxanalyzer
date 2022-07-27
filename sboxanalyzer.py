# -*- coding: utf-8 -*-
#!/usr/local/bin/sage
"""
SA: S-box Analyzer

AUTHORS:

- Hosein Hadipour (2022-05-26)

EXAMPLES::

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
"""

#*****************************************************************************
#       Copyright (C) 2022 Hosein Hadipour <hsn.hadipour@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import subprocess
import os
import time
from sage.all import *
from sage.crypto.sboxes import SBox

# ESPRESO_BIN_PATH = os.path.join(os.environ['SAGE_ROOT'], 'local/bin/espresso')
ESPRESO_BIN_PATH = os.path.join(os.getcwd(), 'bin', 'espresso')


class SboxAnalyzer(SBox):
    r"""
    This module encodes the DDT, LAT and MPT [HE22]_ of a given S-box with MILP/SAT constraints and 
    then simplifies the extracted constraints using logic minimization tools

    EXAMPLES:
    
    We consider the S-box of block cipher PRINTcipher [KLPR2010]_::

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
    
    AUTHORS:

    - Hosein Hadipour (2022-05-26)

    REFERENCES:

    - [KLPR2010]_
    - [HE22]_

    """

    sbox_counter = 0

    def __init__(self, lookuptable):
        """
        Initialize the lookup table of S-box

        :param lookuptable list: list of integers or hexadecimal numbers specifying the S-box mapping
        """

        super().__init__(lookuptable)
        SboxAnalyzer.sbox_counter += 1
        if not os.path.exists(os.path.join(os.getcwd(), 'tmp')):
            os.makedirs(os.path.join(os.getcwd(), 'tmp'))
        self.truth_table_filename = os.path.join(os.getcwd(), 'tmp', 'tt_' + str(SboxAnalyzer.sbox_counter) + '.txt')
        self.simplified_truth_table_filename = os.path.join(os.getcwd(), 'tmp', 'stt_' + str(SboxAnalyzer.sbox_counter) + '.txt')        
        self.ddt = self.difference_distribution_table()
        self.diff_spectrum = set([self.ddt[i][j] for i in range(2**self.m) for j in range(2**self.n)]) - {0, 2**self.m}
        self.diff_spectrum = sorted(list(self.diff_spectrum))
        self.npp = len(self.diff_spectrum)        
        self.diff_weights = [abs(float(log(d/(2**self.m), 2))) for d in self.diff_spectrum]        
        self.nww = int(max(self.diff_weights))
        self.variables_mapping = "Input:\t{}; a0: msb".format("||".join([f"a{i}" for i in range(self.m)]))
        self.variables_mapping += "\nOutput:\t{}; b0: msb".format("||".join([f"b{i}" for i in range(self.n)]))
        self.cryptosmt_compatible = False
        self.ddt_subtable = None
        # self.lat = self.linear_approximation_table(scale="correlation")
        # TODO
        # self.mpt = self.monomial_prediction_table()
        # TODO

    def compute_star_ddt(self, reverse=1):
        """
        Generate the start DDT (or 0/1 DDT)
        Star DDT is a 2^m*2^n binary array describing the possibility of differential transitions through the S-box
        """
        
        self.star_ddt = [[0 for i in range(2**self.m)] for j in range(2**self.n)]
        for dx in range(2**self.m):
            for dy in range(2**self.n):
                if self.ddt[dx][dy] != 0:
                    self.star_ddt[dx][dy] = reverse ^ 1
                else:
                    self.star_ddt[dx][dy] = reverse
        return self.star_ddt

    def _star_ddt_to_boolean_function(self, reverse=1):
        """
        Convert the star-DDT into a Boolean function
        """
        
        self.compute_star_ddt(reverse=reverse)
        boolean_func = dict()
        for dx in range(2**self.m):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.m))))
            for dy in range(2**self.n):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.m))))
                key = x + y
                boolean_func[key] = self.star_ddt[dx][dy]
        boolean_func[tuple([0]*self.m + [0]*self.n)] = reverse ^ 1
        return boolean_func

    def _ddt_to_boolean_function(self, reverse=1):
        """
        Convert the DDT into a Boolean function
        To encode the probabilities, we define as many new binary variables as the 
        number of different (non-zero and non-one) probabilities in DDT and 
        then encode the whole of DDT as a Boolean function. 
        Let dx, and dy denote the input and output differences, respectively. 
        Assuming that the DDT of S-box includes for example three (non-zero and non-one) elements: 
        e0, e1, and e2, we define three new binary variables p0, p1, and p2 and then 
        encode the DDT as a Boolean function f such that f(Bin(dx) || Bin(dy) || p0 || p1 || p2) = 1
        if and only if:
               DDT[dx][dy] = e0 and (p0, p1, p2) = (1, 0, 0),
            or DDT[dx][dy] = e1 and (p0, p1, p2) = (0, 1, 0),
            or DDT[dx][dy] = e2 and (p0, p1, p2) = (0, 0, 1),
            or dx = dy = 0      and (p0, p1, p2) = (0, 0, 0),
        otherwise, f(Bin(dx) || Bin(dy) || p0 || p1 || p2) = 0.

        If reverse = True, then we compute the complement of the derived Boolean function.
        """
        
        boolean_function = dict()
        complexity = self.m + self.n        
        for dx in range(2**self.m):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.m))))
            for dy in range(2**self.n):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.n))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                if self.ddt[dx][dy] in self.diff_spectrum:
                    p = tuple([int(i == self.ddt[dx][dy]) for i in self.diff_spectrum])
                    key = x + y + p
                    boolean_function[key] = 1
        boolean_function[tuple([0]*self.m + [0]*self.n + [0]*self.npp)] = 1
        if reverse == 1:            
            complexity = self.m + self.n + self.npp            
            for dx in range(2**complexity):
                x = tuple(map(int, list(bin(dx)[2:].zfill(complexity))))
                boolean_function[x] = boolean_function.get(x, 0) ^ 1
        return boolean_function

    def _ddt_to_cryptosmt_compatible_boolean_function(self, reverse=1):
        """
        Encode the DDT as a Boolean function in SAT/SMT-compatible form 
        All transition's weights should be integers in this case.
        To encode the DDT/LAT of S-box in this method we define the auxiliary binary variables p0, p1, ..., pn 
        such that p0 + ... + pn is equal to the weight of the corresponding transition. 
        For example, we encode the DDT of a 4-uniform S-box with an 11-bit Boolean function f(x||y||p), 
        where x and y are the 4-bit input and output differences and p = (p0, p1, p2) such that:
        f(x||y||p) = 0 if DDT[x, y] = 0
        f(x||y||p) = 1 if DDT[x, y] = 2^(-3) and p = (1, 1, 1)
        f(x||y||p) = 0 if DDT[x, y] = 2^(-3) and p != (1, 1, 1)
        f(x||y||p) = 1 if DDT[x, y] = 2^(-2) and p = (0, 1, 1)
        f(x||y||p) = 0 if DDT[x, y] = 2^(-2) and p != (0, 1, 1)
        f(x||y||p) = 1 if DDT[x, y] = 1 and p = (0, 0, 0)
        f(x||y||p) = 0 if DDT[x, y] = 1 and p != (0, 0, 0)
        
        This method is suitable for encoding the DDT of an S-box in CryptoSMT.
        """
        
        if not all([i == int(i) for i in self.diff_weights]):
            raise ValueError("All transition's weights should be integers")

        boolean_function = dict()
        complexity = self.m + self.n
        for dx in range(2**self.m):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.m))))
            for dy in range(2**self.n):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.n))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                if self.ddt[dx][dy] != 0:
                    w = int(abs(float(log(self.ddt[dx][dy]/(2**self.m), 2))))
                    p = tuple([0]*(self.nww - w) + [1]*w)
                    key = x + y + p
                    boolean_function[key] = 1        
        if reverse == 1:            
            complexity = self.m + self.n + self.npp
            for dx in range(2**complexity):
                x = tuple(map(int, list(bin(dx)[2:].zfill(complexity))))
                boolean_function[x] = boolean_function.get(x, 0) ^ 1
        return boolean_function

    def _pddt_to_booleanfunction(self, p, reverse=1):
        """
        Convert the p-DDT into a Boolean function
        """

        boolean_function = dict()
        for dx in range(2**self.m):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.m))))
            for dy in range(2**self.n):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.n))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                key = x + y
                if self.ddt[dx][dy] == p:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function

    def _write_truth_table(self, filename, boolean_function):
        """
        Write the Boolean function encoding the DDT of S-box into 
        a text file according to the input format of
        the ESPRESSO

        :param filename str: an string specifying the filename
        """

        if self.ddt_subtable != None:
            file_contents = ".i %d\n" % (self.m + self.n)
            file_contents += ".o 1\n"
            file_contents += ".ilb " + " ".join("a%d" % i for i in range(self.m)) + " " +\
                                   " ".join("b%d" % i for i in range(self.n)) + "\n"
        else:
            if self.cryptosmt_compatible:
                num_of_p_vars = self.nww
            else:
                num_of_p_vars = self.npp
            file_contents = ".i %d\n" % (self.m + self.n + num_of_p_vars)
            file_contents += ".o 1\n"
            file_contents += ".ilb " + " ".join("a%d" % i for i in range(self.m)) + " " +\
                                    " ".join("b%d" % i for i in range(self.n)) + " " +\
                                    " ".join("p%d" % i for i in range(num_of_p_vars)) + "\n"            
        file_contents += ".ob F\n"
        keys = list(boolean_function.keys())        
        for key in keys:
            line = "".join(map(str, key))
            line += "%s\n" % str(boolean_function[key])
            file_contents += line
        with open(filename, "w") as fileobj:
            fileobj.write(file_contents)
            fileobj.write(".e\n")

    def minimized_diff_constraints(self, mode=6, subtable=None, cryptosmt_compatible=False):
        """
        Given a Boolean function, this method writes its truth table into a file
        following the ESPRESSO input format, and calls ESPRESSO to derive its minimized
        representation. Next, it parses the output of ESPRESSO and translates the derived
        representation to the language of MILP and SAT solvers.

        :param mode list: a set of flags specifying the configuration of ESPRESSO program
        :param booleanfunction dict: a Python dictionary representing the truth table of the given Boolean function
        :rtype: string
        :return: the minimized MILP/SAT representation of the given Boolean function
        """
        self.cryptosmt_compatible = cryptosmt_compatible
        self.ddt_subtable = subtable
                                                                        # If reverse = 1 choose ON-SET and if reverse = 0 choose OFF-SET
        espresso_options = [[],                                         # 0 ON-SET :   Derived constraints exclude point p such that f(p) = 1
                            ["-Dexact", "-estrong", "-s", "-t", "-or"], # 1 OFF-SET:*  Derived constraints exclude point p such that f(p) = 0
                            ["-Dexact", "-estrong", "-s", "-t", "-of"], # 2 ON-SET :   Derived constraints exclude point p such that f(p) = 1 *
                            ["-Dexact", "-efast", "-s", "-t", "-or"],   # 3 OFF-SET:*  Derived constraints exclude point p such that f(p) = 0
                            ["-Dexact", "-efast", "-s", "-t", "-of"]]   # 4 ON-SET :   Derived constraints exclude point p such that f(p) = 1
        # ATTENTION!
        # -epos: complements the boolean function
        # Therefore, instead of reverse=1 we can use -epos
        # If so, we should make sure that reverse=0 in this case
        espresso_options += [["-Dexact", "-estrong", "-epos", "-s", "-t", "-of"],   # 5 ON-SET :*  Derived constraints exclude point p such that f(p) = 1]
                            ["-Dmany", "-estrong", "-epos", "-s", "-t", "-of"],    # 6 ON-SET :*  Derived constraints exclude point p such that f(p) = 1]
                            ["-Dmany", "-efast", "-epos", "-s", "-t", "-of"]]      # 7 ON-SET :*  Derived constraints exclude point p such that f(p) = 1]

        # Best options (not always), in terms of optimality: 5, then 6, then 7, then 1 (6 is good and fast enough in cases 5 is too slow)
                
        assert(mode in [0, 1, 2, 3, 4, 5, 6, 7])
        assert(subtable in ["star", None] + list(self.diff_spectrum))
        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1
        
        self.diff_objective = ""
        if subtable == "star":
            boolean_function = self._star_ddt_to_boolean_function(reverse=reverse)
        elif subtable in self.diff_spectrum:
            boolean_function = self._pddt_to_booleanfunction(p=subtable, reverse=reverse)
        elif cryptosmt_compatible:
            boolean_function = self._ddt_to_cryptosmt_compatible_boolean_function(reverse=reverse)
            self.diff_objective = ["p{:d}".format(i) for i in range(self.nww)]
            self.diff_objective = "\nWeight: {}".format(" + ".join(self.diff_objective))
        else:
            boolean_function = self._ddt_to_boolean_function(reverse=reverse)
            self.diff_objective = ["{:0.04f} p{:d}".format(self.diff_weights[i], i) for i in range(self.npp)]
            self.diff_objective = "\nWeight: {}".format(" + ".join(self.diff_objective)     )

        self._write_truth_table(filename=self.truth_table_filename, boolean_function=boolean_function)
        starting_time = time.time()
        print("Simplifying the MILP/SAT constraints ...")
        with open(self.simplified_truth_table_filename, 'w') as fileobj:
            subprocess.call([ESPRESO_BIN_PATH, *espresso_options[mode], self.truth_table_filename], stdout=fileobj)
        elapsed_time = time.time() - starting_time
        with open(self.simplified_truth_table_filename, 'r') as fileobj:
            espresso_output = fileobj.readlines()              
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # Parse the output of ESPRESSO
        if self.ddt_subtable == "star" or self.ddt_subtable in self.diff_spectrum:
            alphabet = ['a%d' % i for i in range(self.m)] + \
                       ['b%d' % i for i in range(self.n)]
        else:
            alphabet = ['a%d' % i for i in range(self.m)] + \
                       ['b%d' % i for i in range(self.n)]
            if self.cryptosmt_compatible:
                alphabet += ['p%d' % i for i in range(self.nww)]
            else:
                alphabet += ['p%d' % i for i in range(self.npp)]
        bd = len(alphabet)
        milp_constraints = []
        sat_clauses = []
        starting_point = 0
        end_point = 0
        for i in range(len(espresso_output)):
            if ".p" in espresso_output[i]:
                starting_point = i + 1
                number_of_constraints = espresso_output[i].split(" ")[1]
            if ".e" in espresso_output[i]:
                end_point = i
        for l in espresso_output[starting_point:end_point]:
            line = l[0:bd]
            orclause = []
            lp = []
            lp_rhs = 0
            for i in range(bd):
                if line[i] == '0':
                    orclause.append(alphabet[i])
                    lp.append(" + {}".format(alphabet[i]))
                elif line[i] == '1':
                    orclause.append("~{}".format(alphabet[i]))
                    lp.append(" - {}".format(alphabet[i]))
                    lp_rhs += 1
            sat_clauses.append("({})".format(' | '.join(orclause)))
            lp_rhs = -(lp_rhs - 1)
            lp_constraint = "".join(lp) + " >= {}".format(lp_rhs)
            if lp_constraint[0:3] == " + ":
                lp_constraint = lp_constraint[3:]
            else:
                lp_constraint = lp_constraint[1:]
            milp_constraints.append("{}".format(lp_constraint))
        sat_clauses = ' & '.join(sat_clauses)        
        print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {}".format(number_of_constraints.rstrip()))
        print("{}".format(self.variables_mapping + self.diff_objective))
        return sat_clauses, milp_constraints

if __name__ == "__main__":
    from sage.crypto.sboxes import SKINNY_8 as sb
    sa = SboxAnalyzer(sb)
    print("Spectrum of DDT entries: {}".format(sa.diff_spectrum))    
    cnf, milp = sa.minimized_diff_constraints(subtable=2, mode=5)
    with open("milp.lp", "w") as fileobj:
        fileobj.write("\nsubject to\n")
        fileobj.write("\n".join(milp))
        fileobj.write("\nend")
    with open("cnf.txt", "w") as fileobj:
        fileobj.write(cnf)

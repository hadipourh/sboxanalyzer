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

# MIT License
# Copyright (c) 2022 Hosein Hadipour

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#*****************************************************************************

import subprocess
import os
import time
from sage.all import *
from sage.crypto.sboxes import SBox
import itertools

# ESPRESO_BIN_PATH = os.path.join(os.environ['SAGE_ROOT'], 'local/bin/espresso')
ESPRESO_BIN_PATH = os.path.join(os.getcwd(), 'espresso', 'build', 'espresso')

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
        self.variables_mapping = "Input:\t{}; a0: msb".format("||".join([f"a{i}" for i in range(self.input_size())]))
        self.variables_mapping += "\nOutput:\t{}; b0: msb".format("||".join([f"b{i}" for i in range(self.output_size())]))

        # A flag to check if the data required for differential analysis are present in memory
        self._data_required_for_differential_analysis = None
        # A flag to check if the data required for linear analysis are present in memory    
        self._data_required_for_linear_analysis = None
        # A flag to check if the data required for integral analysis are presenet in memory
        self._data_required_for_integral_analysis = None
        # A flag to check if the data required for differential-linear analysis are presenet in memory
        self._data_required_for_difflin_analysis = None

        # define a dictionalry to encode deterministic behavior
        self.unknown, self.zero, self.one = -1, 0, 1
        self.deterministic_mask = {self.zero: {0}, self.one:{1}, self.unknown:{0, 1}}
    
    @staticmethod
    def print_table(table):
        """
        Prints the table in a nice format
        """
        
        column_widths = [max(len(str(row[i])) for row in table) for i in range(len(table[0]))]
        for row in table:
            formatted_row = [str(value).rjust(width) for value, width in zip(row, column_widths)]
            print(" ".join(formatted_row))

    def monomial_prediction_table(self):
        """
        Compute the Monomial Prediction Table (MPT) based on [HE22]:
        https://tosc.iacr.org/index.php/ToSC/article/view/9715
        """
        
        """
        SageMath's naming convention for Sbox:
        x0: MSB

        SageMath's naming convention for the ANF of Boolean function:
        x0: LSB
        Example:

        x1, x0 | y
        -----------
        0 , 0  | 1
        0 , 1  | 0
        1 , 0  | 1
        1 , 1  | 0

        y = 1 + x0
        ----------
        from sage.crypto.boolean_function import BooleanFunction
        f = BooleanFunction([1, 0, 1, 0])
        anf = f.algebraic_normal_form()
        Sage output: x0 + 1
        """
        
        mpt = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        BPR = BooleanPolynomialRing(n=self.input_size(), names="x")        
        
        output_components = [0]*self.output_size()
        for i in range(self.output_size()):
            shift_value = 1 << (self.output_size() - i - 1)
            output_mask = self.to_bits(x=shift_value, n=self.output_size())
            output_components[i] = self.component_function(output_mask).algebraic_normal_form()
            
        input_vars = list(BPR.gens())
        input_vars.reverse()        
        for u in range(2**self.input_size()):
            input_mask = self.to_bits(x=u, n=self.input_size())
            monomial = product([input_vars[i]**input_mask[i] for i in range(self.input_size())])           
            for v in range(2**self.output_size()):
                output_mask = self.to_bits(x=v, n=self.output_size())
                f = product([output_components[i]**output_mask[i] for i in range(self.output_size())])
                if monomial in f.monomials():
                    mpt[u][v] = 1
        return mpt

    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    #  ___         _                __                     _           _____  ____   ____   ____   _____  ____  ____    ___  
    # |_ _| _ __  | |_  ___  _ __  / _|  __ _   ___  ___  | |_  ___   | ____|/ ___| |  _ \ |  _ \ | ____|/ ___|/ ___|  / _ \ 
    #  | | | '_ \ | __|/ _ \| '__|| |_  / _` | / __|/ _ \ | __|/ _ \  |  _|  \___ \ | |_) || |_) ||  _|  \___ \\___ \ | | | |
    #  | | | | | || |_|  __/| |   |  _|| (_| || (__|  __/ | |_| (_) | | |___  ___) ||  __/ |  _ < | |___  ___) |___) || |_| |
    # |___||_| |_| \__|\___||_|   |_|   \__,_| \___|\___|  \__|\___/  |_____||____/ |_|    |_| \_\|_____||____/|____/  \___/                                                                                                                                                                 
    # Interface to ESPRESSO
    ###############################################################################################################

    def _write_truth_table(self, filename, boolean_function, input_output_variables=None):
        """
        Write the Boolean function encoding the DDT of S-box into 
        a text file according to the input format of
        the ESPRESSO

        :param filename str: an string specifying the filename
        """

        if input_output_variables is None:
            no_of_input_vars = len(next(iter(boolean_function.keys())))
            file_contents = f".i {no_of_input_vars}\n"
            file_contents += ".o 1\n"
            file_contents += ".ilb {0}\n".format(" ".join(f"x{i}" for i in range(no_of_input_vars)))
        else:         
            file_contents = f".i {len(input_output_variables)}\n"
            file_contents += ".o 1\n"
            file_contents += ".ilb {0}\n".format(" ".join(input_output_variables))         
        file_contents += ".ob F\n"        
        for key, value in boolean_function.items():
            file_contents += "{0}{1}\n".format("".join(map(str, key)), str(value))            
        with open(filename, "w") as fileobj:
            fileobj.write(file_contents)
            fileobj.write(".e\n")
    
    def simplify_by_espresso(self, input_file, output_file, mode):
        """
        Simplify the CNF formula using the ESPRESSO
        """
        valid_values_for_mode = list(range(8))
        if mode not in valid_values_for_mode:
            raise ValueError("Invalid value for mode! mode must be in [0, 1, 2, 3, 4, 5, 6, 7].") 

                                                                        # If reverse = 1 choose ON-SET and if reverse = 0 choose OFF-SET
        self.espresso_options = [[],                                    # 0 ON-SET :   Derived constraints exclude point p such that f(p) = 1
                            ["-Dexact", "-estrong", "-s", "-t", "-or"], # 1 OFF-SET:*  Derived constraints exclude point p such that f(p) = 0
                            ["-Dexact", "-estrong", "-s", "-t", "-of"], # 2 ON-SET :   Derived constraints exclude point p such that f(p) = 1 *
                            ["-Dexact", "-efast", "-s", "-t", "-or"],   # 3 OFF-SET:*  Derived constraints exclude point p such that f(p) = 0
                            ["-Dexact", "-efast", "-s", "-t", "-of"]]   # 4 ON-SET :   Derived constraints exclude point p such that f(p) = 1
        # ATTENTION!
        # -epos: complements the boolean function
        # Therefore, instead of reverse=1 we can use -epos
        # If so, we should make sure that reverse=0 in this case
        self.espresso_options += [["-Dexact", "-estrong", "-epos", "-s", "-t", "-of"],   # 5 ON-SET :*  Derived constraints exclude point p such that f(p) = 1]
                            ["-Dmany", "-estrong", "-epos", "-s", "-t", "-of"],    # 6 ON-SET :*  Derived constraints exclude point p such that f(p) = 1]
                            ["-Dmany", "-efast", "-epos", "-s", "-t", "-of"]]      # 7 ON-SET :*  Derived constraints exclude point p such that f(p) = 1]

        # Best options (not always), in terms of optimality: 5, then 6, then 7, then 1 (6 is good and fast enough in cases 5 is too slow)
        with open(output_file, 'w') as fileobj:
            subprocess.call([ESPRESO_BIN_PATH, *self.espresso_options[mode], input_file], stdout=fileobj)
        
    def _parse_the_output_of_espresso(self, filename, alphabet):
        """
        Parse the output of ESPRESSO
        """
        
        with open(filename, 'r') as fileobj:
            espresso_output = fileobj.readlines()        
        milp_constraints = []
        sat_clauses = []
        starting_point = 0
        end_point = 0
        for i in range(len(espresso_output)):
            if ".p" in espresso_output[i]:
                starting_point = i + 1
                # number_of_constraints = espresso_output[i].split(" ")[1]
            if ".e" in espresso_output[i]:
                end_point = i
        for l in espresso_output[starting_point:end_point]:
            line = l[0:len(alphabet)]
            orclause = []
            lp = []
            lp_rhs = 0
            for i in range(len(alphabet)):
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
        return (sat_clauses, milp_constraints)
        
    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    #  ____   _   __   __                          _    _         _      _                   _              _        ____               _   
    # |  _ \ (_) / _| / _|  ___  _ __  ___  _ __  | |_ (_)  __ _ | |    / \    _ __    __ _ | | _   _  ___ (_) ___  |  _ \  __ _  _ __ | |_ 
    # | | | || || |_ | |_  / _ \| '__|/ _ \| '_ \ | __|| | / _` || |   / _ \  | '_ \  / _` || || | | |/ __|| |/ __| | |_) |/ _` || '__|| __|
    # | |_| || ||  _||  _||  __/| |  |  __/| | | || |_ | || (_| || |  / ___ \ | | | || (_| || || |_| |\__ \| |\__ \ |  __/| (_| || |   | |_ 
    # |____/ |_||_|  |_|   \___||_|   \___||_| |_| \__||_| \__,_||_| /_/   \_\|_| |_| \__,_||_| \__, ||___/|_||___/ |_|    \__,_||_|    \__|
    #                                                                                           |___/                                       
    # Differential Analysis
    ###############################################################################################################

    def _compute_data_for_differential_analysis(self):
        """
        Compute the data required for differential analysis
        """
        
        self.ddt = self.difference_distribution_table()
        self._diff_spectrum = set([self.ddt[i][j] for i in range(2**self.input_size()) for j in range(2**self.output_size())]) - {0, 2**self.input_size()}
        self._diff_spectrum = sorted(list(self._diff_spectrum))
        self._len_diff_spectrum = len(self._diff_spectrum)
        self._diff_weights = [abs(float(log(d/(2**self.input_size()), 2))) for d in self._diff_spectrum]        
        if len(self._diff_spectrum) > 0:
            self._max_diff_weights = int(max(self._diff_weights))
        else:
            self._max_diff_weights = 0
        self._data_required_for_differential_analysis = "Data for differential analysis are computed and stored in memory."
    
    def get_differential_spectrum(self):
        """
        Return the differential spectrum
        """
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        return self._diff_spectrum

    def get_star_ddt(self):
        """
        Generate the start DDT (or 0/1 DDT)
        Star DDT is a 2^m*2^n binary array describing the possibility of differential transitions through the S-box
        """
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        star_ddt = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for dx in range(2**self.input_size()):
            for dy in range(2**self.output_size()):
                if self.ddt[dx][dy] != 0:
                    star_ddt[dx][dy] = 1 #reverse ^ 1
                else:
                    star_ddt[dx][dy] = 0 #reverse
        return star_ddt

    def _star_ddt_to_boolean_function(self, reverse=1):
        """
        Convert the star-DDT into a Boolean function
        """
    
        self.star_ddt = self.get_star_ddt()
        boolean_func = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for dy in range(2**self.output_size()):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.star_ddt[dx][dy] ^ reverse                
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

        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()        
        boolean_function = dict()
        complexity = self.input_size() + self.output_size()        
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for dy in range(2**self.output_size()):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.output_size()))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                if self.ddt[dx][dy] in self._diff_spectrum:
                    p = tuple([int(i == self.ddt[dx][dy]) for i in self._diff_spectrum])
                    key = x + y + p
                    boolean_function[key] = 1
                elif self.ddt[dx][dy] == 2**self.input_size():
                    key = x + y + tuple([0]*self._len_diff_spectrum)
                    boolean_function[key] = 1        
        if reverse == 1:            
            complexity = self.input_size() + self.output_size() + self._len_diff_spectrum            
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
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        if any([i != int(i) for i in self._diff_weights]):
            raise ValueError("All transition's weights should be integers")
        boolean_function = dict()
        complexity = self.input_size() + self.output_size()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for dy in range(2**self.output_size()):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.output_size()))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                if self.ddt[dx][dy] != 0:
                    w = int(abs(float(log(self.ddt[dx][dy]/(2**self.input_size()), 2))))
                    p = tuple([0]*(self._max_diff_weights - w) + [1]*w)
                    key = x + y + p
                    boolean_function[key] = 1        
        if reverse == 1:            
            complexity = self.input_size() + self.output_size() + self._len_diff_spectrum
            for dx in range(2**complexity):
                x = tuple(map(int, list(bin(dx)[2:].zfill(complexity))))
                boolean_function[x] = boolean_function.get(x, 0) ^ 1
        return boolean_function

    def _pddt_to_booleanfunction(self, p, reverse=1):
        """
        Convert the p-DDT into a Boolean function
        """
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        boolean_function = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for dy in range(2**self.output_size()):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.output_size()))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                key = x + y
                if self.ddt[dx][dy] == p:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function

    def minimized_diff_constraints(self, mode=6, subtable=None, cryptosmt_compatible=False):
        """
        This method takes a given Boolean function and records its truth table in a file, 
        adhering to the ESPRESSO input format. It then invokes ESPRESSO to obtain a minimized 
        representation of the function. Following that, it interprets ESPRESSO's output and 
        converts the simplified representation into the language recognized by MILP and SAT solvers.
        """
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        valid_values_for_subtable = ["star", None] + list(self._diff_spectrum)
        if subtable not in valid_values_for_subtable:
            raise ValueError("Invalid value for subtable! subtable must be in {0}.".format(list(map(str, valid_values_for_subtable))))
    
        self.cryptosmt_compatible = cryptosmt_compatible
        self.ddt_subtable = subtable

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1
        
        self.diff_objective = ""
        if subtable == "star":
            boolean_function = self._star_ddt_to_boolean_function(reverse=reverse)
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())]
        elif subtable in self._diff_spectrum:
            boolean_function = self._pddt_to_booleanfunction(p=subtable, reverse=reverse)
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())]
        elif cryptosmt_compatible:
            boolean_function = self._ddt_to_cryptosmt_compatible_boolean_function(reverse=reverse)
            if self._max_diff_weights == 0:
                self.diff_objective = "0"
            else:
                self.diff_objective = ["p{:d}".format(i) for i in range(self._max_diff_weights)]
            self.diff_objective = "\nWeight: {}".format(" + ".join(self.diff_objective))
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())] + \
                                     [f"p{i}" for i in range(self._max_diff_weights)]
        else:
            boolean_function = self._ddt_to_boolean_function(reverse=reverse)
            if self._max_diff_weights == 0:
                self.diff_objective = "0"
            else:
                self.diff_objective = ["{:0.04f} p{:d}".format(self._diff_weights[i], i) for i in range(self._len_diff_spectrum)]
            self.diff_objective = "\nWeight: {}".format(" + ".join(self.diff_objective))
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())] + \
                                     [f"p{i}" for i in range(self._len_diff_spectrum)]
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        print("{}".format(self.variables_mapping + self.diff_objective))
        return sat_clauses, milp_constraints

    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    #  _      _                                _                   _              _        ____               _   
    # | |    (_) _ __    ___   __ _  _ __     / \    _ __    __ _ | | _   _  ___ (_) ___  |  _ \  __ _  _ __ | |_ 
    # | |    | || '_ \  / _ \ / _` || '__|   / _ \  | '_ \  / _` || || | | |/ __|| |/ __| | |_) |/ _` || '__|| __|
    # | |___ | || | | ||  __/| (_| || |     / ___ \ | | | || (_| || || |_| |\__ \| |\__ \ |  __/| (_| || |   | |_ 
    # |_____||_||_| |_| \___| \__,_||_|    /_/   \_\|_| |_| \__,_||_| \__, ||___/|_||___/ |_|    \__,_||_|    \__|
    #                                                                 |___/                                       
    # Linear Analysis 
    ###############################################################################################################

    def _compute_data_for_linear_analysis(self):
        """
        Compute the data required for linear analysis
        """

        input_size = self.input_size()
        output_size = self.output_size()
        self.lat = self.linear_approximation_table(scale='correlation')        
        self._squared_lat = [[x**2 for x in y] for y in self.lat]
        self.lat_scaled_by_absolute_correlation = [[2**input_size * lxy for lxy in ly] for ly in self.lat]
        self._squared_correlation_spectrum = sorted(list(set(flatten(self._squared_lat)) - {0, 1}))
        self._linear_weights = [abs(float(log(x, 2))) for x in self._squared_correlation_spectrum]
        self._len_linear_weights = len(self._linear_weights)
        if self._len_linear_weights > 0:            
            self._max_linear_weights = int(max(self._linear_weights))
        else:
            self._max_linear_weights = 0
        self._data_required_for_linear_analysis = "Data for linear analysis are computed and stored in memory."
    
    def get_squared_lat(self):
        """
        Return the squared LAT
        """
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        return self._squared_lat

    def get_squared_correlation_spectrum(self):
        """
        Return the squared correlation spectrum
        """
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        return self._squared_correlation_spectrum

    def get_star_lat(self):
        """
        Generate the start LAT (or 0/1 LAT)
        Star LAT is a 2^m*2^n binary array describing the possibility of linear transitions through the S-box
        """
        
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        star_lat = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for lx in range(2**self.input_size()):
            for ly in range(2**self.output_size()):
                if self.lat[lx][ly] != 0:
                    star_lat[lx][ly] = 1 #reverse ^ 1
                else:
                    star_lat[lx][ly] = 0 #reverse
        return star_lat
    
    def _star_lat_to_boolean_function(self, reverse=1):
        """
        Convert the star-LAT into a Boolean function
        """
 
        self.star_lat = self.get_star_lat()
        boolean_func = dict()
        for lx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(lx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.star_lat[lx][ly] ^ reverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse
        return boolean_func

    def _sqlat_to_boolean_function(self, reverse=1):
        """
        Convert the squared LAT into a Boolean function
        To encode the squared correlations, we define as many new binary variables as the 
        number of different (non-zero and non-one) entries in squared LAT and 
        then encode the whole of LAT as a Boolean function. 
        Let lx, and ly denote the input and output masks, respectively. 
        Assuming that the squared LAT of S-box includes for example three (non-zero and non-one) elements: 
        e0, e1, and e2, we define three new binary variables p0, p1, and p2 and then 
        encode the squared LAT as a Boolean function f such that f(Bin(lx) || Bin(ly) || p0 || p1 || p2) = 1
        if and only if:
               DDT[lx][ly] = e0 and (p0, p1, p2) = (1, 0, 0),
            or DDT[lx][ly] = e1 and (p0, p1, p2) = (0, 1, 0),
            or DDT[lx][ly] = e2 and (p0, p1, p2) = (0, 0, 1),
            or lx = ly = 0      and (p0, p1, p2) = (0, 0, 0),
        otherwise, f(Bin(lx) || Bin(ly) || p0 || p1 || p2) = 0.

        If reverse = True, then we compute the complement of the derived Boolean function.
        """

        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()        
        boolean_function = dict()
        complexity = self.input_size() + self.output_size()        
        for lx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(lx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                if self._squared_lat[lx][ly] in self._squared_correlation_spectrum:
                    p = tuple([int(i == self._squared_lat[lx][ly]) for i in self._squared_correlation_spectrum])
                    key = x + y + p
                    boolean_function[key] = 1
                elif self._squared_lat[lx][ly] == 1:
                    key = x + y + tuple([0]*self._len_linear_weights)
                    boolean_function[key] = 1        
        if reverse == 1:            
            complexity = self.input_size() + self.output_size() + self._len_linear_weights            
            for lx in range(2**complexity):
                x = tuple(map(int, list(bin(lx)[2:].zfill(complexity))))
                boolean_function[x] = boolean_function.get(x, 0) ^ 1
        return boolean_function

    def _psqlat_to_booleanfunction(self, p, reverse=1):
        """
        Convert the p-SquaredLAT into a Boolean function
        """
        
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        boolean_function = dict()
        for lx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(lx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                key = x + y
                if self._squared_lat[lx][ly] == p:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function

    def minimized_linear_constraints(self, mode=6, subtable=None):
        """
        This method takes a given Boolean function and records its truth table in a file, 
        adhering to the ESPRESSO input format. It then invokes ESPRESSO to obtain a minimized 
        representation of the function. Following that, it interprets ESPRESSO's output and 
        converts the simplified representation into the language recognized by MILP and SAT solvers.
        """
        
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        valid_values_for_subtable = ["star", None] + list(self._squared_correlation_spectrum)
        if subtable not in valid_values_for_subtable:
            raise ValueError("Invalid value for subtable! subtable must be in {0}.".format(list(map(str, valid_values_for_subtable))))

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1
        
        self.linear_objective = ""
        if subtable == "star":
            boolean_function = self._star_lat_to_boolean_function(reverse=reverse)
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())]
        elif subtable in self._squared_correlation_spectrum:
            boolean_function = self._psqlat_to_booleanfunction(p=subtable, reverse=reverse)
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())]
        else:
            boolean_function = self._sqlat_to_boolean_function(reverse=reverse)
            if self._len_linear_weights != []:                
                self.linear_objective = ["{:0.04f} p{:d}".format(self._linear_weights[i], i) for i in range(self._len_linear_weights)]
            else:
                self.linear_objective = "0"
            self.linear_objective = "\nWeight: {}".format(" + ".join(self.linear_objective))
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())] + \
                                     [f"p{i}" for i in range(self._len_linear_weights)]
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        print("{}".format(self.variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints
    
    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    #  ___         _                            _      _                   _              _        ____               _   
    # |_ _| _ __  | |_  ___   __ _  _ __  __ _ | |    / \    _ __    __ _ | | _   _  ___ (_) ___  |  _ \  __ _  _ __ | |_ 
    #  | | | '_ \ | __|/ _ \ / _` || '__|/ _` || |   / _ \  | '_ \  / _` || || | | |/ __|| |/ __| | |_) |/ _` || '__|| __|
    #  | | | | | || |_|  __/| (_| || |  | (_| || |  / ___ \ | | | || (_| || || |_| |\__ \| |\__ \ |  __/| (_| || |   | |_ 
    # |___||_| |_| \__|\___| \__, ||_|   \__,_||_| /_/   \_\|_| |_| \__,_||_| \__, ||___/|_||___/ |_|    \__,_||_|    \__|
    #                        |___/                                            |___/                                          
    # Integral Analysis 
    ###############################################################################################################

    def _compute_data_for_integral_analysis(self):
        """
        Compute the data required for integral analysis 
        """

        # compute monomial prediction table
        self.mpt = self.monomial_prediction_table()
    
    def _mpt_to_boolean_function(self, reverse=1):
        """
        Convert the star-LAT into a Boolean function
        """
 
        if self._data_required_for_integral_analysis is None:
            self._compute_data_for_integral_analysis()
        boolean_func = dict()
        for mx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(mx)[2:].zfill(self.input_size()))))
            for my in range(2**self.output_size()):
                y = tuple(map(int, list(bin(my)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.mpt[mx][my] ^ reverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse
        return boolean_func

    def minimized_integral_constraints(self, mode=6):
        """
        This method takes a given Boolean function and records its truth table in a file, 
        adhering to the ESPRESSO input format. It then invokes ESPRESSO to obtain a minimized 
        representation of the function. Following that, it interprets ESPRESSO's output and 
        converts the simplified representation into the language recognized by MILP and SAT solvers.
        """
        
        if self._data_required_for_integral_analysis is None:
            self._compute_data_for_integral_analysis()        

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1
        
        self.integral_objective = ""
        input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                 [f"b{i}" for i in range(self.output_size())]
        boolean_function = self._mpt_to_boolean_function(reverse=reverse)
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        print("{}".format(self.variables_mapping + self.integral_objective))
        return sat_clauses, milp_constraints
    
    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    #  ____         _                          _         _       _    _         ____         _                    _              
    # |  _ \   ___ | |_  ___  _ __  _ __ ___  (_) _ __  (_) ___ | |_ (_)  ___  | __ )   ___ | |__    __ _ __   __(_)  ___   _ __ 
    # | | | | / _ \| __|/ _ \| '__|| '_ ` _ \ | || '_ \ | |/ __|| __|| | / __| |  _ \  / _ \| '_ \  / _` |\ \ / /| | / _ \ | '__|
    # | |_| ||  __/| |_|  __/| |   | | | | | || || | | || |\__ \| |_ | || (__  | |_) ||  __/| | | || (_| | \ V / | || (_) || |   
    # |____/  \___| \__|\___||_|   |_| |_| |_||_||_| |_||_||___/ \__||_| \___| |____/  \___||_| |_| \__,_|  \_/  |_| \___/ |_|   
    # 
    # Deteministic Behavior
    ###############################################################################################################
    
    def truncated_to_binvectors(self, input_vector):
        """
        Converts a truncated vector to a list of binary vectors
        """
        
        output = list(map(list, (itertools.product(*[self.deterministic_mask[i] for i in input_vector]))))
        return output

    def binvectors_to_truncated(self, input_list):
        """
        Converts a list of binary vectors to a truncated vector
        """
        
        vector_size = len(input_list[0])
        # Transpose the list so we can check each bit position across all numbers
        transposed_list = list(zip(*input_list))
        output = [self.zero]*vector_size
        for i in range(vector_size):
            if (-1) in list(map(int, transposed_list[i])):
                # If -1 is present, set the corresponding output coordinate to unknown
                output[i] = self.unknown
            elif 0 in transposed_list[i] and 1 in transposed_list[i]:
                # If both 0 and 1 are present, set the corresponding output coordinate to unknown
                output[i] = self.unknown
            elif 1 in transposed_list[i]:
                # If only 1 is present, set the corresponding output coordinate to 1
                output[i] = self.one
            # If only 0 is present, the output coordinate remains 0
        return output

    def encode_deterministic_differential_behavior(self):
        """
        Encodes the deterministic differential behavior of the S-box
        """

        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()                
        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()
        propagation_dictionary = {}
        for truncated_input_diff in itertools.product(list(self.deterministic_mask.keys()), repeat=sbox_input_size):        
            binary_input_diffs = self.truncated_to_binvectors(truncated_input_diff)
            possible_output_diffs = []
            for input_diff in binary_input_diffs:
                input_diff = int("".join(str(x) for x in input_diff), base=2)
                possible_output_diffs += [self.to_bits(x=dy, n=sbox_output_size) for dy in range(2**sbox_output_size) if self.ddt[input_diff][dy] != 0]
            truncated_output_diff = self.binvectors_to_truncated(possible_output_diffs)
            if any(x != self.unknown for x in truncated_output_diff):
                propagation_dictionary[truncated_input_diff] = truncated_output_diff
        return propagation_dictionary

    def encode_deterministic_linear_behavior(self):
        """
        Encodes the deterministic linear behavior of the S-box        
        """

        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()                
        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()
        propagation_dictionary = {}
        for truncated_input_mask in itertools.product(list(self.deterministic_mask.keys()), repeat=sbox_input_size):        
            binary_input_masks = self.truncated_to_binvectors(truncated_input_mask)
            possible_output_masks = []
            for input_mask in binary_input_masks:
                input_mask = int("".join(str(x) for x in input_mask), base=2)
                output_mask = [self.to_bits(x=ly, n=sbox_output_size) for ly in range(2**sbox_output_size) if self.lat[input_mask][ly] != 0]
                if output_mask == []:
                    output_mask = [[-1]*sbox_output_size]
                possible_output_masks += output_mask            
            truncated_output_mask = self.binvectors_to_truncated(possible_output_masks)
            if any(x != self.unknown for x in truncated_output_mask):
                propagation_dictionary[truncated_input_mask] = truncated_output_mask            
        return propagation_dictionary

    def encode_deterministic_integral_forward(self):
        """
        Encodes the deterministic integral behavior of the S-box        
        """

        if self._data_required_for_integral_analysis is None:
            self._compute_data_for_integral_analysis()
        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()
        propagation_dictionary = {}
        for truncated_input_mask in itertools.product(list(self.deterministic_mask.keys()), repeat=sbox_input_size):
            binary_input_masks = self.truncated_to_binvectors(truncated_input_mask)
            possible_output_masks = []
            for input_mask in binary_input_masks:
                input_mask = int("".join(str(x) for x in input_mask), base=2)
                output_mask = [self.to_bits(x=my, n=self.output_size()) for my in range(2**sbox_output_size) if self.mpt[input_mask][my] != 0]
                if output_mask == []:
                    possible_output_masks += [[self.unknown]*sbox_output_size]
                else:
                    possible_output_masks += output_mask
            truncated_output_mask = self.binvectors_to_truncated(possible_output_masks)
            if any(x != self.unknown for x in truncated_output_mask):
                propagation_dictionary[truncated_input_mask] = truncated_output_mask
        return propagation_dictionary

    def encode_deterministic_integral_backward(self):
        """
        Encodes the deterministic integral behavior of the S-box        
        """

        if self._data_required_for_integral_analysis is None:
            self._compute_data_for_integral_analysis()
        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()
        propagation_dictionary = {}
        for truncated_output_mask in itertools.product(list(self.deterministic_mask.keys()), repeat=sbox_output_size):
            binary_output_masks = self.truncated_to_binvectors(truncated_output_mask)
            possible_input_masks = []
            for output_mask in binary_output_masks:
                output_mask = int("".join(str(x) for x in output_mask), base=2)
                input_mask = [self.to_bits(x=mx, n=self.input_size()) for mx in range(2**sbox_input_size) if self.mpt[mx][output_mask] != 0]
                if input_mask == []:
                    possible_input_masks += [[self.unknown]*sbox_input_size]
                else:
                    possible_input_masks += input_mask
            truncated_input_mask = self.binvectors_to_truncated(possible_input_masks)
            if any(x != self.unknown for x in truncated_input_mask):
                propagation_dictionary[truncated_output_mask] = truncated_input_mask
        return propagation_dictionary
    
    def generate_cp_constraints(self, propagation_dictionary):
        """
        Generates the constraints for the CP model
        """

        inputs = list(propagation_dictionary.keys())
        outputs = list(propagation_dictionary.values())
        m = len(inputs[0])
        n = len(outputs[0])
        # a = [f"a{m - i - 1}" for i in range(m)]
        # b = [f"b{n - i - 1}" for i in range(n)]        
        a = [f"a{i}" for i in range(m)]
        b = [f"b{i}" for i in range(n)]
        last_condition = " /\\ ".join([f"b{i} = -1" for i in range(n)])
        constraints = ""
        for i, (input, output) in enumerate(propagation_dictionary.items()):
            input_str = [f"{a[i]} == {input[i]}" for i in range(m)]
            input_str = " /\\ ".join(input_str)
            output_str = [f"{b[i]} = {output[i]}" for i in range(n)]
            output_str = " /\\ ".join(output_str)
            if i == 0:
                constraints += f"if ({input_str}) then ({output_str})\n"
            else:
                constraints += f"elseif ({input_str}) then ({output_str})\n"
        constraints += f"else ({last_condition})\nendif"
        print(self.variables_mapping)
        return constraints

    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    #  ____   _   __   __                          _    _         _         _      _                                _                   _              _      
    # |  _ \ (_) / _| / _|  ___  _ __  ___  _ __  | |_ (_)  __ _ | |       | |    (_) _ __    ___   __ _  _ __     / \    _ __    __ _ | | _   _  ___ (_) ___ 
    # | | | || || |_ | |_  / _ \| '__|/ _ \| '_ \ | __|| | / _` || | _____ | |    | || '_ \  / _ \ / _` || '__|   / _ \  | '_ \  / _` || || | | |/ __|| |/ __|
    # | |_| || ||  _||  _||  __/| |  |  __/| | | || |_ | || (_| || ||_____|| |___ | || | | ||  __/| (_| || |     / ___ \ | | | || (_| || || |_| |\__ \| |\__ \
    # |____/ |_||_|  |_|   \___||_|   \___||_| |_| \__||_| \__,_||_|       |_____||_||_| |_| \___| \__,_||_|    /_/   \_\|_| |_| \__,_||_| \__, ||___/|_||___/
    #                                                                                                                                      |___/                  
    # Differential-Linear (DL) Analysis
    ###############################################################################################################    

    @staticmethod
    def dot_product(x, y):
        # Compute bitwise AND of x and y
        bitwise_and = x & y

        # Count the number of set bits in bitwise_and
        count = bin(bitwise_and).count('1')

        # Return the result modulo 2
        return count % 2

    def _compute_data_for_difflin_analysis(self):
        """
        Compute the data required for differential-linear analysis
        """        

        input_size = self.input_size()
        output_size = self.output_size()

        dlct = [[0 for y in range(2**output_size)] for x in range(2**input_size)]
        for diff in range(2**input_size):
            for mask in range(2**output_size):
                counter_0 = 0
                counter_1 = 0
                for data in range(2**input_size):
                    actual_output_diff = self[data] ^ self[data ^ diff]
                    if self.dot_product(mask, actual_output_diff) == 0:
                        counter_0 += 1 
                    else:
                        counter_1 += 1
                dlct[diff][mask] = counter_0 - counter_1
        self._dlct = dlct
        self._dlct_spectrum = sorted(list(set(flatten(self._dlct))))
        self._dlct_weights = [abs(float(log(abs(x), 2))) for x in self._dlct_spectrum if x != 0]
        self._len_dlct_weights = len(self._dlct_weights)
        self._data_required_for_difflin_analysis = "Data for linear analysis are computed and stored in memory."

    def differential_linear_connectivity_table(self):
        """
        Compute the differential-linear connectivity table
        """
        
        if self._data_required_for_difflin_analysis == None:
            self._compute_data_for_difflin_analysis()
        return self._dlct

            
    def upper_differential_linear_connectivity_table(self):
        """
        Compute the upper differential-linear connectivity table        
        """

        input_size = self.input_size()
        output_size = self.output_size()
        self._upper_differential_linear_connectivity_table = [[[0 for _ in range(2**output_size)] for _ in range(2**output_size)] for _ in range(2**input_size)]        
        for input_diff in range(2**input_size):
            for output_diff in range(2**output_size):
                for output_mask in range(2**output_size):
                    counter_0 = 0
                    counter_1 = 0
                    dotproduct = self.dot_product(output_mask, output_diff)
                    for data in range(2**input_size):                        
                        actual_output_diff = self[data] ^ self[data ^ input_diff]
                        if dotproduct == 0 and actual_output_diff == output_diff:
                            counter_0 += 1
                        elif dotproduct == 1 and actual_output_diff == output_diff:
                            counter_1 += 1
                    self._upper_differential_linear_connectivity_table[input_diff][output_diff][output_mask] = counter_0 - counter_1
        return self._upper_differential_linear_connectivity_table

    def lower_differential_linear_connectivity_table(self):
        """
        Compute the lower differential-linear connectivity table        
        """

        input_size = self.input_size()
        output_size = self.output_size()
        self._lower_differential_linear_connectivity_table = [[[0 for _ in range(2**output_size)] for _ in range(2**input_size)] for _ in range(2**input_size)]
        for input_diff in range(2**input_size):
            for input_mask in range(2**input_size):
                dotproduct_input = self.dot_product(input_mask, input_diff)
                for output_mask in range(2**output_size):
                    counter_0 = 0
                    counter_1 = 0                    
                    for data in range(2**input_size):
                        actual_output_diff = self[data] ^ self[data ^ input_diff]
                        dotproduct_output = self.dot_product(output_mask, actual_output_diff)
                        if dotproduct_input ^ dotproduct_output == 0:
                            counter_0 += 1
                        elif dotproduct_input ^ dotproduct_output == 1:
                            counter_1 += 1                        
                    self._lower_differential_linear_connectivity_table[input_diff][input_mask][output_mask] = counter_0 - counter_1
        return self._lower_differential_linear_connectivity_table

    def double_differential_linear_connectivity_table(self):
        """
        Compute the double differential-linear connectivity table
        """
        
        input_size = self.input_size()
        output_size = self.output_size()
        if self._data_required_for_difflin_analysis == None:
            self._compute_data_for_difflin_analysis()
        if self._data_required_for_differential_analysis == None:
            self._compute_data_for_differential_analysis()
        if self._data_required_for_linear_analysis == None:
            self._compute_data_for_linear_analysis()
        dlct = self.differential_linear_connectivity_table()
        ddt = self.ddt
        self._ddlct = [[0 for y in range(2**output_size)] for x in range(2**input_size)]
        for diff in range(2**input_size):
            for mask in range(2**output_size):
                for diff_middle in range(2**input_size):
                    self._ddlct[diff][mask] += ddt[diff][diff_middle] * dlct[diff_middle][mask]        
        return self._ddlct
    
    def check_hadipour_theorem(self):
        """
        Check the Hadipour et al. theorem
        """

        input_size = self.input_size()
        output_size = self.output_size()
        if self._data_required_for_difflin_analysis == None:
            self._compute_data_for_difflin_analysis()
        if self._data_required_for_differential_analysis == None:
            self._compute_data_for_differential_analysis()
        if self._data_required_for_linear_analysis == None:
            self._compute_data_for_linear_analysis()
        correctness_right = True
        correctness_left = True
        ddt = self.ddt
        dlct = self.differential_linear_connectivity_table()
        udlct = self.upper_differential_linear_connectivity_table()
        ldlct = self.lower_differential_linear_connectivity_table()
        squared_lat = [[a**2 for a in column] for column in self.lat_scaled_by_absolute_correlation]
        for input_diff in range(2**input_size):
            for output_mask in range(2**output_size):
                sum_0 = 0
                for diff_middle in range(2**input_size):                    
                    for mask_middle in range(2**output_size):
                        sum_0 += udlct[input_diff][diff_middle][mask_middle] * ldlct[diff_middle][mask_middle][output_mask]
                sum_0 = sum_0//(2**input_size)
                sum_1 = 0
                for diff_middle in range(2**input_size):
                    sum_1 += ddt[input_diff][diff_middle] * dlct[diff_middle][output_mask]
                sum_2 = 0
                for mask_middle in range(2**output_size):
                    sum_2 += dlct[input_diff][mask_middle] * squared_lat[mask_middle][output_mask]
                sum_2 = sum_2//(2**input_size)
                if sum_1 != sum_2:
                    correctness_right = False
                if sum_0 != sum_1 and sum_1 == sum_2:
                    correctness_left = False
        if correctness_left and correctness_right:
            print("The Hadipour et al.'s theorem is satisfied.")
            return True
        else:
            print("The Hadipour et al.'s theorem is not satisfied.")
            return False
    
    def get_dlct_spectrum(self):
        """
        Compute the set of different entries in DLCT
        """

        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        return self._dlct_spectrum

    def compute_star_dlct(self, reverse=1):
        """
        Generate the start DLCT (or 0/1 DLCT)
        Star DLCT is a 2^m*2^n binary array describing the possibility of differential-linear transitions through the S-box
        """
        
        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        self.star_dlct = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for dx in range(2**self.input_size()):
            for ly in range(2**self.output_size()):
                if self._dlct[dx][ly] != 0:
                    self.star_dlct[dx][ly] = reverse ^ 1
                else:
                    self.star_dlct[dx][ly] = reverse
        return self.star_dlct

    def _dlct_to_booleanfunction(self, corr, reverse=1):
        """
        Convert a subtable of DLCT into a Boolean function
        """
        
        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        boolean_function = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                key = x + y
                if self._dlct[dx][ly] == corr:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function

    def _star_dlct_to_boolean_function(self, reverse=1):
        """
        Convert the star-DLCT into a Boolean function
        """
 
        self.compute_star_dlct(reverse=reverse)
        boolean_func = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.star_dlct[dx][ly] ^ reverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse
        return boolean_func

    def _star_dlct_to_characteristic_boolean_function(self, reverse=1, inverse=0):
        """
        Convert the star-inverse-DLCT into a Boolean function
        """
 
        self.compute_star_dlct(reverse=reverse)
        boolean_func = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                for s in range(2):
                    key = x + y + (s, )
                    if self.star_dlct[dx][ly] == s:
                        boolean_func[key] = 1 ^ reverse ^ inverse
                    else:
                        boolean_func[key] = 0 ^ reverse ^ inverse
        return boolean_func
    
    def minimized_differential_linear_constraints(self, mode=6, subtable=None):
        """
        This method takes a given Boolean function and records its truth table in a file, 
        adhering to the ESPRESSO input format. It then invokes ESPRESSO to obtain a minimized 
        representation of the function. Following that, it interprets ESPRESSO's output and 
        converts the simplified representation into the language recognized by MILP and SAT solvers.
        """
        
        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        valid_values_for_subtable = ["star", "star_inverse", None] + list(self._dlct_spectrum)
        if subtable not in valid_values_for_subtable:
            raise ValueError("Invalid value for subtable! subtable must be in {0}.".format(list(map(str, valid_values_for_subtable))))

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1
        
        self.linear_objective = ""
        if subtable == "star_inverse":
            boolean_function = self._star_dlct_to_characteristic_boolean_function(reverse=reverse, inverse=1)
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())] + ["s"]

        elif subtable == "star":
            boolean_function = self._star_dlct_to_boolean_function(reverse=reverse)
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())]
        elif subtable in self._dlct_spectrum:
            boolean_function = self._dlct_to_booleanfunction(corr=subtable, reverse=reverse)
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())]
        else:
            boolean_function = self._dlct_to_booleanfunction(reverse=reverse)
            if self._len_dlct_weights != []:                
                self.linear_objective = ["{:0.04f} p{:d}".format(self._dlct_weights[i], i) for i in range(self._dlct_weights)]
            else:
                self.linear_objective = "0"
            self.linear_objective = "\nWeight: {}".format(" + ".join(self.linear_objective))
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())] + \
                                     [f"p{i}" for i in range(self._len_linear_weights)]
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        print("{}".format(self.variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints

###############################################################################################################
###############################################################################################################
###############################################################################################################
#  _____           _     _____                     _    _                       _  _  _          
# |_   _|___  ___ | |_  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __    __ _ | |(_)| |_  _   _ 
#   | | / _ \/ __|| __| | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \  / _` || || || __|| | | |
#   | ||  __/\__ \| |_  |  _| | |_| || | | || (__ | |_ | || (_) || | | || (_| || || || |_ | |_| |
#   |_| \___||___/ \__| |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_| \__,_||_||_| \__| \__, |
#                                                                                          |___/ 
############################################################################################################### 

if __name__ == "__main__":
    line_seperator = '\n' + '#' * 46 + '\n'    
    # example for differential analysis
    from sage.crypto.sboxes import SBox
    from sage.crypto.sboxes import CLEFIA_S0 as sb
    # sb = SBox([sb(i ^ 0x8) ^ sb(i) for i in range(16)])    
    # sb = SBox([0, 1, 1, 0])
    # Orthros's Sbox    
    # sb = SBox([1, 0, 2, 4, 3, 8, 6, 0xd, 9, 0xa, 0xb, 0xe, 0xf, 0xc, 7, 5]) 
    # Spongent's Sbox
    # sb = SBox([0xe, 0xd, 0xb, 0, 2, 1, 4, 0xf, 7, 0xa, 8, 5, 9, 0xc, 3, 6]) 
    # SPEEDY's S-box
    # sb = SBox([0x08, 0x00, 0x09, 0x03, 0x38, 0x10, 0x29, 0x13, 0x0C, 0x0D, 0x04, 0x07, 0x30, 0x01, 0x20, 0x23, 0x1A, 0x12, 0x18, 0x32, 0x3E, 0x16, 0x2C, 0x36, 0x1C, 0x1D, 0x14, 0x37, 0x34, 0x05, 0x24, 0x27, 0x02, 0x06, 0x0B, 0x0F, 0x33, 0x17, 0x21, 0x15, 0x0A, 0x1B, 0x0E, 0x1F, 0x31, 0x11, 0x25, 0x35, 0x22, 0x26, 0x2A, 0x2E, 0x3A, 0x1E, 0x28, 0x3C, 0x2B, 0x3B, 0x2F, 0x3F, 0x39, 0x19, 0x2D, 0x3D])    
    # sb = sb.inverse()
    sa = SboxAnalyzer(sb)
    print(sa)
    # diff_spec = sa.get_differential_spectrum()
    # print(f"differential spectrum: {diff_spec}")
    # cnf, milp = sa.minimized_diff_constraints()
    # print("\nencode DDT")
    # pretty_print(milp)

    # star_ddt  = sa.get_star_ddt()
    # print("\nStar DDT:")
    # pretty_print(star_ddt)
    # cnf, milp = sa.minimized_diff_constraints(subtable='star')
    # print("\nencode *-DDT")
    # pretty_print(milp)

    # with open("milp.lp", "w") as fileobj:
    #     fileobj.write("\nsubject to\n")
    #     fileobj.write("\n".join(milp))
    #     fileobj.write("\nend")
    # with open("cnf.txt", "w") as fileobj:
    #     fileobj.write(cnf)    
    
    # print(line_seperator)

    # # example for linear analysis
    # squared_lat = sa.get_squared_lat()
    # print("\nSquared LAT:")
    # pretty_print(squared_lat)
    # cnf, milp = sa.minimized_linear_constraints()
    # print("\nencode Squared-LAT")
    # pretty_print(milp)
    # star_lat = sa.get_star_lat()
    # print("\nStar LAT:")
    # cnf, milp = sa.minimized_linear_constraints(subtable="star", mode=6)
    # print("\nencode LAT")
    # pretty_print(milp)
    lin_spec = sa.get_squared_correlation_spectrum()
    stn = 5
    print(lin_spec)
    print("\nencode {0}-LAT".format(lin_spec[stn]))
    print(float(log(lin_spec[stn], 2)))
    cnf, milp = sa.minimized_linear_constraints(subtable=lin_spec[stn], mode=7)

    with open("milp.lp", "w") as fileobj:
        # fileobj.write("\nsubject to\n")
        fileobj.write("\n".join(milp))
        # fileobj.write("\nend")
    # with open("cnf.txt", "w") as fileobj:
    #     fileobj.write(cnf) 
    # print(line_seperator)

    # # compute monomial prediction table
    # print("\ncompute monomial prediction table")
    # start_time = time.time()
    # mpt = sa.monomial_prediction_table()
    # elapsed_time = time.time() - start_time
    # print("\nTime to compute monomial prediction table: {0:0.2f} seconds".format(elapsed_time))
    # pretty_print(mpt)
    # # example for integral analysis
    # print("\nencode MPT")
    # start_time = time.time()
    # cnf, milp = sa.minimized_integral_constraints(mode=2)
    # pretty_print(milp)
    # elapsed_time = time.time() - start_time
    # print("\nTime to encode MPT: {0:0.2f} seconds".format(elapsed_time))

    # with open("milp.lp", "w") as fileobj:
    #     fileobj.write("\nsubject to\n")
    #     fileobj.write("\n".join(milp))
    #     fileobj.write("\nend")
    # with open("cnf.txt", "w") as fileobj:
    #     fileobj.write(cnf) 

    # example for deterministic differential analysis
    # print("\nEncoding deterministic differential behavior ...\n")
    # start_time = time.time()
    # diff_propagation = sa.encode_deterministic_differential_behavior()
    # elapsed_time = time.time() - start_time
    # print("\nTime to generate differential behavior: {0:0.2f} seconds".format(elapsed_time))
    # print("\nDifferential behavior:")
    # pretty_print(diff_propagation)
    # cp_constraints = sa.generate_cp_constraints(diff_propagation)
    # print("\nencode deterministic differential behavior")
    # print(cp_constraints)
    # print(line_seperator)

    # # example for deterministic linear analysis
    # print("\nEncoding deterministic linear behavior ...\n")
    # start_time = time.time()
    # lin_propagation = sa.encode_deterministic_linear_behavior()
    # elapsed_time = time.time() - start_time
    # print("\nLinear behavior:\n")
    # pretty_print(lin_propagation)
    # print("\nTime to generate linear behavior: {0:0.2f} seconds".format(elapsed_time))
    # cp_constraints = sa.generate_cp_constraints(lin_propagation)
    # print("\nencode deterministic linear behavior")
    # print(cp_constraints)
    # print(line_seperator)

    # # example for deterministic integral analysis
    # print("\nEncoding deterministic integral behavior ...\n")
    # start_time = time.time()
    # int_propagation_forward = sa.encode_deterministic_integral_forward()
    # int_propagation_backward = sa.encode_deterministic_integral_backward()
    # elapsed_time = time.time() - start_time
    # print("\n Deterministic Integral Propagation (Forward):\n")
    # pretty_print(int_propagation_forward)
    # print("\n Deterministic Integral Propagation (Backward):\n")
    # pretty_print(int_propagation_backward)
    # print("\nTime to generate integral behavior: {0:0.2f} seconds".format(elapsed_time))

    # # example for differential-linear analysis
    # print(line_seperator)
    # print("\nDLCT")
    # start_time = time.time()
    # dlct = sa.differential_linear_connectivity_table()
    # sa.print_table(dlct)
    # elapsed_time = time.time() - start_time
    # print("\nTime to compute differential-linear connectivity table: {0:0.2f} seconds".format(elapsed_time))
    # test = sa.get_dlct_spectrum()
    # print("\nDLCT spectrum: {0}".format(test))
    # cnf, milp = sa.minimized_differential_linear_constraints(subtable='star_inverse', mode=6)
    # with open("milp.lp", "w") as fileobj:
    #     fileobj.write("\nsubject to\n")
    #     fileobj.write("\n".join(milp))
    #     fileobj.write("\nend")
    # pretty_print(milp)
    
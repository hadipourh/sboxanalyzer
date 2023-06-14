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
    
    def compute_monomial_prediction_table(self):
        """
        Compute the Monomial Prediction Table (MPT) based on [HE22]:
        https://tosc.iacr.org/index.php/ToSC/article/view/9715
        """

        raise NotImplementedError("This function has not been implemented yet.")

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
        self._max_diff_weights = int(max(self._diff_weights))        
        self._data_required_for_differential_analysis = "Data for differential analysis are computed and stored in memory."
    
    def get_differential_spectrum(self):
        """
        Return the differential spectrum
        """
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        return self._diff_spectrum

    def compute_star_ddt(self, reverse=1):
        """
        Generate the start DDT (or 0/1 DDT)
        Star DDT is a 2^m*2^n binary array describing the possibility of differential transitions through the S-box
        """
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        self.star_ddt = [[0 for i in range(2**self.input_size())] for j in range(2**self.output_size())]
        for dx in range(2**self.input_size()):
            for dy in range(2**self.output_size()):
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
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for dy in range(2**self.output_size()):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.input_size()))))
                key = x + y
                boolean_func[key] = self.star_ddt[dx][dy]
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = reverse ^ 1
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

        :param mode list: a set of flags specifying the configuration of ESPRESSO program
        :param booleanfunction dict: a Python dictionary representing the truth table of the given Boolean function
        :rtype: string
        :return: the minimized MILP/SAT representation of the given Boolean function
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
            self.diff_objective = ["p{:d}".format(i) for i in range(self._max_diff_weights)]
            self.diff_objective = "\nWeight: {}".format(" + ".join(self.diff_objective))
            input_output_variables = [f"a{i}" for i in range(self.input_size())] + \
                                     [f"b{i}" for i in range(self.output_size())] + \
                                     [f"p{i}" for i in range(self._max_diff_weights)]
        else:
            boolean_function = self._ddt_to_boolean_function(reverse=reverse)
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
        self.lat = self.linear_approximation_table(scale='correlation')
        self._squared_lat = [[x**2 for x in y] for y in self.lat]
        self._squared_correlation_spectrum = sorted(list(set(flatten(self._squared_lat)) - {0, 1}))
        self._linear_weights = [abs(float(log(x, 2))) for x in self._squared_correlation_spectrum]
        self._len_linear_weights = len(self._linear_weights)         
        self._max_linear_weights = int(max(self._linear_weights))
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

    def compute_star_lat(self, reverse=1):
        """
        Generate the start LAT (or 0/1 LAT)
        Star LAT is a 2^m*2^n binary array describing the possibility of linear transitions through the S-box
        """
        
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        self.star_lat = [[0 for i in range(2**self.input_size())] for j in range(2**self.output_size())]
        for lx in range(2**self.input_size()):
            for ly in range(2**self.output_size()):
                if self.lat[lx][ly] != 0:
                    self.star_lat[lx][ly] = reverse ^ 1
                else:
                    self.star_lat[lx][ly] = reverse
        return self.star_lat
    
    def _star_lat_to_boolean_function(self, reverse=1):
        """
        Convert the star-LAT into a Boolean function
        """
 
        self.compute_star_lat(reverse=reverse)
        boolean_func = dict()
        for lx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(lx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.input_size()))))
                key = x + y
                boolean_func[key] = self.star_lat[lx][ly]
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = reverse ^ 1
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

        :param mode list: a set of flags specifying the configuration of ESPRESSO program
        :param booleanfunction dict: a Python dictionary representing the truth table of the given Boolean function
        :rtype: string
        :return: the minimized MILP/SAT representation of the given Boolean function
        """
        
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        valid_values_for_subtable = ["star", None] + list(self._squared_correlation_spectrum)
        if subtable not in valid_values_for_subtable:
            raise ValueError("Invalid value for subtable! subtable must be in {0}.".format(list(map(str, valid_values_for_subtable))))

        self.lat_subtable = subtable

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
            self.linear_objective = ["{:0.04f} p{:d}".format(self._linear_weights[i], i) for i in range(self._len_linear_weights)]
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
        self.mpt = self.compute_monomial_prediction_table()
        self._data_required_for_integral_analysis = "Data for integral analysis are computed and stored in memory."
        raise NotImplementedError("This function has not been implemented yet.")    

if __name__ == "__main__":
    # example for differential analysis
    from sage.crypto.sboxes import Ascon as sb
    sa = SboxAnalyzer(sb)
    diff_spec = sa.get_differential_spectrum()
    print(f"differential spectrum: {diff_spec}")
    cnf, milp = sa.minimized_diff_constraints()
    print(f"\nencode {diff_spec[0]}-DDT:")
    cnf2, milp2 = sa.minimized_diff_constraints(subtable=diff_spec[0], mode=5)
    
    # example for linear analysis
    print("\nencode *-LAT")
    cnf, milp = sa.minimized_linear_constraints(subtable="star", mode=5)
    print("\nencode LAT")
    cnf, milp = sa.minimized_linear_constraints()
    lin_spec = sa.get_squared_correlation_spectrum()
    print("\nencode {0}-LAT".format(lin_spec[0]))    
    cnf, milp = sa.minimized_linear_constraints(subtable=lin_spec[0], mode=6)


    with open("milp.lp", "w") as fileobj:
        fileobj.write("\nsubject to\n")
        fileobj.write("\n".join(milp))
        fileobj.write("\nend")
    with open("cnf.txt", "w") as fileobj:
        fileobj.write(cnf)    

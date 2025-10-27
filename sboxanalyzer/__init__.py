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
    sage: cnf, milp, cp = sa.minimized_diff_constraints()

    Simplifying the MILP/SAT constraints ...
    Time used to simplify the constraints: 0.00 seconds
    Number of constraints: 17
    Input: a0||a1||a2; msb: a0
    Output: b0||b1||b2; msb: b0
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
#*****************************************************************************
# Copyright (c) 2022 Hosein Hadipour <hsn.hadipour@gmail.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# 
# In case you use this tool please include the above copyright informations (name, contact, license)
#*****************************************************************************

__version__ = "1.0.2"

import subprocess
import os
import time
from sage.all import *
from sage.crypto.sboxes import SBox
import itertools
from .espresso_loader import get_espresso_path_cached

# Get espresso binary path (uses pre-built binary or builds from source)
ESPRESO_BIN_PATH = get_espresso_path_cached()

class SboxAnalyzer(SBox):
    r"""
    This module encodes the DDT, LAT and MPT [HE22]_ of a given S-box with MILP/SAT constraints and 
    then simplifies the extracted constraints using logic minimization tools

    EXAMPLES:
    
    We consider the S-box of block cipher PRINTcipher [KLPR2010]_::

        sage: from sboxanalyzer import *                                                
        sage: from sage.crypto.sboxes import PRINTcipher as sb                          
        sage: sa = SboxAnalyzer(sb)                                                     
        sage: cnf, milp, cp = sa.minimized_diff_constraints()
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
        # A flag to check if the data required for differential analysis are present in memory
        self._data_required_for_differential_analysis = None
        # A flag to check if the data required for linear analysis are present in memory    
        self._data_required_for_linear_analysis = None
        # A flag to check if the data required for integral analysis are presenet in memory
        self._data_required_for_integral_analysis = None
        # A flag to check if the data required for differential-linear analysis are presenet in memory
        self._data_required_for_difflin_analysis = None
        # A flag to check if the data required for boomerang analysis are presenet in memory
        self._data_required_for_boomerang_analysis = None

        # define a dictionalry to encode deterministic behavior
        self.unknown, self.zero, self.one = -1, 0, 1
        self.zero_binary = [0, 0]
        self.one_binary = [0, 1]
        self.unknown_binary = [1, 0]
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

    @staticmethod
    def _write_truth_table(filename, boolean_function, input_output_variables=None):
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
    
    @staticmethod
    def simplify_by_espresso(input_file, output_file, mode):
        """
        Simplify the CNF formula using the ESPRESSO
        """
        valid_values_for_mode = list(range(8))
        if mode not in valid_values_for_mode:
            raise ValueError("Invalid value for mode! mode must be in [0, 1, 2, 3, 4, 5, 6, 7].") 

                                                                        # If reverse = 1 choose ON-SET and if reverse = 0 choose OFF-SET
        espresso_options =  [[],                                        # 0 ON-SET :   Derived constraints exclude point p such that f(p) = 1
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
        with open(output_file, 'w') as fileobj:
            subprocess.call([ESPRESO_BIN_PATH, *espresso_options[mode], input_file], stdout=fileobj)
        
    @staticmethod
    def _parse_the_output_of_espresso(filename, alphabet):
        """
        Parse the output of ESPRESSO

        :param filename str: an string specifying the filename
        :param alphabet list: a list of strings specifying the variables
        :return: a tuple of strings including the SAT clauses, MILP constraints, and CP constraints
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
        cp_constraints =  "\n/\\\n".join(milp_constraints)
        return (sat_clauses, milp_constraints, cp_constraints)

    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    #  ____                 _                              _____                     _    _                                        _   ____         _             __   ____   _                             __     __           _                    
    # | __ )   ___    ___  | |  ___   __ _  _ __    __ _  |  ___|_   _  _ __    ___ | |_ (_)  ___   _ __   ___    __ _  _ __    __| | / ___|   ___ | |_    ___   / _| | __ ) (_) _ __    __ _  _ __  _   _  \ \   / /___   ___ | |_  ___   _ __  ___ 
    # |  _ \  / _ \  / _ \ | | / _ \ / _` || '_ \  / _` | | |_  | | | || '_ \  / __|| __|| | / _ \ | '_ \ / __|  / _` || '_ \  / _` | \___ \  / _ \| __|  / _ \ | |_  |  _ \ | || '_ \  / _` || '__|| | | |  \ \ / // _ \ / __|| __|/ _ \ | '__|/ __|
    # | |_) || (_) || (_) || ||  __/| (_| || | | || (_| | |  _| | |_| || | | || (__ | |_ | || (_) || | | |\__ \ | (_| || | | || (_| |  ___) ||  __/| |_  | (_) ||  _| | |_) || || | | || (_| || |   | |_| |   \ V /|  __/| (__ | |_| (_) || |   \__ \
    # |____/  \___/  \___/ |_| \___| \__,_||_| |_| \__, | |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_||___/  \__,_||_| |_| \__,_| |____/  \___| \__|  \___/ |_|   |____/ |_||_| |_| \__,_||_|    \__, |    \_/  \___| \___| \__|\___/ |_|   |___/
    #                                              |___/                                                                                                                                             |___/                                           
    # Encoding Booleann functions and set of binary vectors into MILP/SAT constraints

    @staticmethod
    def encode_boolean_function(truth_table, input_variables=None, mode=6):
        """
        Encode the Boolean function into MILP/SAT constraints
        """

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1

        if any([i not in {0, 1} for i in truth_table]):
            raise ValueError("Each element of the truth table must be either 0 or 1.")
        
        log2_length = log(len(truth_table), 2)
        if log2_length != int(log2_length):
            raise ValueError("The length of the truth table must be a power of 2.") 
              
        boolean_function = dict()
        for i in range(len(truth_table)):
            key = tuple(map(int, list(bin(i)[2:].zfill(log2_length))))
            boolean_function[key] = truth_table[i] ^ reverse
        input_file = os.path.join(os.getcwd(), 'tt_' + hex(randint(0, 65536))[2:] + '.txt')
        output_file = os.path.join(os.getcwd(), 'stt_' + hex(randint(0, 65536))[2:] + '.txt')
        if input_variables is None:
            input_variables = [f"x{i}" for i in range(log2_length)]
        else:            
            if len(input_variables) != log2_length:
                raise ValueError("The length of input variables must be equal to the number of input variables.")
        SboxAnalyzer._write_truth_table(filename=input_file, 
                                        boolean_function=boolean_function, 
                                        input_output_variables=input_variables)        
        # print("Generateing and simplifying the MILP/SAT constraints ...")
        starting_time = time.time()
        SboxAnalyzer.simplify_by_espresso(input_file=input_file, 
                                          output_file=output_file, 
                                          mode=mode)
        elapsed_time = time.time() - starting_time
        # print("Time used to simplify the constraints: {:.2f} seconds".format(elapsed_time))        
        sat_clauses, milp_constraints, cp_constraints = SboxAnalyzer._parse_the_output_of_espresso(filename=output_file,
                                                                                   alphabet=input_variables)
        print("Number of constraints: {}".format(len(milp_constraints)))
        print("Variables: {}; msb: {}".format("||".join(input_variables), input_variables[0]))
        os.remove(input_file)
        os.remove(output_file)
        return (sat_clauses, milp_constraints, cp_constraints)
    
    @staticmethod
    def encode_set_of_binary_vectors(set_of_binary_vectors, input_variables=None, mode=6):
        """
        Encode the set of binary vectors into MILP/SAT constraints
        """

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1

        set_of_binary_vectors = list(set(set_of_binary_vectors))
        # check if each element of the set has the same length
        if len(set_of_binary_vectors) == 0:
            raise ValueError("The set of binary vectors is empty.")
        num_of_elements = len(set_of_binary_vectors)
        num_of_bits = len(set_of_binary_vectors[0])
        if not all([num_of_bits == len(v) for v in set_of_binary_vectors]):
            raise ValueError("All elements of the set must have the same length.")
        
        boolean_function = dict()
        for v in range(2**num_of_bits):
            key = tuple(map(int, list(bin(v)[2:].zfill(num_of_bits))))
            boolean_function[key] = int(key in set_of_binary_vectors) ^ reverse
        
        input_file = os.path.join(os.getcwd(), 'tt_' + hex(randint(0, 65536))[2:] + '.txt')
        output_file = os.path.join(os.getcwd(), 'stt_' + hex(randint(0, 65536))[2:] + '.txt')
        if input_variables is None:
            input_variables = [f"x{i}" for i in range(num_of_bits)]
        else:
            if len(input_variables) != num_of_bits:
                raise ValueError("The length of input variables must be equal to the length of the bit vectors.")
        SboxAnalyzer._write_truth_table(filename=input_file, 
                                        boolean_function=boolean_function, 
                                        input_output_variables=input_variables)
        # print("Generateing and simplifying the MILP/SAT constraints ...")
        starting_time = time.time()
        SboxAnalyzer.simplify_by_espresso(input_file=input_file, 
                                          output_file=output_file, 
                                          mode=mode)   
        elapsed_time = time.time() - starting_time
        # print("Time used to simplify the constraints: {:.2f} seconds".format(elapsed_time))
        sat_clauses, milp_constraints, cp_constraints = SboxAnalyzer._parse_the_output_of_espresso(filename=output_file,
                                                                                    alphabet=input_variables)
        print("Number of constraints: {}".format(len(milp_constraints)))
        print("Variables: {}; msb: {}".format("||".join(input_variables), input_variables[0]))
        os.remove(input_file)
        os.remove(output_file)
        return (sat_clauses, milp_constraints, cp_constraints)
    
    def encode_sbox(self, input_variables=None, mode=6):
        """
        Encode the S-box by a CNF or a set of integer inequalities
        """
        
        io_set = []
        for i in range(2**self.input_size()):
            io_set.append(tuple(list(map(int, list(bin(i)[2:].zfill(4)))) + list(map(int, list(bin(self(i))[2:].zfill(4))))))
        io_set = set(io_set)
        cnf, milp, cp = self.encode_set_of_binary_vectors(set_of_binary_vectors=io_set, input_variables=input_variables, mode=mode)
        return cnf, milp,cp

        
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
        
        self._ddt = self.difference_distribution_table()
        self._diff_spectrum = set([self._ddt[i][j] for i in range(2**self.input_size()) for j in range(2**self.output_size())]) - {0, 2**self.input_size()}
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
        Generate the *-DDT (or 0/1 DDT)
        Star DDT is a 2^m*2^n binary array describing the possibility of differential transitions through the S-box
        """
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        star_ddt = [[0 for _ in range(2**self.output_size())] for _ in range(2**self.input_size())]
        for dx in range(2**self.input_size()):
            for dy in range(2**self.output_size()):
                if self._ddt[dx][dy] != 0:
                    star_ddt[dx][dy] = 1 #reverse ^ 1
                else:
                    star_ddt[dx][dy] = 0 #reverse
        return star_ddt

    def _star_ddt_to_boolean_function(self, reverse=1, inverse=0):
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
                boolean_func[key] = self.star_ddt[dx][dy] ^ reverse ^ inverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse ^ inverse
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
                if self._ddt[dx][dy] in self._diff_spectrum:
                    p = tuple([int(i == self._ddt[dx][dy]) for i in self._diff_spectrum])
                    key = x + y + p
                    boolean_function[key] = 1
                elif self._ddt[dx][dy] == 2**self.input_size():
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
                if self._ddt[dx][dy] != 0:
                    w = int(abs(float(log(self._ddt[dx][dy]/(2**self.input_size()), 2))))
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
                if self._ddt[dx][dy] == p:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function

    def minimized_diff_constraints(self, mode=6, subtable=None, cryptosmt_compatible=False, input_variables=None, output_variables=None, cost_variables=None):
        """
        This method takes a given Boolean function and records its truth table in a file, 
        adhering to the ESPRESSO input format. It then invokes ESPRESSO to obtain a minimized 
        representation of the function. Following that, it interprets ESPRESSO's output and 
        converts the simplified representation into the language recognized by MILP and SAT solvers.
        """
        
        if self._data_required_for_differential_analysis is None:
            self._compute_data_for_differential_analysis()
        valid_values_for_subtable = ["star", "star_inverse", None] + list(self._diff_spectrum)
        if subtable not in valid_values_for_subtable:
            raise ValueError("Invalid value for subtable! subtable must be in {0}.".format(list(map(str, valid_values_for_subtable))))
    
        self.cryptosmt_compatible = cryptosmt_compatible
        self.ddt_subtable = subtable

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1
        
        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        else:            
            if len(input_variables) != self.input_size():
                raise ValueError("The length of input variables must be equal to the number of input variables.")
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        else:
            if len(output_variables) != self.output_size():
                raise ValueError("The length of output variables must be equal to the number of output variables.")
        if cost_variables is None:
            if cryptosmt_compatible:
                cost_variables = [f"p{i}" for i in range(self._max_diff_weights)]
            else:
                cost_variables = [f"p{i}" for i in range(self._len_diff_spectrum)]
        else:
            if cryptosmt_compatible:
                if len(cost_variables) != self._max_diff_weights:
                    raise ValueError("The length of cost variables mus be equal to self._max_diff_weights.")
            else:
                if len(cost_variables) != self._len_diff_spectrum:
                    raise ValueError("The length of cost variables mus be equal to self._len_diff_spectrum.")
        input_output_variables = input_variables + output_variables
        self.diff_objective = ""
        if subtable == "star":
            boolean_function = self._star_ddt_to_boolean_function(reverse=reverse)
        elif subtable == "star_inverse":
            boolean_function = self._star_ddt_to_boolean_function(reverse=reverse, inverse=1)
        elif subtable in self._diff_spectrum:
            boolean_function = self._pddt_to_booleanfunction(p=subtable, reverse=reverse)
        elif cryptosmt_compatible:
            input_output_variables += cost_variables
            boolean_function = self._ddt_to_cryptosmt_compatible_boolean_function(reverse=reverse)
            if self._max_diff_weights == 0:
                self.diff_objective = "0"
            else:
                self.diff_objective = cost_variables
            self.diff_objective = "\nWeight: {}".format(" + ".join(self.diff_objective))
        else:
            input_output_variables += cost_variables
            boolean_function = self._ddt_to_boolean_function(reverse=reverse)
            if self._max_diff_weights == 0:
                self.diff_objective = "0"
            else:
                self.diff_objective = [f"{self._diff_weights[i]:0.04f} {cost_variables[i]}" for i in range(self._len_diff_spectrum)]
            self.diff_objective = "\nWeight: {}".format(" + ".join(self.diff_objective))
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.diff_objective))
        return sat_clauses, milp_constraints, cp_constraints

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
        self._lat = self.linear_approximation_table(scale='correlation')        
        self._squared_lat = [[x**2 for x in y] for y in self._lat]
        self.lat_scaled_by_absolute_correlation = [[2**input_size * lxy for lxy in ly] for ly in self._lat]
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
        Generate the *-LAT (or 0/1 LAT)
        Star LAT is a 2^m*2^n binary array describing the possibility of linear transitions through the S-box
        """
        
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        star_lat = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for lx in range(2**self.input_size()):
            for ly in range(2**self.output_size()):
                if self._lat[lx][ly] != 0:
                    star_lat[lx][ly] = 1 #reverse ^ 1
                else:
                    star_lat[lx][ly] = 0 #reverse
        return star_lat
    
    def _star_lat_to_boolean_function(self, reverse=1, inverse=0):
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
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse ^ inverse
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

    def minimized_linear_constraints(self, mode=6, subtable=None, input_variables=None, output_variables=None, cost_variables=None):
        """
        This method takes a given Boolean function and records its truth table in a file, 
        adhering to the ESPRESSO input format. It then invokes ESPRESSO to obtain a minimized 
        representation of the function. Following that, it interprets ESPRESSO's output and 
        converts the simplified representation into the language recognized by MILP and SAT solvers.
        """
        
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        valid_values_for_subtable = ["star", "star_inverse", None] + list(self._squared_correlation_spectrum)
        if subtable not in valid_values_for_subtable:
            raise ValueError("Invalid value for subtable! subtable must be in {0}.".format(list(map(str, valid_values_for_subtable))))

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1

        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        else:
            if len(input_variables) != self.input_size():
                raise ValueError("The length of input variables must be equal to the number of input variables.")
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        else:
            if len(output_variables) != self.output_size():
                raise ValueError("The length of output variables must be equal to the number of output variables.")
        if cost_variables is None:
            cost_variables = [f"p{i}" for i in range(self._len_linear_weights)]
        else:
            if len(cost_variables) != self._len_linear_weights:
                raise ValueError("The length of cost variables must be equal to self._len_linear_weights.")
        input_output_variables = input_variables + output_variables
        self.linear_objective = ""
        if subtable == "star_inverse":
            boolean_function = self._star_lat_to_boolean_function(reverse=reverse, inverse=1)
        elif subtable == "star":
            boolean_function = self._star_lat_to_boolean_function(reverse=reverse)
        elif subtable in self._squared_correlation_spectrum:
            boolean_function = self._psqlat_to_booleanfunction(p=subtable, reverse=reverse)
        else:
            input_output_variables += cost_variables
            boolean_function = self._sqlat_to_boolean_function(reverse=reverse)
            if self._len_linear_weights != []:
                self.linear_objective = [f"{self._linear_weights[i]:0.04f} {cost_variables[i]}" for i in range(self._len_linear_weights)]
            else:
                self.linear_objective = "0"
            self.linear_objective = "\nWeight: {}".format(" + ".join(self.linear_objective))
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints, cp_constraints
    
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

    def minimized_integral_constraints(self, mode=6, input_variables=None, output_variables=None):
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
        
        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        else:
            if len(input_variables) != self.input_size():
                raise ValueError("The length of input variables must be equal to the number of input variables.")
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        else:
            if len(output_variables) != self.output_size():
                raise ValueError("The length of output variables must be equal to the number of output variables.")            
        
        self.integral_objective = ""
        input_output_variables = input_variables + output_variables
        boolean_function = self._mpt_to_boolean_function(reverse=reverse)
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.integral_objective))
        return sat_clauses, milp_constraints, cp_constraints
    
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
        
        if input_list == []:
            return []
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
                possible_output_diffs += [self.to_bits(x=dy, n=sbox_output_size) for dy in range(2**sbox_output_size) if self._ddt[input_diff][dy] != 0]
            truncated_output_diff = self.binvectors_to_truncated(possible_output_diffs)
            if any(x != self.unknown for x in truncated_output_diff):
                propagation_dictionary[truncated_input_diff] = truncated_output_diff
        return propagation_dictionary

    def encode_deterministic_linear_behavior(self, transpose=False):
        """
        Encodes the deterministic linear behavior of the S-box        
        """
        
        if self._data_required_for_linear_analysis is None:
            self._compute_data_for_linear_analysis()
        if transpose:
            temp_lat = self._lat.transpose()
        else:
            temp_lat = self._lat          
        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()
        propagation_dictionary = {}
        for truncated_input_mask in itertools.product(list(self.deterministic_mask.keys()), repeat=sbox_input_size):        
            binary_input_masks = self.truncated_to_binvectors(truncated_input_mask)
            possible_output_masks = []
            for input_mask in binary_input_masks:
                input_mask = int("".join(str(x) for x in input_mask), base=2)
                output_mask = [self.to_bits(x=ly, n=sbox_output_size) for ly in range(2**sbox_output_size) if temp_lat[input_mask][ly] != 0]
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
    
    def generate_cp_constraints(self, propagation_dictionary, input_variables=None, output_variables=None):
        """
        Generates the constraints for the CP model
        """

        inputs = list(propagation_dictionary.keys())
        outputs = list(propagation_dictionary.values())
        m = len(inputs[0])
        n = len(outputs[0])
        if input_variables is None:
            input_variables = [f"a{i}" for i in range(m)]
        else:            
            if len(input_variables) != m:
                raise ValueError(f"The size of input variables should be {m}")
            input_variables = input_variables
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(n)]
        else:
            if len(output_variables) != n:
                raise ValueError(f"The size of output variables should be {n}")
            output_variables = output_variables
        last_condition = " /\\ ".join([f"{output_variables[i]} = -1" for i in range(n)])
        constraints = ""
        for i, (input, output) in enumerate(propagation_dictionary.items()):
            input_str = [f"{input_variables[i]} == {input[i]}" for i in range(m)]
            input_str = " /\\ ".join(input_str)
            output_str = [f"{output_variables[i]} = {output[i]}" for i in range(n)]
            output_str = " /\\ ".join(output_str)
            if i == 0:
                constraints += f"if ({input_str}) then ({output_str})\n"
            else:
                constraints += f"elseif ({input_str}) then ({output_str})\n"
        constraints += f"else ({last_condition})\nendif"
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print(variables_mapping)
        return constraints
    
    def minimized_deterministic_diff_constraints(self, mode=6, input_variables=None, output_variables=None):
        """
        Generate and minimize the MILP/SAT constraints for the deterministic differential trails
        """

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1
        
        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()

        if input_variables is None:
            input_variables = [f"a{i}" for i in range(2*self.input_size())]
        elif len(input_variables) != 2*sbox_input_size:
            raise ValueError("The length of input variables must be equal to the number of input variables.")
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(2*self.output_size())]
        elif len(output_variables) != 2*sbox_output_size:
            raise ValueError("The length of output variables must be equal to the number of output variables.")        
        binary_variables = input_variables + output_variables
        propagation_dictionary = self.encode_deterministic_differential_behavior()
        informative_inputs = propagation_dictionary.keys()
        boolean_function = dict()
        function_size = 2*(sbox_input_size + sbox_output_size)
        for input_vector in itertools.product([self.unknown, self.zero, self.one], repeat=sbox_input_size):
            if input_vector in informative_inputs:
                output_vector = propagation_dictionary[input_vector]
                key = input_vector + tuple(output_vector)
                key = tuple(flatten([self.unknown_binary if x == self.unknown else self.zero_binary if x == self.zero else self.one_binary for x in key]))
            else:
                output_vector = [self.unknown]*sbox_output_size
                key = input_vector + tuple(output_vector)
                key = tuple(flatten([self.unknown_binary if x == self.unknown else self.zero_binary if x == self.zero else self.one_binary for x in key]))
            boolean_function[key] = 1
        
        if reverse == 1:
            for input_vector in itertools.product([self.zero, self.one ], repeat=function_size):
                boolean_function[input_vector] = boolean_function.get(input_vector, 0) ^ 1
    
        self._write_truth_table(filename=self.truth_table_filename,
                                boolean_function=boolean_function,
                                input_output_variables=binary_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=binary_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))  
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping_binary = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping_binary += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping_binary))
        return sat_clauses, milp_constraints, cp_constraints
    
    def minimized_deterministic_lin_constraints(self, mode=6, input_variables=None, output_variables=None):
        """
        Generate and minimize the MILP/SAT constraints for the deterministic linear trails
        """

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1        

        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()

        if input_variables is None:
            input_variables = [f"a{i}" for i in range(2*self.input_size())]
        elif len(input_variables) != 2*sbox_input_size:
            raise ValueError("The length of input variables must be equal to the number of input variables.")
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(2*self.output_size())]
        elif len(output_variables) != 2*sbox_output_size:
            raise ValueError("The length of output variables must be equal to the number of output variables.")
        binary_variables = input_variables + output_variables
        propagation_dictionary = self.encode_deterministic_linear_behavior()
        informative_inputs = propagation_dictionary.keys()
        boolean_function = dict()
        function_size = 2*(sbox_input_size + sbox_output_size)
        for input_vector in itertools.product([self.unknown, self.zero, self.one], repeat=sbox_input_size):
            if input_vector in informative_inputs:
                output_vector = propagation_dictionary[input_vector]
                key = input_vector + tuple(output_vector)
                key = tuple(flatten([self.unknown_binary if x == self.unknown else self.zero_binary if x == self.zero else self.one_binary for x in key]))
            else:
                output_vector = [self.unknown]*sbox_output_size
                key = input_vector + tuple(output_vector)
                key = tuple(flatten([self.unknown_binary if x == self.unknown else self.zero_binary if x == self.zero else self.one_binary for x in key]))
            boolean_function[key] = 1
        
        if reverse == 1:
            for input_vector in itertools.product([self.zero, self.one ], repeat=function_size):
                boolean_function[input_vector] = boolean_function.get(input_vector, 0) ^ 1
    
        self._write_truth_table(filename=self.truth_table_filename,
                                boolean_function=boolean_function,
                                input_output_variables=binary_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=binary_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # print("Time used to simplify the constraints: {:0.02} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping_binary = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping_binary += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping_binary))
        return sat_clauses, milp_constraints, cp_constraints

    def minimized_determinisitc_integral_forward_constraints(self, mode=6):
        """
        Generate and minimize the MILP/SAT constraints for the deterministic integral trails
        """

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1

        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()

        if input_variables is None:
            input_variables = [f"a{i}" for i in range(2*self.input_size())]
        elif len(input_variables) != 2*sbox_input_size:
            raise ValueError("The length of input variables must be equal to the number of input variables.")
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(2*self.output_size())]
        elif len(output_variables) != 2*sbox_output_size:
            raise ValueError("The length of output variables must be equal to the number of output variables.")
        binary_variables = input_variables + output_variables
        propagation_dictionary = self.encode_deterministic_integral_forward()
        informative_inputs = propagation_dictionary.keys()
        boolean_function = dict()
        function_size = 2*(sbox_input_size + sbox_output_size)
        for input_vector in itertools.product([self.unknown, self.zero, self.one], repeat=sbox_input_size):
            if input_vector in informative_inputs:
                output_vector = propagation_dictionary[input_vector]
                key = input_vector + tuple(output_vector)
                key = tuple(flatten([self.unknown_binary if x == self.unknown else self.zero_binary if x == self.zero else self.one_binary for x in key]))
            else:
                output_vector = [self.unknown]*sbox_output_size
                key = input_vector + tuple(output_vector)
                key = tuple(flatten([self.unknown_binary if x == self.unknown else self.zero_binary if x == self.zero else self.one_binary for x in key]))
            boolean_function[key] = 1
        
        if reverse == 1:
            for input_vector in itertools.product([self.zero, self.one ], repeat=function_size):
                boolean_function[input_vector] = boolean_function.get(input_vector, 0) ^ 1
    
        self._write_truth_table(filename=self.truth_table_filename,
                                boolean_function=boolean_function,
                                input_output_variables=binary_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=binary_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # print("Time used to simplify the constraints: {:0.02} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping_binary = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping_binary += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping_binary))
        return sat_clauses, milp_constraints, cp_constraints

    def minimized_determinisitc_integral_backward_constraints(self, mode=6):
        """
        Generate and minimize the MILP/SAT constraints for the deterministic integral trails
        """

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1

        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()

        if input_variables is None:
            input_variables = [f"a{i}" for i in range(2*self.input_size())]
        elif len(input_variables) != 2*sbox_input_size:
            raise ValueError("The length of input variables must be equal to the number of input variables.")
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(2*self.output_size())]
        elif len(output_variables) != 2*sbox_output_size:
            raise ValueError("The length of output variables must be equal to the number of output variables.")
        binary_variables = input_variables + output_variables
        propagation_dictionary = self.encode_deterministic_integral_backward()
        informative_inputs = propagation_dictionary.keys()
        boolean_function = dict()
        function_size = 2*(sbox_input_size + sbox_output_size)
        for input_vector in itertools.product([self.unknown, self.zero, self.one], repeat=sbox_output_size):
            if input_vector in informative_inputs:
                output_vector = propagation_dictionary[input_vector]
                key = input_vector + tuple(output_vector)
                key = tuple(flatten([self.unknown_binary if x == self.unknown else self.zero_binary if x == self.zero else self.one_binary for x in key]))
            else:
                output_vector = [self.unknown]*sbox_input_size
                key = input_vector + tuple(output_vector)
                key = tuple(flatten([self.unknown_binary if x == self.unknown else self.zero_binary if x == self.zero else self.one_binary for x in key]))
            boolean_function[key] = 1
        
        if reverse == 1:
            for input_vector in itertools.product([self.zero, self.one ], repeat=function_size):
                boolean_function[input_vector] = boolean_function.get(input_vector, 0) ^ 1
    
        self._write_truth_table(filename=self.truth_table_filename,
                                boolean_function=boolean_function,
                                input_output_variables=binary_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=binary_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # print("Time used to simplify the constraints: {:0.02} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping_binary = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping_binary += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping_binary))
        return sat_clauses, milp_constraints, cp_constraints    

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
        self._data_required_for_difflin_analysis = "Data for differential-linear analysis are computed and stored in memory."

    def differential_linear_connectivity_table(self):
        """
        Compute the differential-linear connectivity table
        """
        
        if self._data_required_for_difflin_analysis == None:
            self._compute_data_for_difflin_analysis()
        return self._dlct
    

    def get_dlct_spectrum(self):
        """
        Compute the set of different entries in DLCT
        """

        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        return self._dlct_spectrum
            
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
        ddt = self._ddt
        self._ddlct = [[0 for y in range(2**output_size)] for x in range(2**input_size)]
        for diff in range(2**input_size):
            for mask in range(2**output_size):
                for diff_middle in range(2**input_size):
                    self._ddlct[diff][mask] += ddt[diff][diff_middle] * dlct[diff_middle][mask]
        self._ddlct_spectrum = sorted(list(set(flatten(self._ddlct))))
        self._ddlct_weights = [abs(float(log(abs(x), 2))) for x in self._ddlct_spectrum if x != 0]
        self._len_ddlct_weights = len(self._ddlct_weights)   
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
        ddt = self._ddt
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

    def compute_star_dlct(self, reverse=1):
        """
        Generate the *-DLCT (or 0/1 DLCT)
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

    def _star_dlct_to_boolean_function(self, reverse=1, inverse=0):
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
                boolean_func[key] = self.star_dlct[dx][ly] ^ reverse ^ inverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse ^ inverse
        return boolean_func
    
    def minimized_differential_linear_constraints(self, mode=6, subtable=None, input_variables=None, output_variables=None):
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

        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        elif len(input_variables) != self.input_size():
            raise ValueError("The size of input variables should be {0}".format(self.input_size()))
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        elif len(output_variables) != self.output_size():
            raise ValueError("The size of output variables should be {0}".format(self.output_size()))
        input_output_variables = input_variables + output_variables
        self.linear_objective = ""
        if subtable == "star_inverse":
            boolean_function = self._star_dlct_to_boolean_function(reverse=reverse, inverse=1)
        elif subtable == "star":
            boolean_function = self._star_dlct_to_boolean_function(reverse=reverse)
        elif subtable in self._dlct_spectrum:
            boolean_function = self._dlct_to_booleanfunction(corr=subtable, reverse=reverse)
        else:
            boolean_function = self._dlct_to_booleanfunction(reverse=reverse)
            if self._len_dlct_weights != []:                
                self.linear_objective = ["{:0.04f} p{:d}".format(self._dlct_weights[i], i) for i in range(self._dlct_weights)]
            else:
                self.linear_objective = "0"
            self.linear_objective = "\nWeight: {}".format(" + ".join(self.linear_objective))
            input_output_variables.extend([f"p{i}" for i in range(self._len_linear_weights)])
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints, cp_constraints

    def compute_star_double_dlct(self, reverse=1):
        """
        Generate the *-DDLCT (or 0/1 DDLCT)
        Star DDLCT is a 2^m*2^n binary array describing the possibility of differential-linear transitions 
        through two consecutive S-boxes with key-addition or a linear layer in between. 
        """

        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        ddlct = self.double_differential_linear_connectivity_table()
        self.star_ddlct = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for dx in range(2**self.input_size()):
            for ly in range(2**self.output_size()):
                if ddlct[dx][ly] != 0:
                    self.star_ddlct[dx][ly] = reverse ^ 1
                else:
                    self.star_ddlct[dx][ly] = reverse
        return self.star_ddlct
    
    def _ddlct_to_booleanfunction(self, corr, reverse=1):
        """
        Convert a subtable of DDLCT into a Boolean function
        """
        
        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        
        boolean_function = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                if self._ddlct[dx][ly] == corr:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function
    
    def _star_ddlct_to_boolean_function(self, reverse=1, inverse=0):
        """
        Convert the *-DDLCT into a Boolean function
        """
 
        self.compute_star_double_dlct(reverse=reverse)
        boolean_func = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.star_ddlct[dx][ly] ^ reverse ^ inverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse ^ inverse
        return boolean_func
    
    def minimized_double_differential_linear_constraints(self, mode=6, subtable=None, input_variables=None, output_variables=None):
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
        
        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        elif len(input_variables) != self.input_size():
            raise ValueError("The size of input variables should be {0}".format(self.input_size()))
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        elif len(output_variables) != self.output_size():
            raise ValueError("The size of output variables should be {0}".format(self.output_size()))
        input_output_variables = input_variables + output_variables
        
        self.linear_objective = ""
        if subtable == "star_inverse":
            boolean_function = self._star_ddlct_to_boolean_function(reverse=reverse, inverse=1)
        elif subtable == "star":
            boolean_function = self._star_ddlct_to_boolean_function(reverse=reverse)
        elif subtable in self._dlct_spectrum:
            boolean_function = self._ddlct_to_booleanfunction(corr=subtable, reverse=reverse)
        else:
            boolean_function = self._ddlct_to_booleanfunction(reverse=reverse)
            if self._len_dlct_weights != []:                
                self.linear_objective = ["{:0.04f} p{:d}".format(self._dlct_weights[i], i) for i in range(self._len_dlct_weights)]
            else:
                self.linear_objective = "0"
            self.linear_objective = "\nWeight: {}".format(" + ".join(self.linear_objective))
            input_output_variables.extend([f"p{i}" for i in range(self._len_dlct_weights)])
        self._write_truth_table(filename=self.truth_table_filename,
                                boolean_function=boolean_function,
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints, cp_constraints
    
    def triple_differential_linear_connectivity_table(self):
        """
        Compute the triple differential-linear connectivity table
        """

        input_size = self.input_size()
        output_size = self.output_size()
        if self._data_required_for_difflin_analysis == None:
            self._compute_data_for_difflin_analysis()
        if self._data_required_for_differential_analysis == None:
            self._compute_data_for_differential_analysis()
        if self._data_required_for_linear_analysis == None:
            self._compute_data_for_linear_analysis()
        ddlct = self.double_differential_linear_connectivity_table()        
        ddt = self._ddt
        self._tdlct = [[0 for _ in range(2**output_size)] for _ in range(2**input_size)]
        for dx in range(2**input_size):
            for ly in range(2**output_size):
                for dm in range(2**output_size):
                    self._tdlct[dx][ly] += ddt[dx][dm] * ddlct[dm][ly]
        return self._tdlct
    
    def _tdlct_to_booleanfunction(self, corr, reverse=1):
        """
        Convert a subtable of TDLCT into a Boolean function
        """
        
        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        tdlct = self.triple_differential_linear_connectivity_table()
        boolean_function = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                if tdlct[dx][ly] == corr:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function

    def compute_star_triple_dlct(self, reverse=1):
        """
        Generate the *-TDLCT (or 0/1 TDLCT)
        Star TDLCT is a 2^m*2^n binary array describing the possibility of differential-linear transitions 
        through three consecutive S-boxes with key-addition or a linear layer in between. 
        """

        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        tdlct = self.triple_differential_linear_connectivity_table()
        self.star_tdlct = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for dx in range(2**self.input_size()):
            for ly in range(2**self.output_size()):
                if tdlct[dx][ly] != 0:
                    self.star_tdlct[dx][ly] = reverse ^ 1
                else:
                    self.star_tdlct[dx][ly] = reverse
        return self.star_tdlct
        
    def _star_tdlct_to_boolean_function(self, reverse=1, inverse=0):
        """
        Convert the *-TDLCT into a Boolean function
        """
 
        self.compute_star_triple_dlct(reverse=reverse)
        boolean_func = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.star_tdlct[dx][ly] ^ reverse ^ inverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse ^ inverse
        return boolean_func
    
    def minimized_triple_differential_linear_constraints(self, mode=6, subtable=None, input_variables=None, output_variables=None):
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
        
        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        elif len(input_variables) != self.input_size():
            raise ValueError("The size of input variables should be {0}".format(self.input_size()))
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        elif len(output_variables) != self.output_size():
            raise ValueError("The size of output variables should be {0}".format(self.output_size()))
        input_output_variables = input_variables + output_variables
        
        self.linear_objective = ""
        if subtable == "star_inverse":
            boolean_function = self._star_tdlct_to_boolean_function(reverse=reverse, inverse=1)
        elif subtable == "star":
            boolean_function = self._star_tdlct_to_boolean_function(reverse=reverse)
        elif subtable in self._dlct_spectrum:
            boolean_function = self._tdlct_to_booleanfunction(corr=subtable, reverse=reverse)
        else:
            boolean_function = self._tdlct_to_booleanfunction(reverse=reverse)
            if self._len_dlct_weights != []:                
                self.linear_objective = ["{:0.04f} p{:d}".format(self._dlct_weights[i], i) for i in range(self._len_dlct_weights)]
            else:
                self.linear_objective = "0"
            self.linear_objective
            input_output_variables.extend([f"p{i}" for i in range(self._len_dlct_weights)])
        self._write_truth_table(filename=self.truth_table_filename,
                                boolean_function=boolean_function,
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints, cp_constraints
    
    def quadruple_differential_linear_connectivity_table(self):
        """
        Compute the quadruple differential-linear connectivity table (4-DLCT)
        """

        input_size = self.input_size()
        output_size = self.output_size()
        if self._data_required_for_difflin_analysis == None:
            self._compute_data_for_difflin_analysis()
        if self._data_required_for_differential_analysis == None:
            self._compute_data_for_differential_analysis()
        if self._data_required_for_linear_analysis == None:
            self._compute_data_for_linear_analysis()
        tdlct = self.triple_differential_linear_connectivity_table()
        ddt = self._ddt
        self._qdlct = [[0 for _ in range(2**output_size)] for _ in range(2**input_size)]
        for dx in range(2**input_size):
            for ly in range(2**output_size):
                for dm in range(2**output_size):
                    self._qdlct[dx][ly] += ddt[dx][dm] * tdlct[dm][ly]
        return self._qdlct

    def _qdlct_to_booleanfunction(self, corr, reverse=1):
        """
        Convert a subtable of 4-DLCT into a Boolean function
        """
        
        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        qdlct = self.quadruple_differential_linear_connectivity_table()
        boolean_function = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                if qdlct[dx][ly] == corr:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function
    
    def compute_star_quadruple_dlct(self, reverse=1):
        """
        Generate the *-4-DLCT (or 0/1 QDLCT)
        Star QDLCT is a 2^m*2^n binary array describing the possibility of differential-linear transitions 
        through four consecutive S-boxes with key-addition or a linear layer in between. 
        """

        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        quaddlct = self.quadruple_differential_linear_connectivity_table()
        self.star_quaddlct = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for dx in range(2**self.input_size()):
            for ly in range(2**self.output_size()):
                if quaddlct[dx][ly] != 0:
                    self.star_quaddlct[dx][ly] = reverse ^ 1
                else:
                    self.star_quaddlct[dx][ly] = reverse
        return self.star_quaddlct
    
    def _star_quaddlct_to_boolean_function(self, reverse=1, inverse=0):
        """
        Convert the *-4-DLCT into a Boolean function
        """
 
        self.compute_star_quadruple_dlct(reverse=reverse)
        boolean_func = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.star_quaddlct[dx][ly] ^ reverse ^ inverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse ^ inverse
        return boolean_func
    
    def minimized_quadruple_differential_linear_constraints(self, mode=6, subtable=None, input_variables=None, output_variables=None):
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
        
        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        elif len(input_variables) != self.input_size():
            raise ValueError("The size of input variables should be {0}".format(self.input_size()))
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        elif len(output_variables) != self.output_size():
            raise ValueError("The size of output variables should be {0}".format(self.output_size()))
        input_output_variables = input_variables + output_variables
        
        self.linear_objective = ""
        if subtable == "star_inverse":
            boolean_function = self._star_quaddlct_to_boolean_function(reverse=reverse, inverse=1)
        elif subtable == "star":
            boolean_function = self._star_quaddlct_to_boolean_function(reverse=reverse)
        elif subtable in self._dlct_spectrum:
            boolean_function = self._qdlct_to_booleanfunction(corr=subtable, reverse=reverse)
        else:
            boolean_function = self._qdlct_to_booleanfunction(reverse=reverse)
            if self._len_dlct_weights != []:                
                self.linear_objective = ["{:0.04f} p{:d}".format(self._dlct_weights[i], i) for i in range(self._len_dlct_weights)]
            else:
                self.linear_objective = "0"
            self.linear
            input_output_variables.extend([f"p{i}" for i in range(self._len_dlct_weights)])
        self._write_truth_table(filename=self.truth_table_filename,
                                boolean_function=boolean_function,
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints, cp_constraints
    
    def quintuple_differential_linear_connectivity_table(self):
        """
        Compute the quintuple differential-linear connectivity table (5-DLCT)
        """

        input_size = self.input_size()
        output_size = self.output_size()
        if self._data_required_for_difflin_analysis == None:
            self._compute_data_for_difflin_analysis()
        if self._data_required_for_differential_analysis == None:
            self._compute_data_for_differential_analysis()
        if self._data_required_for_linear_analysis == None:
            self._compute_data_for_linear_analysis()
        quaddlct = self.quadruple_differential_linear_connectivity_table()
        ddt = self._ddt
        self._quindlct = [[0 for _ in range(2**output_size)] for _ in range(2**input_size)]
        for dx in range(2**input_size):
            for ly in range(2**output_size):
                for dm in range(2**output_size):
                    self._quindlct[dx][ly] += ddt[dx][dm] * quaddlct[dm][ly]
        return self._quindlct
    
    def _quindlct_to_booleanfunction(self, corr, reverse=1):
        """
        Convert a subtable of 5-DLCT into a Boolean function
        """
        
        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        quindlct = self.quintuple_differential_linear_connectivity_table()
        boolean_function = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                if quindlct[dx][ly] == corr:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function
    
    def compute_star_quintuple_dlct(self, reverse=1):
        """
        Generate the *-5-DLCT (or 0/1 QDLCT)
        Star QDLCT is a 2^m*2^n binary array describing the possibility of differential-linear transitions 
        through five consecutive S-boxes with key-addition or a linear layer in between. 
        """

        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        quindlct = self.quintuple_differential_linear_connectivity_table()
        self.star_quindlct = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for dx in range(2**self.input_size()):
            for ly in range(2**self.output_size()):
                if quindlct[dx][ly] != 0:
                    self.star_quindlct[dx][ly] = reverse ^ 1
                else:
                    self.star_quindlct[dx][ly] = reverse
        return self.star_quindlct

    def _star_quindlct_to_boolean_function(self, reverse=1, inverse=0):
        """
        Convert the *-5-DLCT into a Boolean function
        """
 
        self.compute_star_quintuple_dlct(reverse=reverse)
        boolean_func = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for ly in range(2**self.output_size()):
                y = tuple(map(int, list(bin(ly)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.star_quindlct[dx][ly] ^ reverse ^ inverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse ^ inverse
        return boolean_func
    
    def minimized_quintuple_differential_linear_constraints(self, mode=6, subtable=None, input_variables=None, output_variables=None):
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
        
        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        elif len(input_variables) != self.input_size():
            raise ValueError("The size of input variables should be {0}".format(self.input_size()))
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        elif len(output_variables) != self.output_size():
            raise ValueError("The size of output variables should be {0}".format(self.output_size()))
        input_output_variables = input_variables + output_variables
        
        self.linear_objective = ""
        if subtable == "star_inverse":
            boolean_function = self._star_quindlct_to_boolean_function(reverse=reverse, inverse=1)
        elif subtable == "star":
            boolean_function = self._star_quindlct_to_boolean_function(reverse=reverse)
        elif subtable in self._dlct_spectrum:
            boolean_function = self._quindlct_to_booleanfunction(corr=subtable, reverse=reverse)
        else:
            boolean_function = self._quindlct_to_booleanfunction(reverse=reverse)
            if self._len_dlct_weights != []:                
                self.linear_objective = ["{:0.04f} p{:d}".format(self._dlct_weights[i], i) for i in range(self._len_dlct_weights)]
            else:
                self.linear_objective = "0"
            self.linear_objective
            input_output_variables.extend([f"p{i}" for i in range(self._len_dlct_weights)])
        self._write_truth_table(filename=self.truth_table_filename,
                                boolean_function=boolean_function,
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints, cp_constraints

    def encode_deterministic_difflin_behavior(self, inverse=0):
        """
        Encodes the deterministic differential-linear behavior of the S-box
        """

        if self._data_required_for_difflin_analysis is None:
            self._compute_data_for_difflin_analysis()
        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()
        propagation_dictionary = {}
        for truncated_input_diff in itertools.product(list(self.deterministic_mask.keys()), repeat=sbox_input_size):        
            binary_input_diffs = self.truncated_to_binvectors(truncated_input_diff)
            possible_output_masks = []
            for input_diff in binary_input_diffs:
                input_diff = int("".join(str(x) for x in input_diff), base=2)
                possible_output_masks += [self.to_bits(x=dy, n=sbox_output_size) for dy in range(2**sbox_output_size) if ((self._dlct[input_diff][dy] != 0) ^ inverse)]
            truncated_output_mask = self.binvectors_to_truncated(possible_output_masks)
            if any(x != self.unknown for x in truncated_output_mask):
                propagation_dictionary[truncated_input_diff] = truncated_output_mask
        return propagation_dictionary

    ###############################################################################################################
    ###############################################################################################################
    ###############################################################################################################
    #  ____                                                                _                   _              _      
    # | __ )   ___    ___   _ __ ___    ___  _ __  __ _  _ __    __ _     / \    _ __    __ _ | | _   _  ___ (_) ___ 
    # |  _ \  / _ \  / _ \ | '_ ` _ \  / _ \| '__|/ _` || '_ \  / _` |   / _ \  | '_ \  / _` || || | | |/ __|| |/ __|
    # | |_) || (_) || (_) || | | | | ||  __/| |  | (_| || | | || (_| |  / ___ \ | | | || (_| || || |_| |\__ \| |\__ \
    # |____/  \___/  \___/ |_| |_| |_| \___||_|   \__,_||_| |_| \__, | /_/   \_\|_| |_| \__,_||_| \__, ||___/|_||___/
    #                                                           |___/                             |___/              
    # Boomerang Analysis
    ###############################################################################################################    
    def _compute_data_for_boomerang_analysis(self):
        """
        Compute the data required for boomerang analysis
        """

        self._bct = [[x for x in row] for row in self.boomerang_connectivity_table()]
        self._bct_spectrum = sorted(list(set(flatten(self._bct))))
        self._bct_weights = [abs(float(log(abs(x), 2))) for x in self._bct_spectrum if x != 0]
        self._len_bct_weights = len(self._bct_weights)
        self._data_required_for_boomerang_analysis = "Data for boomerang analysis are computed and stored in memory."

    def get_bct_spectrum(self):
        """
        Compute the set of different entries in BCT
        """

        if self._data_required_for_boomerang_analysis is None:
            self._compute_data_for_boomerang_analysis()
        return self._bct_spectrum
    
    def compute_star_bct(self, reverse=1):
        """
        Generate the *-BCT (or 0/1 BCT)
        Star BCT is a 2^m*2^n binary array describing the possibility of boomerang transitions through the S-box
        """

        if self._data_required_for_boomerang_analysis is None:
            self._compute_data_for_boomerang_analysis()
        self.star_bct = [[0 for i in range(2**self.output_size())] for j in range(2**self.input_size())]
        for dx in range(2**self.input_size()):
            for dy in range(2**self.output_size()):
                if self._bct[dx][dy] != 0:
                    self.star_bct[dx][dy] = reverse ^ 1
                else:
                    self.star_bct[dx][dy] = reverse
        return self.star_bct

    def _bct_to_booleanfunction(self, pr, reverse=1):
        """
        Convert a subtable of BCT into a Boolean function
        """
        if self._data_required_for_boomerang_analysis is None:
            self._compute_data_for_boomerang_analysis()
        boolean_function = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for dy in range(2**self.output_size()):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.output_size()))))
                # Specifying 0 points is not necessary at the input of ESPRESSO
                key = x + y
                if self._bct[dx][dy] == pr:
                    boolean_function[key] = reverse ^ 1
                else:
                    boolean_function[key] = reverse
        return boolean_function

    def _star_bct_to_boolean_function(self, reverse=1, inverse=0):
        """
        Convert the star-BCT into a Boolean function
        """

        self.compute_star_bct(reverse=reverse)
        boolean_func = dict()
        for dx in range(2**self.input_size()):
            x = tuple(map(int, list(bin(dx)[2:].zfill(self.input_size()))))
            for dy in range(2**self.output_size()):
                y = tuple(map(int, list(bin(dy)[2:].zfill(self.output_size()))))
                key = x + y
                boolean_func[key] = self.star_bct[dx][dy] ^ reverse ^ inverse
        boolean_func[tuple([0]*self.input_size() + [0]*self.output_size())] = 1 ^ reverse ^ inverse
        return boolean_func
    
    def minimized_boomerang_constraints(self, mode=6, subtable=None, input_variables=None, output_variables=None):
        """
        This method takes a given Boolean function and records its truth table in a file, 
        adhering to the ESPRESSO input format. It then invokes ESPRESSO to obtain a minimized 
        representation of the function. Following that, it interprets ESPRESSO's output and 
        converts the simplified representation into the language recognized by MILP and SAT solvers.
        """
        
        if self._data_required_for_boomerang_analysis is None:
            self._compute_data_for_boomerang_analysis()
        valid_values_for_subtable = ["star", "star_inverse", None] + list(self._bct_spectrum)
        if subtable not in valid_values_for_subtable:
            raise ValueError("Invalid value for subtable! subtable must be in {0}.".format(list(map(str, valid_values_for_subtable))))

        if mode in [1, 3, 5, 6, 7]:
            reverse = 0
        else:
            reverse = 1

        if input_variables is None:
            input_variables = [f"a{i}" for i in range(self.input_size())]
        elif len(input_variables) != self.input_size():
            raise ValueError("The size of input variables should be {0}".format(self.input_size()))
        if output_variables is None:
            output_variables = [f"b{i}" for i in range(self.output_size())]
        elif len(output_variables) != self.output_size():
            raise ValueError("The size of output variables should be {0}".format(self.output_size()))
        input_output_variables = input_variables + output_variables
        self.linear_objective = ""
        if subtable == "star_inverse":
            boolean_function = self._star_bct_to_boolean_function(reverse=reverse, inverse=1)
        elif subtable == "star":
            boolean_function = self._star_bct_to_boolean_function(reverse=reverse)
        elif subtable in self._bct_spectrum:
            boolean_function = self._bct_to_booleanfunction(corr=subtable, reverse=reverse)
        else:
            boolean_function = self._bct_to_booleanfunction(reverse=reverse)
            if self._len_bct_weights != []:                
                self.linear_objective = ["{:0.04f} p{:d}".format(self._bct_weights[i], i) for i in range(self._bct_weights)]
            else:
                self.linear_objective = "0"
            self.linear_objective = "\nWeight: {}".format(" + ".join(self.linear_objective))
            input_output_variables.extend([f"p{i}" for i in range(self._len_linear_weights)])
        self._write_truth_table(filename=self.truth_table_filename, 
                                boolean_function=boolean_function, 
                                input_output_variables=input_output_variables)
        starting_time = time.time()
        # print("Simplifying the MILP/SAT constraints ...")
        self.simplify_by_espresso(input_file=self.truth_table_filename, output_file=self.simplified_truth_table_filename, mode=mode)            
        elapsed_time = time.time() - starting_time
        sat_clauses, milp_constraints, cp_constraints = self._parse_the_output_of_espresso(filename=self.simplified_truth_table_filename, alphabet=input_output_variables)
        os.remove(self.simplified_truth_table_filename)
        os.remove(self.truth_table_filename)        
        # print("Time used to simplify the constraints: {:0.02f} seconds".format(elapsed_time))
        print("Number of constraints: {0}".format(len(milp_constraints)))
        variables_mapping = "Input: {0}; msb: {1}".format("||".join(input_variables), input_variables[0])
        variables_mapping += "\nOutput: {0}; msb: {1}".format("||".join(output_variables), output_variables[0])
        print("{}".format(variables_mapping + self.linear_objective))
        return sat_clauses, milp_constraints, cp_constraints

    def encode_deterministic_boomerang_behavior(self, inverse=0):
        """
        Encodes the deterministic boomerang behavior of the S-box
        """

        if self._data_required_for_boomerang_analysis is None:
            self._compute_data_for_boomerang_analysis()                
        sbox_input_size = self.input_size()
        sbox_output_size = self.output_size()
        propagation_dictionary = {}
        for truncated_input_diff in itertools.product(list(self.deterministic_mask.keys()), repeat=sbox_input_size):        
            binary_input_diffs = self.truncated_to_binvectors(truncated_input_diff)
            possible_output_diffs = []
            for input_diff in binary_input_diffs:
                input_diff = int("".join(str(x) for x in input_diff), base=2)
                possible_output_diffs += [self.to_bits(x=dy, n=sbox_output_size) for dy in range(2**sbox_output_size) if ((self._bct[input_diff][dy] != 0) ^ inverse)]
            truncated_output_diff = self.binvectors_to_truncated(possible_output_diffs)
            if any(x != self.unknown for x in truncated_output_diff):
                propagation_dictionary[truncated_input_diff] = truncated_output_diff
        return propagation_dictionary
    
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
    # from sage.crypto.sboxes import SBox
    # from sage.crypto.sboxes import CLEFIA_S0 as sb
    # sb = SBox([sb(i ^ 0x8) ^ sb(i) for i in range(16)])    
    # sb = SBox([0, 1, 1, 0])
    # Orthros's Sbox    
    # sb = SBox([1, 0, 2, 4, 3, 8, 6, 0xd, 9, 0xa, 0xb, 0xe, 0xf, 0xc, 7, 5]) 
    # Spongent's Sbox
    # sb = SBox([0xe, 0xd, 0xb, 0, 2, 1, 4, 0xf, 7, 0xa, 8, 5, 9, 0xc, 3, 6]) 
    # SPEEDY's S-box
    # sb = SBox([0x08, 0x00, 0x09, 0x03, 0x38, 0x10, 0x29, 0x13, 0x0C, 0x0D, 0x04, 0x07, 0x30, 0x01, 0x20, 0x23, 0x1A, 0x12, 0x18, 0x32, 0x3E, 0x16, 0x2C, 0x36, 0x1C, 0x1D, 0x14, 0x37, 0x34, 0x05, 0x24, 0x27, 0x02, 0x06, 0x0B, 0x0F, 0x33, 0x17, 0x21, 0x15, 0x0A, 0x1B, 0x0E, 0x1F, 0x31, 0x11, 0x25, 0x35, 0x22, 0x26, 0x2A, 0x2E, 0x3A, 0x1E, 0x28, 0x3C, 0x2B, 0x3B, 0x2F, 0x3F, 0x39, 0x19, 0x2D, 0x3D])    
    # sb = sb.inverse()
    # sa = SboxAnalyzer(sb)
    # print(sa)
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
    # lin_spec = sa.get_squared_correlation_spectrum()
    # stn = 5
    # print(lin_spec)
    # print("\nencode {0}-LAT".format(lin_spec[stn]))
    # print(float(log(lin_spec[stn], 2)))
    # cnf, milp = sa.minimized_linear_constraints(subtable=lin_spec[stn], mode=7)

    # with open("milp.lp", "w") as fileobj:
        # fileobj.write("\nsubject to\n")
        # fileobj.write("\n".join(milp))
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
    # reset()
    # from sboxanalyzer import SboxAnalyzer
    # from sage.crypto.sboxes import CRAFT as sb
    # sa = SboxAnalyzer(sb)
    # print("\nEncoding deterministic differential behavior ...\n")
    # diff_propagation = sa.encode_deterministic_differential_behavior()
    # print("\nDifferential behavior:")
    # pretty_print(diff_propagation)
    # cp_constraints = sa.generate_cp_constraints(diff_propagation)
    # print("\nencode deterministic differential behavior")
    # print(cp_constraints)
    # cnf, milp = sa.minimized_deterministic_diff_constraints()
    # with open("milp.lp", "w") as fileobj:
    #     fileobj.write("\nsubject to\n")
    #     fileobj.write("\n".join(milp))
    #     fileobj.write("\nbin\n")
    #     for var in sa.variables_binary:
    #         fileobj.write(f"{var}\n")            
    #     fileobj.write("end")
    # pretty_print(milp)


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
    # reset()
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

    # # example of encoding a Boolean function
    # reset()
    # print("\nEncoding a Boolean function")
    # truth_table = [0, 0, 0, 0, 0, 0, 1, 1]
    # input_variables = ["a0", "a1", "a2"]
    # from sboxanalyzer import SboxAnalyzer as sa
    # cnf, milp = sa.encode_boolean_function(truth_table=truth_table, input_variables=input_variables, mode=2)
    # with open("milp.lp", "w") as fileobj:
    #     fileobj.write("\nsubject to\n")
    #     fileobj.write("\n".join(milp))
    #     fileobj.write("\nbin\n")
    #     for var in input_variables:
    #         fileobj.write(f"{var}\n")
    #     fileobj.write("end")
    # pretty_print(milp)

    # # example of encoding a set of binary vectors
    # reset()
    # print("\nEncoding a set of binary vectors")
    # binary_vectors = [(1, 1, 1, 1, 0, 0, 0, 1), (1, 1, 1, 1, 1, 0, 0, 0), (1, 1, 1, 1, 0, 0, 1, 1), (1, 1, 1, 1, 0, 0, 1, 0)]
    # input_variables = ["a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7"]
    # from sboxanalyzer import SboxAnalyzer as sa
    # cnf, milp = sa.encode_set_of_binary_vectors(set_of_binary_vectors=binary_vectors, input_variables=input_variables, mode=6)
    # with open("milp.lp", "w") as fileobj:
    #     fileobj.write("\nsubject to\n")
    #     fileobj.write("\n".join(milp))
    #     fileobj.write("\nbin\n")
    #     for var in input_variables:
    #         fileobj.write(f"{var}\n")
    #     fileobj.write("end")
    # pretty_print(milp)

    # # compute the ddlct of GIFT-64 
    # reset()
    # from sage.crypto.sboxes import GIFT as sb
    # from sboxanalyzer import SboxAnalyzer
    # sa = SboxAnalyzer(sb)
    # dlct = sa.double_differential_linear_connectivity_table()
    # ddt = sa.difference_distribution_table()
    # group_maps = [12, 1, 6, 11, 8, 13, 2, 7, 4, 9, 14, 3, 0, 5, 10, 15]
    # input_diffs = [0x1 << i for i in range(16)]
    # output_masks = [0x1 << i for i in range(16)]
    # ddlct = dict()
    # for di in input_diffs:
    #     for lo in output_masks:
    #         ddlct[(di, lo)] = 0
    # for di in input_diffs:
    #     # split di into 4 chunkes of 4 bits each 
    #     di0 = (di >> 12) & 0xf
    #     di1 = (di >> 8) & 0xf
    #     di2 = (di >> 4) & 0xf
    #     di3 = di & 0xf
    #     for lo in output_masks:
    #         # split li into 4 chunkes of 4 bits each
    #         lo0 = (lo >> 12) & 0xf
    #         lo1 = (lo >> 8) & 0xf
    #         lo2 = (lo >> 4) & 0xf
    #         lo3 = lo & 0xf
    #         for dm in range(2**16):
    #             dm0 = (dm >> 12) & 0xf
    #             dm1 = (dm >> 8) & 0xf
    #             dm2 = (dm >> 4) & 0xf
    #             dm3 = dm & 0xf
    #             if ddt[di0][dm0] * ddt[di1][dm1] * ddt[di2][dm2] * ddt[di3][dm3] > 0:                    
    #                 dm = list(bin(dm)[2:].zfill(16))
    #                 # apply the group map such that position i in dm is mapped to group_maps[i]
    #                 dm = [dm[group_maps[i]] for i in range(16)]
    #                 dm = int("".join(dm), 2)                
    #                 dm0n = (dm >> 12) & 0xf
    #                 dm1n = (dm >> 8) & 0xf
    #                 dm2n = (dm >> 4) & 0xf
    #                 dm3n = dm & 0xf
    #                 if dlct[dm0n][lo0]*dlct[dm1n][lo1]*dlct[dm2n][lo2]*dlct[dm3n][lo3] > 0:
    #                     ddlct[(di, lo)] = ddlct.get((di, lo), 0) +\
    #                                       ddt[di0][dm0]*dlct[dm0n][lo0]*\
    #                                       ddt[di1][dm1]*dlct[dm1n][lo1]*\
    #                                       ddt[di2][dm2]*dlct[dm2n][lo2]*\
    #                                       ddt[di3][dm3]*dlct[dm3n][lo3]
    # for di in input_diffs:
    #     cnt = 0
    #     for lo in output_masks:
    #         if ddlct[(di, lo)] != 0:
    #             cnt += 1
    #         print(f"ddlct[{di}][{lo}] = {ddlct[(di, lo)]}")
    #     print(f"Number of non-zero entries for input difference {di}: {cnt}")

    # test sbox encoder
    reset()
    from sage.crypto.sboxes import PRESENT as sb
    from sboxanalyzer import SboxAnalyzer as SA
    sa = SA(sb)
    cnf, milp, cp = sa.encode_sbox()
    ## compare it with cnf generator of SageMath
    cnf_sage = sa.cnf()
    print("Length of CNF generated by S-box Analzyer: {:d}".format(len(cnf.split("&"))))
    print("Length of CNF generated by SageMath      : {:d}".format(len(cnf_sage)))


   # comment 
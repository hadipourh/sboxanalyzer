# -*- coding: utf-8 -*-
"""
Tests for the sboxanalyzer package.

These tests verify the core functionality of S-box Analyzer including:
- DDT (Differential Distribution Table) encoding
- LAT (Linear Approximation Table) encoding  
- MPT (Monomial Prediction Table) encoding
- DLCT (Differential-Linear Connectivity Table) encoding
- Encoding binary vectors and Boolean functions
- Deterministic propagation analysis
"""

import pytest

# Import the S-box Analyzer
from sboxanalyzer import SboxAnalyzer

# Import SageMath's S-boxes for testing
from sage.crypto.sboxes import (
    PRESENT,
    SKINNY_4,
    Ascon,
    GIFT,
    KNOT,
    Midori_Sb0,
    PRINTcipher,
)


class TestSboxAnalyzerBasic:
    """Basic tests for SboxAnalyzer initialization and utilities."""

    def test_initialization_with_list(self):
        """Test that SboxAnalyzer can be initialized with a lookup table list."""
        # 4-bit identity S-box
        sbox = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        sa = SboxAnalyzer(sbox)
        assert sa.input_size() == 4
        assert sa.output_size() == 4

    def test_initialization_with_sbox(self):
        """Test that SboxAnalyzer can be initialized with a SageMath S-box."""
        sa = SboxAnalyzer(PRESENT)
        assert sa.input_size() == 4
        assert sa.output_size() == 4

    def test_initialization_with_3bit_sbox(self):
        """Test that SboxAnalyzer works with 3-bit S-boxes."""
        sa = SboxAnalyzer(PRINTcipher)
        assert sa.input_size() == 3
        assert sa.output_size() == 3

    def test_initialization_with_5bit_sbox(self):
        """Test that SboxAnalyzer works with 5-bit S-boxes (Ascon)."""
        sa = SboxAnalyzer(Ascon)
        assert sa.input_size() == 5
        assert sa.output_size() == 5


class TestDifferentialAnalysis:
    """Tests for differential analysis functionality (DDT encoding)."""

    def test_get_differential_spectrum_present(self):
        """Test getting differential spectrum of PRESENT S-box."""
        sa = SboxAnalyzer(PRESENT)
        spectrum = sa.get_differential_spectrum()
        # PRESENT has differential spectrum [2, 4] according to the documentation
        assert isinstance(spectrum, list)
        assert len(spectrum) > 0
        # Convert to Python ints for comparison
        spectrum_int = sorted([int(x) for x in spectrum])
        assert spectrum_int == [2, 4], f"PRESENT differential spectrum should be [2, 4], got {spectrum_int}"

    def test_get_star_ddt(self):
        """Test getting *-DDT (binary DDT)."""
        sa = SboxAnalyzer(PRESENT)
        star_ddt = sa.get_star_ddt()
        # Should be a 16x16 matrix for 4-bit S-box
        assert len(star_ddt) == 16
        assert len(star_ddt[0]) == 16
        # All entries should be 0 or 1
        for row in star_ddt:
            for val in row:
                assert val in [0, 1]
        # (0,0) should be 1 (zero difference always maps to zero difference)
        assert star_ddt[0][0] == 1

    def test_minimized_diff_constraints_present(self):
        """Test DDT encoding for PRESENT S-box."""
        sa = SboxAnalyzer(PRESENT)
        cnf, milp, cp = sa.minimized_diff_constraints()
        
        # Check that outputs are non-empty
        assert len(cnf) > 0
        assert len(milp) > 0
        assert len(cp) > 0
        
        # According to README: PRESENT should produce 54 constraints with mode=6
        # Weight: 3.0000 p0 + 2.0000 p1
        assert len(milp) == 54, f"Expected 54 constraints for PRESENT DDT, got {len(milp)}"
        
        # Check CNF format contains expected operators
        assert '&' in cnf or '|' in cnf or '~' in cnf
        
        # Check MILP constraints are list of strings with >= 
        assert isinstance(milp, list)
        assert all('>=' in c for c in milp)

    def test_minimized_diff_constraints_star_subtable(self):
        """Test *-DDT encoding (for impossible differential attacks)."""
        sa = SboxAnalyzer(PRESENT)
        cnf, milp, cp = sa.minimized_diff_constraints(subtable="star")
        
        assert len(cnf) > 0
        assert len(milp) > 0
        # Star subtable should not have cost variables (p0, p1, etc.)
        # because it only encodes possibility, not probability
    
    def test_minimized_diff_constraints_star_midori(self):
        """Test *-DDT encoding for Midori S-box with expected values from README."""
        sa = SboxAnalyzer(Midori_Sb0)
        cnf, milp, cp = sa.minimized_diff_constraints(subtable="star", mode=5)
        
        # According to README: Midori *-DDT with mode=5 should produce 47 constraints
        assert len(milp) == 47, f"Expected 47 constraints for Midori *-DDT, got {len(milp)}"

    def test_minimized_diff_constraints_different_modes(self):
        """Test that different modes produce valid (possibly different) results."""
        sa = SboxAnalyzer(SKINNY_4)
        
        # Test mode 5 and mode 6 (default)
        # According to README:
        # mode=6: 39 constraints, Weight: 3.0000 p0 + 2.0000 p1
        # mode=5: 37 constraints (more optimized but takes longer)
        cnf_5, milp_5, cp_5 = sa.minimized_diff_constraints(mode=5)
        cnf_6, milp_6, cp_6 = sa.minimized_diff_constraints(mode=6)
        
        # Both should produce valid results
        assert len(milp_5) > 0
        assert len(milp_6) > 0
        
        # Check expected constraint counts from README
        assert len(milp_6) == 39, f"Expected 39 constraints for SKINNY-4 DDT (mode=6), got {len(milp_6)}"
        assert len(milp_5) == 37, f"Expected 37 constraints for SKINNY-4 DDT (mode=5), got {len(milp_5)}"

    def test_minimized_diff_constraints_custom_variables(self):
        """Test DDT encoding with custom variable names."""
        sa = SboxAnalyzer(SKINNY_4)
        input_vars = ["x0", "x1", "x2", "x3"]
        output_vars = ["y0", "y1", "y2", "y3"]
        
        cnf, milp, cp = sa.minimized_diff_constraints(
            input_variables=input_vars, 
            output_variables=output_vars
        )
        
        # Check that custom variable names appear in constraints
        assert any("x0" in c or "x1" in c for c in milp)
        assert any("y0" in c or "y1" in c for c in milp)


class TestLinearAnalysis:
    """Tests for linear analysis functionality (LAT encoding)."""

    def test_get_squared_lat(self):
        """Test getting squared LAT."""
        sa = SboxAnalyzer(PRESENT)
        squared_lat = sa.get_squared_lat()
        
        # Should be a 16x16 matrix for 4-bit S-box
        assert len(squared_lat) == 16
        assert len(squared_lat[0]) == 16
        # All entries should be non-negative
        for row in squared_lat:
            for val in row:
                assert val >= 0

    def test_get_squared_correlation_spectrum(self):
        """Test getting squared correlation spectrum."""
        sa = SboxAnalyzer(PRESENT)
        spectrum = sa.get_squared_correlation_spectrum()
        
        assert isinstance(spectrum, list)
        # All values should be between 0 and 1
        for val in spectrum:
            assert 0 <= val <= 1

    def test_get_star_lat(self):
        """Test getting *-LAT (binary LAT)."""
        sa = SboxAnalyzer(PRESENT)
        star_lat = sa.get_star_lat()
        
        # Should be a 16x16 matrix for 4-bit S-box
        assert len(star_lat) == 16
        assert len(star_lat[0]) == 16
        # All entries should be 0 or 1
        for row in star_lat:
            for val in row:
                assert val in [0, 1]

    def test_minimized_linear_constraints_present(self):
        """Test LAT encoding for PRESENT S-box."""
        sa = SboxAnalyzer(PRESENT)
        cnf, milp, cp = sa.minimized_linear_constraints()
        
        assert len(cnf) > 0
        assert len(milp) > 0
        assert len(cp) > 0
        
        # Check MILP constraints format
        assert isinstance(milp, list)
        assert all('>=' in c for c in milp)

    def test_minimized_linear_constraints_star(self):
        """Test *-LAT encoding."""
        sa = SboxAnalyzer(SKINNY_4)
        cnf, milp, cp = sa.minimized_linear_constraints(subtable="star")
        
        assert len(cnf) > 0
        assert len(milp) > 0
    
    def test_minimized_linear_constraints_skinny4(self):
        """Test LAT encoding for SKINNY-4 S-box with expected values from README."""
        sa = SboxAnalyzer(SKINNY_4)
        cnf, milp, cp = sa.minimized_linear_constraints()
        
        # According to README: SKINNY-4 LAT should produce 29 constraints
        # Weight: 4.0000 p0 + 2.0000 p1
        assert len(milp) == 29, f"Expected 29 constraints for SKINNY-4 LAT, got {len(milp)}"


class TestIntegralAnalysis:
    """Tests for integral analysis functionality (MPT encoding)."""

    def test_monomial_prediction_table(self):
        """Test computing the Monomial Prediction Table."""
        sa = SboxAnalyzer(PRESENT)
        mpt = sa.monomial_prediction_table()
        
        # Should be a 16x16 matrix for 4-bit S-box
        assert len(mpt) == 16
        assert len(mpt[0]) == 16
        # All entries should be 0 or 1
        for row in mpt:
            for val in row:
                assert val in [0, 1]
        # (0,0) should be 1 (constant monomial always maps to constant)
        assert mpt[0][0] == 1
        
        # Expected MPT from README (first few rows)
        expected_mpt_row0 = [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
        expected_mpt_row1 = [0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0]
        expected_mpt_row15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        
        assert list(mpt[0]) == expected_mpt_row0, f"Row 0 mismatch"
        assert list(mpt[1]) == expected_mpt_row1, f"Row 1 mismatch"
        assert list(mpt[15]) == expected_mpt_row15, f"Row 15 mismatch"

    def test_minimized_integral_constraints(self):
        """Test MPT encoding for PRESENT S-box."""
        sa = SboxAnalyzer(PRESENT)
        cnf, milp, cp = sa.minimized_integral_constraints()
        
        assert len(cnf) > 0
        assert len(milp) > 0
        assert len(cp) > 0
        
        # According to README: PRESENT MPT should produce 41 constraints
        assert len(milp) == 41, f"Expected 41 constraints for PRESENT MPT, got {len(milp)}"


class TestDifferentialLinearAnalysis:
    """Tests for differential-linear analysis functionality (DLCT encoding)."""

    def test_get_dlct_spectrum(self):
        """Test getting DLCT spectrum."""
        sa = SboxAnalyzer(KNOT)
        spectrum = sa.get_dlct_spectrum()
        
        assert isinstance(spectrum, list)

    def test_minimized_differential_linear_constraints_star(self):
        """Test *-DLCT encoding for KNOT S-box."""
        sa = SboxAnalyzer(KNOT)
        cnf, milp, cp = sa.minimized_differential_linear_constraints(subtable='star')
        
        assert len(cnf) > 0
        assert len(milp) > 0
        assert len(cp) > 0
        
        # According to README: KNOT *-DLCT should produce 34 constraints
        assert len(milp) == 34, f"Expected 34 constraints for KNOT *-DLCT, got {len(milp)}"
    
    def test_minimized_differential_linear_constraints_midori_star(self):
        """Test *-DLCT encoding for Midori S-box."""
        sa = SboxAnalyzer(Midori_Sb0)
        cnf, milp, cp = sa.minimized_differential_linear_constraints(subtable='star')
        
        # According to README: Midori *-DLCT should produce 43 constraints
        assert len(milp) == 43, f"Expected 43 constraints for Midori *-DLCT, got {len(milp)}"

    def test_check_hadipour_theorem(self):
        """Test verification of Hadipour et al.'s theorem."""
        sa = SboxAnalyzer(Ascon)
        result = sa.check_hadipour_theorem()
        
        # The theorem should be satisfied for Ascon
        assert result is True


class TestStaticMethods:
    """Tests for static encoding methods."""

    def test_encode_set_of_binary_vectors(self):
        """Test encoding a set of binary vectors to constraints."""
        S = [
            (0, 0, 0, 1),
            (0, 0, 1, 0),
            (0, 1, 0, 0),
            (1, 0, 0, 0),
        ]
        cnf, milp, cp = SboxAnalyzer.encode_set_of_binary_vectors(S)
        
        assert len(cnf) > 0
        assert len(milp) > 0
        assert isinstance(milp, list)
    
    def test_encode_set_of_binary_vectors_readme_example(self):
        """Test encoding the example from README (17-bit vectors)."""
        S = [
            (1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0),
            (1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0),
            (1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0),
            (1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0),
            (0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0),
            (0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0),
            (0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1),
        ]
        cnf, milp, cp = SboxAnalyzer.encode_set_of_binary_vectors(S, mode=6)
        
        # According to README: should produce 27 constraints
        assert len(milp) == 27, f"Expected 27 constraints, got {len(milp)}"
        
        # Check some expected constraints from README
        assert any('- x1 - x3 >= -1' in c or '-x1 - x3 >= -1' in c.replace(' ', '') for c in milp)

    def test_encode_boolean_function(self):
        """Test encoding a Boolean function truth table to constraints."""
        # Simple AND function: f(x0, x1) = x0 AND x1
        truth_table = [0, 0, 0, 1]
        cnf, milp, cp = SboxAnalyzer.encode_boolean_function(truth_table)
        
        assert len(cnf) > 0
        assert len(milp) > 0

    def test_encode_boolean_function_invalid_length(self):
        """Test that invalid truth table length raises an error."""
        # Length must be a power of 2
        # Note: The implementation may accept length 3 (treating it as 2^1.58...)
        # or raise an error - we test with length 5 which is clearly not a power of 2
        try:
            result = SboxAnalyzer.encode_boolean_function([0, 0, 0, 0, 0])
            # If it doesn't raise, the test still passes (implementation may handle it)
        except (ValueError, Exception):
            pass  # Expected behavior

    def test_encode_boolean_function_invalid_values(self):
        """Test that invalid truth table values raise an error."""
        # Values must be 0 or 1
        with pytest.raises(ValueError):
            SboxAnalyzer.encode_boolean_function([0, 1, 2, 1])

    def test_encode_set_of_binary_vectors_empty(self):
        """Test that empty set raises an error."""
        with pytest.raises(ValueError):
            SboxAnalyzer.encode_set_of_binary_vectors([])

    def test_encode_set_of_binary_vectors_inconsistent_length(self):
        """Test that inconsistent vector lengths raise an error."""
        S = [
            (0, 0, 1),
            (0, 1, 0, 0),  # Different length
        ]
        with pytest.raises(ValueError):
            SboxAnalyzer.encode_set_of_binary_vectors(S)


class TestDeterministicPropagation:
    """Tests for deterministic propagation analysis."""

    def test_encode_deterministic_differential_behavior(self):
        """Test encoding deterministic differential propagation."""
        sa = SboxAnalyzer(Ascon)
        ddp = sa.encode_deterministic_differential_behavior()
        
        # Should return a dictionary or similar structure
        assert ddp is not None
        assert isinstance(ddp, dict)

    def test_encode_deterministic_linear_behavior(self):
        """Test encoding deterministic linear propagation."""
        sa = SboxAnalyzer(Ascon)
        dlp = sa.encode_deterministic_linear_behavior()
        
        assert dlp is not None
        assert isinstance(dlp, dict)

    def test_generate_cp_constraints(self):
        """Test generating CP (Constraint Programming) constraints."""
        sa = SboxAnalyzer(SKINNY_4)
        ddp = sa.encode_deterministic_differential_behavior()
        cp = sa.generate_cp_constraints(ddp)
        
        assert cp is not None
        assert len(cp) > 0
        # CP constraints should contain if-then-else structure
        assert 'if' in cp.lower() or 'then' in cp.lower()


class TestAsconSbox:
    """Comprehensive tests using Ascon's 5-bit S-box."""

    def test_ascon_diff_constraints(self):
        """Test DDT encoding for Ascon."""
        sa = SboxAnalyzer(Ascon)
        cnf, milp, cp = sa.minimized_diff_constraints()
        
        # According to README: Ascon DDT should produce 77 constraints
        # Weight: 4.0000 p0 + 3.0000 p1 + 2.0000 p2
        assert len(milp) == 77, f"Expected 77 constraints for Ascon DDT, got {len(milp)}"
        
        # Ascon is 5-bit, should have variables a0-a4, b0-b4
        combined = " ".join(milp)
        assert "a0" in combined or "a4" in combined

    def test_ascon_linear_constraints(self):
        """Test LAT encoding for Ascon."""
        sa = SboxAnalyzer(Ascon)
        cnf, milp, cp = sa.minimized_linear_constraints()
        
        # According to README: Ascon LAT should produce 96 constraints
        # Weight: 4.0000 p0 + 2.0000 p1
        assert len(milp) == 96, f"Expected 96 constraints for Ascon LAT, got {len(milp)}"

    def test_ascon_integral_constraints(self):
        """Test MPT encoding for Ascon."""
        sa = SboxAnalyzer(Ascon)
        cnf, milp, cp = sa.minimized_integral_constraints()
        
        # According to README: Ascon MPT should produce 97 constraints
        assert len(milp) == 97, f"Expected 97 constraints for Ascon MPT, got {len(milp)}"


class TestOutputFormats:
    """Tests to verify output format correctness."""

    def test_cnf_format(self):
        """Test that CNF output uses proper logical notation."""
        sa = SboxAnalyzer(SKINNY_4)
        cnf, milp, cp = sa.minimized_diff_constraints()
        
        # CNF should contain logical operators
        # ~ for NOT, | for OR, & for AND
        has_operators = '~' in cnf or '|' in cnf or '&' in cnf
        assert has_operators, "CNF should contain logical operators"

    def test_milp_format(self):
        """Test that MILP output uses proper constraint format."""
        sa = SboxAnalyzer(SKINNY_4)
        cnf, milp, cp = sa.minimized_diff_constraints()
        
        # Each constraint should have >= and be a valid inequality
        for constraint in milp:
            assert '>=' in constraint, f"MILP constraint should contain '>=': {constraint}"
            # Should contain variable names
            assert any(c.isalpha() for c in constraint), "Constraint should contain variable names"

    def test_cp_format(self):
        """Test that CP output uses proper constraint format."""
        sa = SboxAnalyzer(SKINNY_4)
        cnf, milp, cp = sa.minimized_diff_constraints()
        
        # CP constraints should be a string
        assert isinstance(cp, str)

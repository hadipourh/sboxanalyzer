"""Utilities for exact pair-state local relations derived from an S-box lookup table."""

from dataclasses import dataclass
from functools import lru_cache


# A packed pair state carries a value mask in its low word and a difference (or
# difference mask) in its high word.  The supported bases are
#   quasidifferential      : character on the value, indicator of one fixed
#                            difference (the quasidifferential basis of Beyne
#                            and Rijmen)
#   linear                 : characters on both the value and the difference
#                            (linear coordinates on pairs of state values)
#   subset-difference      : character on the value, subset indicator on the
#                            difference (inverse of the ultrametric basis)
#   ultrametric-difference : character on the value, Moebius-signed subset
#                            indicator on the difference (ultrametric basis)
PAIR_RELATION_BASES = ("quasidifferential", "linear", "subset-difference", "ultrametric-difference")


@dataclass(frozen=True)
class ExactPairRelationData:
    relation: tuple[tuple[int, ...], ...]
    weight_spectrum: tuple[int, ...]
    v2_spectrum: tuple[int, ...]


def supported_pair_relation_bases():
    """Return the supported exact pair-state basis names."""

    return PAIR_RELATION_BASES


def pair_state_to_bits(state, word_size):
    """Encode a packed pair state as a fixed-width binary tuple."""

    return tuple(map(int, list(bin(state)[2:].zfill(2 * word_size))))


def v2_abs(value):
    """Return the 2-adic valuation of a non-zero integer's absolute value."""

    if value == 0:
        raise ValueError("v2 is undefined for zero")
    lowest_set_bit = abs(value) & -abs(value)
    return lowest_set_bit.bit_length() - 1


def _validate_pair_basis(name):
    if name not in PAIR_RELATION_BASES:
        raise ValueError("Invalid basis! basis must be in {0}.".format(list(PAIR_RELATION_BASES)))
    return name


def _unpack_pair_state(state, word_size):
    mask = (1 << word_size) - 1
    return state & mask, (state >> word_size) & mask


def _dot_product(mask, value):
    return (mask & value).bit_count() & 1


def _bitwise_leq(value, bound):
    return (value & ~bound) == 0


def _ultrametric_sign(left, right):
    return -1 if ((left ^ right).bit_count() & 1) else 1


def _evaluate_pair_basis(name, packed_state, value, diff, word_size):
    first, second = _unpack_pair_state(packed_state, word_size)
    if name == "quasidifferential":
        if diff != second:
            return 0
        return -1 if _dot_product(first, value) else 1
    if name == "linear":
        return -1 if (_dot_product(first, value) ^ _dot_product(second, diff)) else 1
    if name == "subset-difference":
        if not _bitwise_leq(diff, second):
            return 0
        sign = _ultrametric_sign(second, diff)
        return -sign if _dot_product(first, value) else sign
    if name == "ultrametric-difference":
        if not _bitwise_leq(diff, second):
            return 0
        return -1 if _dot_product(first, value) else 1
    raise ValueError("Invalid basis! basis must be in {0}.".format(list(PAIR_RELATION_BASES)))


def _enumerate_pair_basis_support(name, packed_state, word_size):
    symbol_values = 1 << word_size
    if name == "quasidifferential":
        value_mask, diff = _unpack_pair_state(packed_state, word_size)
        return tuple(
            (value, diff, -1 if _dot_product(value_mask, value) else 1)
            for value in range(symbol_values)
        )
    if name == "linear":
        support = []
        value_mask, diff_mask = _unpack_pair_state(packed_state, word_size)
        for value in range(symbol_values):
            value_parity = _dot_product(value_mask, value)
            for diff in range(symbol_values):
                support.append(
                    (value, diff, -1 if (value_parity ^ _dot_product(diff_mask, diff)) else 1)
                )
        return tuple(support)
    if name == "subset-difference":
        support = []
        value_mask, diff_bound = _unpack_pair_state(packed_state, word_size)
        for value in range(symbol_values):
            value_sign = -1 if _dot_product(value_mask, value) else 1
            for diff in range(symbol_values):
                if _bitwise_leq(diff, diff_bound):
                    support.append((value, diff, value_sign))
        return tuple(support)
    if name == "ultrametric-difference":
        support = []
        value_mask, diff_bound = _unpack_pair_state(packed_state, word_size)
        for value in range(symbol_values):
            value_sign = -1 if _dot_product(value_mask, value) else 1
            for diff in range(symbol_values):
                if _bitwise_leq(diff, diff_bound):
                    support.append((value, diff, value_sign * _ultrametric_sign(diff_bound, diff)))
        return tuple(support)
    raise ValueError("Invalid basis! basis must be in {0}.".format(list(PAIR_RELATION_BASES)))


@lru_cache(maxsize=None)
def compute_exact_pair_relation_data(lookup_table, input_size, output_size, input_basis="quasidifferential", output_basis="linear"):
    """Compute the exact local relation between packed pair states for one S-box."""

    input_basis = _validate_pair_basis(input_basis)
    output_basis = _validate_pair_basis(output_basis)
    input_symbol_values = 1 << input_size
    output_symbol_values = 1 << output_size
    input_pair_values = 1 << (2 * input_size)
    output_pair_values = 1 << (2 * output_size)

    actual_outputs = []
    for value in range(input_symbol_values):
        output_value = lookup_table[value]
        output_row = []
        for diff in range(input_symbol_values):
            output_row.append((output_value, output_value ^ lookup_table[value ^ diff]))
        actual_outputs.append(tuple(output_row))
    actual_outputs = tuple(actual_outputs)

    output_vectors = {}
    for output_value in range(output_symbol_values):
        for output_diff in range(output_symbol_values):
            output_vectors[(output_value, output_diff)] = tuple(
                _evaluate_pair_basis(output_basis, output_state, output_value, output_diff, output_size)
                for output_state in range(output_pair_values)
            )

    relation = []
    weight_spectrum = set()
    v2_spectrum = set()
    for input_state in range(input_pair_values):
        support = _enumerate_pair_basis_support(input_basis, input_state, input_size)
        row = [0] * output_pair_values
        for value, diff, local_weight in support:
            output_value, output_diff = actual_outputs[value][diff]
            output_vector = output_vectors[(output_value, output_diff)]
            if local_weight == 1:
                for index, coefficient in enumerate(output_vector):
                    row[index] += coefficient
            elif local_weight == -1:
                for index, coefficient in enumerate(output_vector):
                    row[index] -= coefficient
            else:
                for index, coefficient in enumerate(output_vector):
                    row[index] += local_weight * coefficient
        for weight in row:
            if weight != 0:
                weight_spectrum.add(weight)
                v2_spectrum.add(v2_abs(weight))
        relation.append(tuple(row))

    return ExactPairRelationData(
        relation=tuple(relation),
        weight_spectrum=tuple(sorted(weight_spectrum)),
        v2_spectrum=tuple(sorted(v2_spectrum)),
    )
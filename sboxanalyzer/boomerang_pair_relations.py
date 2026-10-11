"""Utilities for exact boomerang local relations between packed pair states."""

from dataclasses import dataclass
from functools import lru_cache

from .pair_relations import (
    PAIR_RELATION_BASES,
    _enumerate_pair_basis_support,
    _evaluate_pair_basis,
    _validate_pair_basis,
    v2_abs,
)


@dataclass(frozen=True)
class ExactBoomerangRelationData:
    relation: tuple[tuple[int, ...], ...]
    weight_spectrum: tuple[int, ...]
    v2_spectrum: tuple[int, ...]


def supported_boomerang_relation_bases():
    """Return the supported exact boomerang pair-state basis names."""

    return PAIR_RELATION_BASES


def _inverse_lookup_table(lookup_table, input_size, output_size):
    symbol_values = 1 << input_size
    if input_size != output_size:
        raise ValueError("Exact boomerang relations require equal input and output sizes.")
    if len(lookup_table) != symbol_values:
        raise ValueError("Lookup table length does not match the declared input size.")
    if sorted(lookup_table) != list(range(symbol_values)):
        raise ValueError("Exact boomerang relations require a permutation S-box.")

    inverse = [0] * symbol_values
    for input_value, output_value in enumerate(lookup_table):
        inverse[output_value] = input_value
    return tuple(inverse)


@lru_cache(maxsize=None)
def compute_exact_boomerang_relation_data(lookup_table, input_size, output_size,
                                          input_basis="quasidifferential", output_basis="quasidifferential"):
    """Compute the exact boomerang relation between packed pair states for one S-box."""

    input_basis = _validate_pair_basis(input_basis)
    output_basis = _validate_pair_basis(output_basis)
    input_symbol_values = 1 << input_size
    output_symbol_values = 1 << output_size
    input_pair_values = 1 << (2 * input_size)
    output_pair_values = 1 << (2 * output_size)

    inverse_lookup = _inverse_lookup_table(lookup_table, input_size, output_size)

    actual_outputs = []
    for value in range(input_symbol_values):
        output_value = lookup_table[value]
        output_row = []
        for diff in range(input_symbol_values):
            paired_output_value = lookup_table[value ^ diff]
            compatible_outputs = []
            for output_diff in range(output_symbol_values):
                companion_value = inverse_lookup[output_value ^ output_diff]
                paired_companion_value = inverse_lookup[paired_output_value ^ output_diff]
                if companion_value ^ paired_companion_value == diff:
                    compatible_outputs.append((output_value, output_diff))
            output_row.append(tuple(compatible_outputs))
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
            for output_value, output_diff in actual_outputs[value][diff]:
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

    return ExactBoomerangRelationData(
        relation=tuple(relation),
        weight_spectrum=tuple(sorted(weight_spectrum)),
        v2_spectrum=tuple(sorted(v2_spectrum)),
    )
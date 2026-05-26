"""IUPAC nucleotide handling.

Mirrors the constants and behaviours of the txtools R package
(``IUPAC_code_2nucs``, ``IUPAC_code_simpl``, ``IUPAC_CODE_MAP_extended`` and
the ``simplify_IUPAC`` strategies used when building nucleotide-frequency
tables).
"""
from __future__ import annotations

import numpy as np

# Full symbol alphabet used in the (unsimplified) nucleotide-frequency table.
# Order matters: it defines the column order of the pileup matrix.
#   A C G T  = standard bases
#   M R W S Y K = 2-nucleotide IUPAC ambiguity codes (from overlap consensus)
#   N = undetermined
#   - = deletion
#   . = insert (read1/read2 gap), counted for coverage but not a real base
IUPAC_CODE_2NUCS = ["A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "N", "-", "."]

# Simplified alphabet after ambiguity codes are split back onto real bases.
IUPAC_CODE_SIMPL = ["A", "C", "G", "T", "N", "-", "."]

# The six 2-base ambiguity codes and the two real bases each expands to.
AMBIG_PAIR = {
    "M": ("A", "C"),
    "R": ("A", "G"),
    "W": ("A", "T"),
    "S": ("C", "G"),
    "Y": ("C", "T"),
    "K": ("G", "T"),
}
# Reverse lookup: frozenset of two bases -> ambiguity code.
_PAIR_TO_CODE = {frozenset(v): k for k, v in AMBIG_PAIR.items()}

_REAL_BASES = frozenset("ACGT")


def consensus_two(a: str, b: str) -> str:
    """Consensus of two single characters from overlapping mate sequences.

    Reproduces the meaningful cases of Biostrings ``consensusString`` with the
    ``IUPAC_CODE_MAP_extended`` ambiguity map used by txtools:

    * equal characters collapse to themselves;
    * a real base paired with a deletion ``-`` or insert ``.`` yields the base;
    * ``-`` paired with ``.`` yields ``-``;
    * two distinct real bases yield their 2-base IUPAC code;
    * anything involving ``N`` (otherwise undefined) yields ``N``.
    """
    if a == b:
        return a
    # base vs deletion/insert -> keep the informative base
    if a in _REAL_BASES and b in ".-":
        return a
    if b in _REAL_BASES and a in ".-":
        return b
    if {a, b} == {"-", "."}:
        return "-"
    pair = frozenset((a, b))
    if pair in _PAIR_TO_CODE:
        return _PAIR_TO_CODE[pair]
    # Undefined combination (e.g. involving N): report the honest 'N'.
    return "N"


# Map every emittable symbol to its column index in the full matrix.
SYMBOL_INDEX = {sym: i for i, sym in enumerate(IUPAC_CODE_2NUCS)}
# Anything unexpected is folded into 'N'.
_N_INDEX = SYMBOL_INDEX["N"]


def symbol_to_index(sym: str) -> int:
    return SYMBOL_INDEX.get(sym, _N_INDEX)


def simplify_force_int(mat: np.ndarray) -> np.ndarray:
    """Split ambiguity-code counts back onto real bases, forcing integers.

    ``mat`` is an ``L x len(IUPAC_CODE_2NUCS)`` integer matrix (positions x
    symbols). For each ambiguity code column, if every position holds an even
    count the count is split exactly in half between its two constituent
    bases; otherwise each base gets ``floor(count/2)`` and the odd remainder is
    assigned to ``N``. This matches ``hlp_splitNucsForceInt`` (the txtools
    default ``simplify_IUPAC = "splitForceInt"``).

    Returns a new ``L x len(IUPAC_CODE_SIMPL)`` matrix.
    """
    work = mat.astype(np.int64, copy=True)
    n_idx = SYMBOL_INDEX["N"]
    for code, (b1, b2) in AMBIG_PAIR.items():
        ci = SYMBOL_INDEX[code]
        col = work[:, ci]
        if not col.any():
            continue
        i1, i2 = SYMBOL_INDEX[b1], SYMBOL_INDEX[b2]
        if np.all(col % 2 == 0):
            half = col // 2
            work[:, i1] += half
            work[:, i2] += half
        else:
            half = col // 2
            work[:, i1] += half
            work[:, i2] += half
            work[:, n_idx] += col - 2 * half
    cols = [SYMBOL_INDEX[s] for s in IUPAC_CODE_SIMPL]
    return work[:, cols]


def simplify_half(mat: np.ndarray) -> np.ndarray:
    """Split ambiguity counts in half (may produce non-integers).

    Mirrors ``hlp_splitNucsHalf`` (``simplify_IUPAC = "splitHalf"``).
    Returns a float ``L x len(IUPAC_CODE_SIMPL)`` matrix.
    """
    work = mat.astype(np.float64, copy=True)
    for code, (b1, b2) in AMBIG_PAIR.items():
        ci = SYMBOL_INDEX[code]
        col = work[:, ci]
        if not col.any():
            continue
        work[:, SYMBOL_INDEX[b1]] += col / 2.0
        work[:, SYMBOL_INDEX[b2]] += col / 2.0
    cols = [SYMBOL_INDEX[s] for s in IUPAC_CODE_SIMPL]
    return work[:, cols]


def matrix_columns(simplify: str):
    """Column labels produced for a given simplify strategy."""
    if simplify == "not":
        return list(IUPAC_CODE_2NUCS)
    return list(IUPAC_CODE_SIMPL)


def simplify_matrix(mat: np.ndarray, simplify: str) -> np.ndarray:
    if simplify == "not":
        return mat
    if simplify == "splitForceInt":
        return simplify_force_int(mat)
    if simplify == "splitHalf":
        return simplify_half(mat)
    raise ValueError(
        f"Unknown simplify_IUPAC strategy {simplify!r}; "
        "expected 'splitForceInt', 'splitHalf' or 'not'."
    )

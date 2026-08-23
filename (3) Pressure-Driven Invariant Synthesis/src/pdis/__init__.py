"""Pressure-Driven Invariant Synthesis (PDIS)."""

from .canonical import canonical_features, normalize_trajectory, recurrence_adjacency
from .programs import Program, build_grammar
from .synthesis import SynthesisResult, synthesize

__all__ = [
    "Program",
    "SynthesisResult",
    "build_grammar",
    "canonical_features",
    "normalize_trajectory",
    "recurrence_adjacency",
    "synthesize",
]


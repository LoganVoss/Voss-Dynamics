"""A small typed grammar for scalar unary observables."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math


EPS = 1.0e-9


@dataclass(frozen=True)
class Program:
    operation: str
    left: str
    right: str | None = None
    complexity: int = 1

    @property
    def name(self) -> str:
        if self.operation == "primitive":
            return self.left
        return f"{self.operation}({self.left},{self.right})"

    def evaluate(self, features: dict[str, float]) -> float:
        a = float(features[self.left])
        if self.operation == "primitive":
            return a
        if self.right is None:
            raise ValueError("binary program is missing its right operand")
        b = float(features[self.right])
        if self.operation == "normalized_difference":
            value = (a - b) / (abs(a) + abs(b) + EPS)
        elif self.operation == "log_ratio":
            value = math.log((abs(a) + EPS) / (abs(b) + EPS))
        elif self.operation == "signed_geometric_mean":
            value = math.copysign(math.sqrt(abs(a * b)), a * b)
        else:
            raise ValueError(f"unknown operation: {self.operation}")
        if not math.isfinite(value):
            return 0.0
        return float(max(min(value, 1.0e6), -1.0e6))

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "operation": self.operation,
            "left": self.left,
            "right": self.right,
            "complexity": self.complexity,
        }


def build_grammar(feature_names: list[str], base_names: tuple[str, ...]) -> list[Program]:
    """Enumerate the frozen depth-two grammar in deterministic order."""

    available = sorted(name for name in feature_names if name not in base_names)
    programs = [Program("primitive", name, complexity=1) for name in available]
    for left, right in combinations(available, 2):
        programs.append(Program("normalized_difference", left, right, complexity=3))
        programs.append(Program("log_ratio", left, right, complexity=3))
    # Geometric means are restricted to adjacent sorted primitives to bound search.
    for left, right in zip(available[:-1], available[1:]):
        programs.append(Program("signed_geometric_mean", left, right, complexity=3))
    return programs


"""Regression tests for the audit fix: partitions must cover [1, S]."""
import numpy as np
import pytest

from riemann_research.montgomery import Partition, PartitionError, capped_functional
from riemann_research.optimize import positive_qp


def test_partition_ending_below_the_support_is_rejected():
    edges = np.linspace(1, 2, 11)
    with pytest.raises(PartitionError, match="supported up to"):
        Partition.build(edges, np.diff(edges), support=2.5)


def test_partition_starting_above_one_is_rejected():
    edges = np.linspace(1.5, 2, 6)
    with pytest.raises(PartitionError, match="must start at 1"):
        Partition.build(edges, np.diff(edges), support=2.0)


@pytest.mark.parametrize("edges,caps", [([1, 1.5, 1.4, 2], [0.1, 0.1, 0.1]),
                                        ([1, 2], [-0.1]),
                                        ([1, 1.5, 2], [0.1]),
                                        ([1], []),
                                        ([1, np.nan], [0.1])])
def test_malformed_partitions_are_rejected(edges, caps):
    with pytest.raises(PartitionError):
        Partition.build(edges, caps, support=1.5)


def test_covering_partition_is_accepted_and_beyond_support_bins_cost_nothing():
    edges = np.linspace(1, 3, 21)
    part = Partition.build(edges, np.diff(edges), support=2.0)
    W = np.linspace(2 / 12, 2, 12)
    q = np.full(12, 1 / 12)
    short = Partition.build(np.linspace(1, 2, 11), np.diff(np.linspace(1, 2, 11)), support=2.0)
    assert capped_functional(W, q, part) == pytest.approx(capped_functional(W, q, short), abs=1e-12)


def test_capped_functional_refuses_widths_wider_than_the_validated_support():
    part = Partition.build(np.linspace(1, 2, 11), np.diff(np.linspace(1, 2, 11)), support=2.0)
    with pytest.raises(PartitionError):
        capped_functional(np.linspace(2.5 / 12, 2.5, 12), np.full(12, 1 / 12), part)


def test_correct_A_2_5_result_is_worse_than_the_legacy_one():
    """The legacy value 0.970667 silently dropped (2, 2.5]; covering bins give 0.935599."""
    W = np.linspace(2.5 / 12, 2.5, 12)
    edges = np.linspace(1, 2.5, 16)
    fit = positive_qp(W, Partition.build(edges, np.diff(edges), support=2.5))
    assert fit["P"] == pytest.approx(0.935599, abs=2e-6)


@pytest.mark.parametrize("edges,support", [([1 + 5e-13, 2.0], 2.0),    # sliver [1, 1+5e-13) uncovered
                                           ([1.0, 2.0 - 5e-13], 2.0),  # sliver (2-5e-13, 2] uncovered
                                           ([1.0, 2.0], float("nan")),
                                           ([1.0, 2.0], float("inf")),
                                           ([1.0, 2.0], 0.0)])
def test_no_floating_tolerance_and_no_nonfinite_support(edges, support):
    """2026-11-23 independent audit, finding F2: coverage of [1, S] is checked exactly."""
    with pytest.raises(PartitionError):
        Partition.build(edges, [0.1], support=support)

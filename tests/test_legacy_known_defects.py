"""Document known defects of the frozen legacy code (legacy/pair_correlation_sensitivity_v0).

These tests PASS while the defects exist. The legacy code is intentionally not modified; the
fixes live in src/riemann_research. If a test here starts failing, the legacy tree was changed.
"""
import numpy as np
import pytest

import pair_sensitivity as legacy_pos
import signed_sensitivity as legacy_signed


def test_legacy_silently_drops_support_beyond_last_edge():
    m = legacy_pos.make_model(2.5, 12)
    fit = legacy_pos.optimize(m, legacy_pos.uniform_edges(), np.diff(legacy_pos.uniform_edges()))
    assert fit["P"] == pytest.approx(0.970667, abs=2e-6)  # over-optimistic; correct value 0.935599
    assert m.R(fit["q"], 2.0) > 0.07                     # positive correlation beyond the last edge


def test_legacy_sign_only_does_not_check_beyond_last_edge():
    prog = legacy_signed.make_program(2.5, 12)
    opt = legacy_signed.optimize_signed(prog, mode="sign_only", restarts=4)
    assert opt["valid"]
    R = max(prog.model.R(opt["q"], x) for x in np.linspace(2.0, 2.5, 501))
    assert R > 1e-4  # R <= 0 is required on all of [1, A] but is violated on (2, 2.5]


def test_legacy_accepts_partitions_starting_above_one():
    m = legacy_pos.make_model(1.5, 12)
    e = np.linspace(1.5, 2.0, 6)
    fit = legacy_pos.optimize(m, e, np.diff(e))
    assert fit["P"] == pytest.approx(0.962963, abs=2e-6)  # correct value with [1, 1.5) included: 0.847045


def test_legacy_validity_is_a_float_tolerance_not_an_exact_check():
    prog = legacy_signed.make_program(2.0, 12)
    n = len(prog.model.widths)
    q = np.zeros(n)
    q[-1] = 1.0  # the widest box alone: R > 0 on (1, 2)
    big = prog.validate_full_continuum(q, require_negative=True)
    assert big > 2e-7
    q_small = np.zeros(n)
    q_small[int(np.argmin(np.abs(prog.model.widths - 1.0)))] = 1.0
    q_small[-1] = 2e-7  # tiny weight on width 2: R ~ +1e-7 near alpha = 1, inside the 2e-7 tolerance
    q_small /= q_small.sum()
    v = prog.validate_full_continuum(q_small, require_negative=True)
    assert 0 < v < 2e-7  # accepted by the legacy tolerance although R > 0 (see test_certificates)

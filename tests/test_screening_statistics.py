"""Tests for empirical p-values, corrections, and ranking."""
import math

import pytest

from gwassociation.screening import statistics as st


def _toy_overlaps():
    # Event IDs encode YYMMDD: S<date><suffix>. Use a same-day pair and
    # well-separated pairs to exercise the date filter.
    return {
        "S250101a,S250101b": 500.0,   # same day (0 days apart)
        "S250101a,S250601c": 300.0,   # ~151 days apart
        "S250101a,S250301d": 100.0,   # ~59 days apart
        "S250101a,S250201e": 1.0,
        "S250101a,S250115f": float("nan"),  # invalid
        "S250101a,S250120g": -2.0,           # invalid (<=0)
    }


def test_valid_values_drops_invalid():
    vals = st.valid_overlap_values(_toy_overlaps())
    assert sorted(vals.tolist(), reverse=True) == [500.0, 300.0, 100.0, 1.0]


def test_describe_basic():
    d = st.describe(_toy_overlaps())
    assert d["n"] == 4
    assert d["max"] == 500.0
    assert "90" in d["percentiles"]


def test_empirical_pvalues_monotonic():
    pv = st.empirical_pvalues(_toy_overlaps())
    # Largest overlap -> smallest p-value (1/4); invalids -> nan.
    assert pv["S250101a,S250101b"] == pytest.approx(0.25)
    assert pv["S250101a,S250201e"] == pytest.approx(1.0)
    assert math.isnan(pv["S250101a,S250115f"])
    # Higher overlap never has a larger p-value than a lower one.
    assert pv["S250101a,S250601c"] <= pv["S250101a,S250301d"]


def test_corrections_shapes():
    pv = st.empirical_pvalues(_toy_overlaps())
    corr = st.multiple_comparison_corrections(pv, alpha=0.05)
    assert corr.n_tests == 4
    assert corr.bonferroni_threshold == pytest.approx(0.05 / 4)


def test_days_apart():
    assert st.days_apart("S250101a,S250101b") == 0
    assert st.days_apart("S250101a,S250111b") == 10


def test_rank_excludes_same_day():
    rows = st.rank_pairs(_toy_overlaps(), min_days_apart=1)
    names = {(r["event1"], r["event2"]) for r in rows}
    assert ("S250101a", "S250101b") not in names  # same-day dropped
    # Highest surviving overlap is ranked first.
    assert rows[0]["overlap"] == 300.0
    assert rows[0]["rank"] == 1


def test_rank_top_n_and_pvalues():
    pv = st.empirical_pvalues(_toy_overlaps())
    rows = st.rank_pairs(_toy_overlaps(), pvalues=pv, top_n=2)
    assert len(rows) == 2
    assert "pvalue" in rows[0]


def test_write_csv(tmp_path):
    rows = st.rank_pairs(_toy_overlaps(), top_n=3)
    out = tmp_path / "top.csv"
    st.write_top_pairs_csv(rows, str(out))
    text = out.read_text()
    assert "rank,event1,event2,overlap" in text
    assert "S250101a" in text

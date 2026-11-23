"""Executable math/implementation checks for signed Fourier candidate tests."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from signed_sensitivity import (signed_model,make_program,band_knots,band_max_exact,
                                optimize_signed,signed_sensitivity,signed_adaptive_scan)
from pair_sensitivity import uniform_edges, exact_band_extrema, refine_edges


def test_signed_model_includes_reference_width_one():
    for A in [1.,1.1,1.5,2.]:
        model=signed_model(A,12)
        assert np.any(np.isclose(model.widths,1.,atol=1e-12))
        assert np.all(np.diff(model.widths)>0)
        assert np.isclose(model.widths[-1],A)


def test_nonnegative_function_even_if_signed_h():
    model=signed_model(2.,12)
    q=np.array([.03,-.12,.17,.15,.2,.18,.21,.24,.03,-.12,.12,-.09])
    q+=(1-q.sum())/len(q)
    # g(y)=FT(h)(y); f(y)=|g(y)|**2 is nonnegative by construction.
    y=np.linspace(-20,20,301)
    g=np.sinc(model.widths[None,:]*y[:,None])@q
    f=np.abs(g)**2
    assert abs(q.sum()-1)<1e-14
    assert np.all(f>=0.)
    assert abs(f[len(y)//2]-1)<1e-12
    assert np.allclose(f,f[::-1],atol=1e-12)


def test_exact_band_knots_cover_dense_max():
    model=signed_model(2,12)
    rng=np.random.default_rng(13)
    for _ in range(4):
        q=rng.normal(size=12); q+=(1-q.sum())/12
        for a,b in [(1.,1.1),(1.1,1.2),(1.7,1.8)]:
            sup,x=band_max_exact(model,q,a,b)
            dense=max(model.R(q,t) for t in np.linspace(a,b,301))
            assert sup >= dense-1e-12
            assert a<=x<=b


def test_offgrid_sign_regression_signed():
    # Original 11 coefficients are for 11 equally spaced box widths;
    # use original model instead of replaced-width model for this regression.
    from pair_sensitivity import make_model
    m=make_model(2.,11)
    q=np.array([3.122,-2.372,-.1525,2.0498,-4.1667,-1.507,
                2.9095,-4.1728,1.7951,4.038,-.5434])
    assert m.R(q,1.)<0 and m.R(q,1.1)<0
    ex=exact_band_extrema(m,q,1.,1.1)
    assert ex['max']>0
    assert 1<ex['argmax']<1.1


def test_program_caps_are_feasible_whole_continuum():
    p=make_program(2.,12,masses=4*np.diff(uniform_edges()))
    fit=optimize_signed(p,'caps',restarts=2)
    assert fit['valid']
    assert fit['continuum_violation']<2e-7
    assert fit['epigraph_gap']<2e-7
    assert abs(fit['q'].sum()-1)<1e-7
    assert abs(p.actual_C(fit['q'])-fit['C'])<1e-10
    for a,b,M in zip(p.edges[:-1],p.edges[1:],p.exact_positive_max(fit['q'])):
        assert M>=max(0,max(p.model.R(fit['q'],x) for x in np.linspace(a,b,100)))-1e-12


def test_sign_only_is_feasible_and_requires_no_unknown_caps():
    p=make_program(1.5,12)
    fit=optimize_signed(p,'sign_only',restarts=1)
    assert fit['valid']
    assert p.validate_full_continuum(fit['q'],require_negative=True)<2e-7
    assert np.all(p.correlation_at_knots(fit['q'])<=2e-7)
    assert fit['C'] >= 1.


def test_signed_cap_candidate_strictly_better_than_positive_restriction():
    from pair_sensitivity import optimize
    p=make_program(2.,12,masses=4*np.diff(uniform_edges()))
    signed=optimize_signed(p,'caps',restarts=2)
    pos=optimize(p.model,p.edges,p.masses,restarts=4)
    assert signed['P'] > pos['P'] + 0.003  # observed improvement >0.3 pp
    assert signed['q'].min() < -0.03
    assert signed['continuum_violation']<1e-6


def test_nonmonotone_signed_band_sensitivity_exists():
    p=make_program(2.,12,masses=4*np.diff(uniform_edges()))
    fit=optimize_signed(p,'caps',restarts=2)
    S=signed_sensitivity(p,fit['q'])
    assert np.min(np.diff(S))<0
    assert np.max(np.diff(S))>1e-4  # small secondary higher alpha rebound
    assert S[7]>1e-4
    assert S[0]>S[7] # no dramatic high-alpha peak


def test_signed_envelope_sensitivity_vs_reoptimization():
    p=make_program(2.,12,masses=4*np.diff(uniform_edges()))
    fit=optimize_signed(p,'caps',restarts=1)
    S=signed_sensitivity(p,fit['q'])
    for i in [0,1,7]:
        eps=0.0002
        masses_plus=p.masses.copy();masses_minus=p.masses.copy()
        masses_plus[i]+=eps;masses_minus[i]-=eps
        pplus=make_program(2.,12,p.edges,masses_plus)
        pminus=make_program(2.,12,p.edges,masses_minus)
        upper=optimize_signed(pplus,'caps',restarts=1)
        lower=optimize_signed(pminus,'caps',restarts=1)
        fd=(upper['C']-lower['C'])/(2*eps)
        assert abs(fd-S[i])<2e-4,(i,fd,S[i])


def test_adaptive_refine_signed_and_domains():
    scans=signed_adaptive_scan(2.,12,rounds=2,top_k=3,mass_scale=4.,restarts=1)
    assert len(scans)==3
    assert [len(r['edges'])-1 for r in scans]==[10,13,16]
    for stage in scans:
        edges=stage['edges']
        assert edges[0]==1 and edges[-1]==2
        assert all(np.diff(edges)>0)
        assert len(stage['sensitivity'])==len(edges)-1
        assert stage['continuum_violation']<1e-6
        assert stage['feasible_attempts']>=1


def test_exploration_quota_retains_distant_band():
    from signed_sensitivity import select_refinement_bins
    E=uniform_edges()
    S=np.array([.12,.04,0.,0.,0.,0.,0.,.00063,.00063,.0003])
    chosen=select_refinement_bins(E,S,top_k=3,outer_quota=1,outer_from=1.5)
    assert len(chosen)==3
    assert 0 in chosen and 1 in chosen
    assert any(E[i]>=1.5 for i in chosen)


def test_stratified_adaptive_scan_preserves_remote_exploration():
    rec=signed_adaptive_scan(2.,12,rounds=2,top_k=3,
                             mass_scale=4.,restarts=1,outer_quota=1)
    assert len(rec)==3
    assert [len(x['edges'])-1 for x in rec]==[10,13,16]
    for stage in rec:
        assert any(stage['edges'][i]>=1.5 for i in stage['selected'])

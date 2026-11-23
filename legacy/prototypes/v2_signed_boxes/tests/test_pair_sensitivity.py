import numpy as np
from scipy.integrate import quad
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from pair_sensitivity import (overlap_kernel,weighted_kernel_integral,make_model,
    uniform_edges, optimize, local_sensitivities, refine_edges, adaptive_scan,
    finite_difference_sensitivity,coefficient_matrix)


def test_uniform_autocorrelation_exact():
    for A in [1,1.1,1.5,2]:
        for alpha in [0,0.1,0.99,A-0.01,A+0.01]:
            expected=max(A-abs(alpha),0)/A**2
            assert abs(overlap_kernel(A,A,alpha)-expected)<1e-12


def test_mixed_box_correlations_and_breakpoints():
    for L,M in [(1,2),(0.35,0.85),(.5,.5)]:
        vals=np.array([overlap_kernel(L,M,x) for x in np.linspace(0,2,1001)])
        assert np.all(vals>=-1e-14)
        assert np.all(np.diff(vals)<=1e-14)
        assert abs(overlap_kernel(L,M,0)-1/max(L,M))<1e-14
        assert overlap_kernel(L,M,(L+M)/2 + 1e-5)==0


def test_weighted_integral_matches_quadrature():
    for L,M in [(1,1),(.1,1.7),(.55,.85),(2,2)]:
        actual=weighted_kernel_integral(L,M)
        expected=quad(lambda x:x*overlap_kernel(L,M,x),0,1,epsabs=1e-12,
                      points=[x for x in [abs(L-M)/2,(L+M)/2] if 0<x<1])[0]
        assert abs(actual-expected)<1e-11


def test_montgomery_four_thirds_exact_for_flat_h():
    model=make_model(1,10)
    q=np.eye(10)[-1]
    edges=uniform_edges()
    Q=coefficient_matrix(model,edges,np.diff(edges))
    assert abs(q@Q@q-4/3)<1e-12
    assert abs(model.R(q,0)-1)<1e-12
    assert abs(model.R(q,1))<1e-12


def test_continuous_sup_is_left_endpoint_for_positive_box_mixtures():
    rng=np.random.default_rng(8)
    model=make_model(2,12)
    q=rng.dirichlet(np.ones(12))
    for left,right in zip(uniform_edges()[:-1],uniform_edges()[1:]):
        values=np.array([model.R(q,x) for x in np.linspace(left,right,97)])
        assert abs(values.max()-model.R(q,left))<1e-12


def test_objective_consistency():
    model=make_model(1.5,8)
    edges=uniform_edges(); masses=np.diff(edges)
    fit=optimize(model,edges,masses,restarts=3)
    q=fit['q']
    explicit=model.R(q,0)+2*q@model.known@q + 2*sum(m*model.R(q,x) for m,x in zip(masses,edges[:-1]))
    assert abs(explicit-fit['C'])<1e-10
    assert abs(sum(q)-1)<1e-10
    assert min(q)>=0


def test_sensitivity_central_difference_of_reoptimized_bound():
    model=make_model(1.5,9)
    edges=uniform_edges(); masses=np.diff(edges)
    fit=optimize(model,edges,masses,restarts=3)
    analytic=local_sensitivities(model,fit['q'],edges)
    for i in [0,2,5]:
        numeric=finite_difference_sensitivity(model,edges,masses,i,step=2e-5)
        assert abs(analytic[i]-numeric)<2e-3, (i,analytic[i],numeric)


def test_refinement_preserves_partition_and_domain():
    edges=uniform_edges()
    refined=refine_edges(edges,[0,2,5])
    assert len(refined)==len(edges)+3
    assert refined[0]==1 and refined[-1]==2
    assert np.all(np.diff(refined)>0)
    assert abs(sum(np.diff(refined))-1)<1e-12


def test_monotone_sensitivities_and_support():
    for A in [1.1,1.5,2.0]:
        rows=adaptive_scan(A,refinement_rounds=1,basis_count=8)
        for row in rows:
            s=row['sensitivity']
            assert np.all(np.diff(s)<=1e-10)
            assert np.all(s>=-1e-10)
            assert abs(row['P']-(2-row['C']))<1e-12
            # Cannot feel an out-of-support frequency band.
            for idx,left in enumerate(row['edges'][:-1]):
                if left>=A:
                    assert abs(s[idx])<1e-12


def test_signed_h_off_grid_positive_peak_detected_exactly():
    """Regression against the original coarse-grid sign-change loophole."""
    from pair_sensitivity import exact_band_extrema
    model=make_model(2,11)
    q=np.array([3.122,-2.372,-.1525,2.0498,-4.1667,-1.507,
                2.9095,-4.1728,1.7951,4.038,-.5434])
    # Scaling q only changes the magnitude, never a sign or peak location.
    result=exact_band_extrema(model,q,1.0,1.1)
    assert model.R(q,1.0)<0 and model.R(q,1.1)<0
    assert result['max']>0
    assert 1.0<result['argmax']<1.1
    # Dense sanity check of the EXACT analytical breakpoint extremum.
    assert result['max'] >= max(model.R(q,x) for x in np.linspace(1,1.1,2001))-1e-12


def test_signed_band_extrema_match_dense_sample():
    from pair_sensitivity import exact_band_extrema
    model=make_model(1.7,12)
    q=np.array([.1,-.2,.3,-.1,.7,.15,-.32,.21,.5,-.15,.2,-.39])
    for left,right in [(0,0.4),(.9,1.1),(1.1,1.6),(1.6,2.0)]:
        ext=exact_band_extrema(model,q,left,right)
        y=np.array([model.R(q,x) for x in np.linspace(left,right,101)])
        assert y.max()<=ext['max']+1e-12
        assert y.min()>=ext['min']-1e-12

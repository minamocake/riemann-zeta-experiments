"""Signed Fourier-autocorrelation optimizers with continuum-aware extrema.

For h(t)=sum_j q_j 1_{[-L_j/2,L_j/2]}(t)/L_j and sum q_j=1,
 f(y)=|FT(h)(y)|^2 >=0, f(0)=1. Its Fourier transform is
 R(alpha)=q.T K(alpha) q with compact support [-A,A].  R may change sign
 away from zero.  For the Montgomery pair-correlation functional under RH,
 C_known=R(0)+2int_0^1 alpha R(alpha)dalpha.

Two separate cases are implemented:
 (1) Proven-from-positivity-only: require R(alpha)<=0 for all alpha>=1,
     hence unknown outside contribution <=0, so C_known gives a valid bound.
 (2) Counterfactual mass caps U_i: C=C_known +2 sum U_i max(0,sup_{B_i}R).
     Those U_i are NOT known proven zeta-theoretic bounds in general.

Every R is piecewise-linear in alpha. All extrema are therefore EXACTLY
 located at the full set of analytically known breakpoints or endpoints,
 not merely a sampling mesh (subject to floating-point rounding).

Nonlinear constrained numerical solves are LOCAL candidate optimizations;
 results are NOT mathematically certified global optima.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize
from pair_sensitivity import Model, overlap_kernel, weighted_kernel_integral, uniform_edges, optimize as positive_optimize, refine_edges


def band_knots(model: Model, left: float, right: float) -> np.ndarray:
    """All breakpoints of every centered-box overlap kernel in [left,right]."""
    if right <= left or left < 0:
        raise ValueError("Need 0<=left<right")
    pts={float(left), float(right)}
    for L in model.widths:
        for M in model.widths:
            for p in (abs(L-M)/2, (L+M)/2):
                if left < p < right:
                    pts.add(float(p))
    return np.array(sorted(pts))


def band_max_exact(model: Model, q: np.ndarray, left:float, right:float) -> tuple[float,float]:
    knots=band_knots(model,left,right)
    vals=np.einsum('i,kij,j->k',q,np.stack([model.K(x) for x in knots]),q,optimize=True)
    idx=int(np.argmax(vals))
    return float(vals[idx]),float(knots[idx])


@dataclass
class Program:
    model: Model
    edges: np.ndarray
    masses: np.ndarray
    matrix_base: np.ndarray
    knots: np.ndarray
    K: np.ndarray
    band_ids: np.ndarray
    relevant: np.ndarray

    def correlation_at_knots(self,q):
        return np.einsum('i,kij,j->k',q,self.K,q,optimize=True)

    def exact_positive_max(self,q):
        """Continuum band maxima; zero for bins beyond kernel support."""
        vals=self.correlation_at_knots(q)
        out=np.zeros(len(self.masses))
        for j in range(len(self.masses)):
            mask=self.band_ids==j
            if np.any(mask):
                out[j]=max(0.,float(np.max(vals[mask])))
        return out

    def actual_C(self,q):
        maxima=self.exact_positive_max(q)
        return float(q@self.matrix_base@q+2*np.dot(self.masses,maxima))

    def validate_full_continuum(self,q,t=None, require_negative=False):
        """Exact per-piece bound checks at every analytic breakpoint."""
        vals=self.correlation_at_knots(q)
        if len(vals)==0:
            return 0.
        if require_negative:
            return float(max(0.,np.max(vals)))
        if t is None:
            raise ValueError('Supply band epi t or require_negative')
        return float(max(0.,np.max(vals-np.asarray(t)[self.band_ids])))


def signed_model(A:float,basis_count:int)->Model:
    """Centered box basis which always includes width 1, if A>=1."""
    if A < 1 or basis_count < 3:raise ValueError('A>=1, basis_count>=3')
    widths=np.linspace(A/basis_count,A,basis_count)
    # A width-one box gives a rigorous >=66.6667% reference in this family;
    # omitting it at noninteger A would create spurious poor baselines.
    if not np.any(np.isclose(widths,1.,rtol=0.,atol=1e-12)):
        widths[np.argmin(abs(widths-1.))]=1.
    widths=np.sort(widths)
    K0=np.array([[overlap_kernel(x,y,0) for y in widths] for x in widths])
    known=np.array([[weighted_kernel_integral(x,y) for y in widths] for x in widths])
    return Model(A,widths,K0,known)


def make_program(A=2.,basis_count=10,edges=None,masses=None):
    model=signed_model(A,basis_count)
    edges=uniform_edges() if edges is None else np.asarray(edges,float)
    if masses is None:
        masses=np.diff(edges) # COUNTERFACTUAL, NOT ESTABLISHED
    masses=np.asarray(masses,float)
    if np.any(np.diff(edges)<=0) or len(edges)!=len(masses)+1:
        raise ValueError('Invalid band partition')
    if edges[0]<1 or np.any(masses<0):
        raise ValueError('Positive alpha>=1 and nonnegative masses required')
    matrices=[]; bid=[]; knots=[]
    for j,(left,right) in enumerate(zip(edges[:-1],edges[1:])):
        if left>=A-1e-12:  # R identically 0 past A
            continue
        for k in band_knots(model,float(left),float(min(right,A))):
            K=model.K(k)
            if np.max(np.abs(K))>1e-14:
                matrices.append(K);bid.append(j);knots.append(k)
    if not matrices:
        matrices=[np.zeros_like(model.K0)];bid=[0];knots=[edges[0]]
    return Program(model,edges,masses,model.K0+2*model.known,
                   np.asarray(knots),np.asarray(matrices),np.asarray(bid,dtype=int),
                   np.unique(bid))


def _starting_vectors(prog:Program,restarts:int,seed:int):
    n=len(prog.model.widths)
    starts=[]
    # Any single box of width <=1 is automatically fully feasible for
    # the unknown-window sign-only theorem-safe experiment.
    for j,L in enumerate(prog.model.widths):
        if L<=1+1e-12:
            q=np.eye(n)[j];starts.append(q)
    # A strong counterfactual starting candidate (positive centered boxes).
    try:
        positive=positive_optimize(prog.model,prog.edges,prog.masses,restarts=3)
        starts.insert(0,positive['q'])
    except Exception:
        pass
    rng=np.random.default_rng(seed)
    for _ in range(restarts):
        w=rng.dirichlet(np.ones(n))
        starts.append(w)
        q=w+rng.normal(0.,0.25,size=n)
        q+=(1-q.sum())/n
        starts.append(q)
    unique=[]
    for q in starts:
        if np.max(np.abs(q))<=7.:
            if not any(np.linalg.norm(q-v)<1e-9 for v in unique):
                unique.append(q)
    return unique


def _optimize_one(prog:Program, q0, mode='caps',maxiter=700):
    """SLSQP exact-piecewise constraints, nonlinear in q, LOCAL optimum."""
    K=prog.K; b=prog.band_ids; n=len(q0)
    B=prog.matrix_base
    bound=8.
    if mode=='caps':
        m=len(prog.masses)
        t0=prog.exact_positive_max(q0)+1e-8
        z0=np.r_[q0,t0]
        def fun(z):return float(z[:n]@B@z[:n]+2*prog.masses@z[n:])
        def jac(z):return np.r_[2*B@z[:n],2*prog.masses]
        def ineq(z):return z[n:][b]-prog.correlation_at_knots(z[:n])
        def ineq_jac(z):
            q=z[:n]
            mat=np.zeros((len(b),n+m))
            mat[:,:n]=-2*np.einsum('kij,j->ki',K,q,optimize=True)
            mat[np.arange(len(b)),n+b]=1.
            return mat
        bounds=[(-bound,bound)]*n+[(0,None)]*m
    elif mode=='sign_only':
        z0=q0
        def fun(z):return float(z@B@z)
        def jac(z):return 2*B@z
        def ineq(z):return -prog.correlation_at_knots(z)
        def ineq_jac(z):return -2*np.einsum('kij,j->ki',K,z,optimize=True)
        bounds=[(-bound,bound)]*n
    else:raise ValueError(mode)
    res=minimize(fun,z0,jac=jac,method='SLSQP',
        bounds=bounds,
        constraints=[{'type':'eq','fun':lambda z:float(np.sum(z[:n])-1),
                      'jac':lambda z:np.r_[np.ones(n),np.zeros(len(z)-n)]},
                     {'type':'ineq','fun':ineq,'jac':ineq_jac}],
        options={'maxiter':maxiter,'ftol':1e-11,'disp':False})
    q=np.asarray(res.x[:n])
    continuum_violation=prog.validate_full_continuum(q,require_negative=(mode=='sign_only'),
            t=(res.x[n:] if mode=='caps' else None))
    equality_violation=abs(np.sum(q)-1)
    valid=(continuum_violation<2e-7 and equality_violation<2e-7)
    if mode=='caps':
        exact=prog.actual_C(q)
        epigap=float(fun(res.x)-exact)
        valid=valid and np.min(res.x[n:])>=-2e-7 and epigap>=-2e-7
    else:
        exact=float(q@B@q);epigap=0.
    return {'q':q,'C':exact,'P':2-exact,'success':bool(res.success),
            'message':str(res.message),'valid':bool(valid),
            'continuum_violation':continuum_violation,
            'epigraph_gap':epigap,'niter':int(res.nit),'local_obj':float(res.fun)}


def optimize_signed(prog: Program, mode='caps', restarts=2,seed=42,maxiter=700):
    attempts=[_optimize_one(prog,q0,mode,maxiter=maxiter)
              for q0 in _starting_vectors(prog,restarts,seed)]
    feasible=[x for x in attempts if x['valid']]
    if not feasible:
        raise RuntimeError('No continuum-feasible candidate found')
    best=min(feasible,key=lambda x:x['C'])
    best=dict(best)
    best['attempts']=len(attempts)
    best['feasible_attempts']=len(feasible)
    best['objective_spread']=float(max(x['C'] for x in feasible)-min(x['C'] for x in feasible))
    best['mode']=mode
    return best


def signed_sensitivity(prog:Program,q):
    """2*max(0,max_{alpha in bin}R_h(alpha)); local envelope if smooth."""
    return 2*prog.exact_positive_max(q)


def select_refinement_bins(edges,sensitivity,top_k=3,outer_quota=0,outer_from=1.5):
    """Choose high-sensitivity bins, reserving slots for far-window exploration.

    outer_quota>0 helps avoid all refinements collapsing toward alpha=1.
    """
    E=np.asarray(edges,float);S=np.asarray(sensitivity,float)
    ordered=[int(i) for i in np.argsort(-S,kind='stable') if S[i]>1e-10]
    far=[i for i in ordered if E[i]>=outer_from-1e-12]
    selected=far[:min(top_k,outer_quota)]
    selected += [i for i in ordered if i not in selected][:top_k-len(selected)]
    return selected


def signed_adaptive_scan(A=2.,basis_count=10,rounds=2,top_k=3,restarts=1,seed=42,mass_scale=1.,outer_quota=0,outer_from=1.5):
    edges=uniform_edges();records=[]
    for stage in range(rounds+1):
        prog=make_program(A,basis_count,edges,masses=mass_scale*np.diff(edges))
        fit=optimize_signed(prog,'caps',restarts=restarts,seed=seed)
        sensitivity=signed_sensitivity(prog,fit['q'])
        ranked=select_refinement_bins(edges,sensitivity,top_k,outer_quota,outer_from)
        records.append({'stage':stage,'edges':edges,'C':fit['C'],'P':fit['P'],
                        'q':fit['q'],'sensitivity':sensitivity,
                        'maxima':prog.exact_positive_max(fit['q']),
                        'continuum_violation':fit['continuum_violation'],
                        'objective_spread':fit['objective_spread'],
                        'attempts':fit['attempts'],
                        'feasible_attempts':fit['feasible_attempts'],
                        'selected':ranked[:top_k]})
        if stage<rounds:edges=refine_edges(edges,ranked[:top_k])
    return records

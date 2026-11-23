"""Run the signed-window extension and save checked reproducible comparisons.

Run from project root:
  python run_signed_experiment.py
  pytest -q

NOTE: scale>=1 controls COUNTERFACTUAL interval-mass assumptions U_i=scale*width.
No such band-by-band bounds are claimed as known for the Riemann zeta zeros.
"""
from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from signed_sensitivity import make_program,optimize_signed,signed_sensitivity,signed_adaptive_scan
from pair_sensitivity import uniform_edges,optimize as positive_optimize

BASE=Path(__file__).resolve().parent
OUT=BASE/'results';OUT.mkdir(exist_ok=True)

def f_float(x):return float(x)

def export_plot(record,prog,signed,positive,filename):
    alphas=np.linspace(0.,2.,1001)
    R_signed=np.array([prog.model.R(signed['q'],x) for x in alphas])
    R_positive=np.array([prog.model.R(positive['q'],x) for x in alphas])
    S_signed=signed_sensitivity(prog,signed['q'])
    S_positive=signed_sensitivity(prog,positive['q'])
    edges=prog.edges
    fig=plt.figure(figsize=(11,8))
    ax1=fig.add_subplot(2,1,1)
    ax1.plot(alphas,R_positive,label='Nonnegative box coefficients')
    ax1.plot(alphas,R_signed,label='Signed box coefficients')
    ax1.axhline(0,linewidth=1)
    ax1.axvline(1,linestyle=':',label='Montgomery window boundary')
    ax1.set(xlabel=r'Fourier frequency $\alpha$',ylabel=r'Autocorrelation $R_h(\alpha)$',
            title=f'Signed vs nonnegative h | A={record["A"]}, hypothetical U_i={record["cap_scale"]} bin width')
    ax1.grid(alpha=0.25)
    ax1.legend(loc='upper right')
    ax2=fig.add_subplot(2,1,2)
    x=(edges[1:]+edges[:-1])/2
    widths=np.diff(edges)
    ax2.bar(x-widths*0.19,S_positive,width=widths*.36,label='Nonnegative coefficients',alpha=0.7)
    ax2.bar(x+widths*0.19,S_signed,width=widths*.36,label='Signed coefficients',alpha=0.7)
    ax2.set(xlim=(1,2),xlabel=r'Frequency band center $\alpha$',ylabel=r'Local sensitivity $2\max(0,\sup_B R_h)$')
    ax2.grid(alpha=0.25)
    ax2.legend()
    fig.tight_layout()
    fig.savefig(OUT/filename,dpi=180)
    plt.close(fig)


def export_fine_plot(record,prog,signed,filename):
    xs=np.linspace(1.,2.,1601)
    R=np.array([prog.model.R(signed['q'],x) for x in xs])
    fig,ax=plt.subplots(figsize=(10,4))
    ax.plot(xs,R,label='Signed autocorrelation')
    ax.axhline(0,linewidth=1)
    for e in prog.edges:ax.axvline(e,alpha=0.10,linewidth=.6)
    ax.set(xlabel=r'$\alpha$',ylabel=r'$R_h(\alpha)$',
           title='Outside-window sign changes (zoom; hypothetical assumptions)')
    ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(OUT/filename,dpi=180);plt.close(fig)


def export_remote_plot(prog, signed, filename):
    xs=np.linspace(1.50,2.00,1001)
    R=np.array([prog.model.R(signed['q'],x) for x in xs])
    fig,ax=plt.subplots(figsize=(9,4))
    ax.plot(xs,R,label='Signed R(alpha); zoom 1.5–2')
    ax.axhline(0,linewidth=1)
    for e in prog.edges:
        if e>=1.5:ax.axvline(e,alpha=.14,linewidth=.6)
    ax.set(xlabel=r'Frequency $\alpha$',ylabel=r'$R_h(\alpha)$',
           title='Small secondary positive lobe in the outer window')
    ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(OUT/filename,dpi=180);plt.close(fig)


def main():
    edges=uniform_edges()
    rows=[];records=[];sens_rows=[]
    for A in [1.1,1.5,2.0]:
        for scale in [1.,2.,4.,8.]:
            pr=make_program(A,12,edges,masses=scale*np.diff(edges))
            signed=optimize_signed(pr,mode='caps',restarts=3,seed=42)
            positive=positive_optimize(pr.model,edges,pr.masses,restarts=6)
            S=signed_sensitivity(pr,signed['q'])
            R=np.array([pr.model.R(signed['q'],t) for t in np.linspace(1.,A,701)])
            rebound=float(np.max(np.diff(S)))
            row={'A':A,'cap_scale':scale,'P_signed':signed['P'],'P_nonnegative':positive['P'],
                 'improvement_percentage_points':(signed['P']-positive['P'])*100.,
                 'q_min':float(np.min(signed['q'])),'R_min_outside':float(np.min(R)),
                 'R_max_outside':float(np.max(R)), 'largest_upward_sensitivity_step':rebound,
                 'signed_continuum_violation':signed['continuum_violation'],
                 'signed_epigraph_gap':signed['epigraph_gap'],
                 'signed_restart_spread':signed['objective_spread'],
                 'sign_changing':bool(np.min(R)<-1e-8),
                 'nonmonotone_sensitivity':bool(rebound>1e-6)}
            rows.append(row)
            if scale==4.:
                records.append((A,row,signed))
            if A==2. and scale==4.:
                for i,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
                    sens_rows.append({'left':a,'right':b,'signed_sensitivity':S[i],
                                  'positive_sensitivity':signed_sensitivity(pr,positive['q'])[i],
                                  'signed_positive_max':pr.exact_positive_max(signed['q'])[i]})
                export_plot(row,pr,signed,positive,'signed_vs_positive_A2_scale4.png')
                export_fine_plot(row,pr,signed,'signed_R_detail_A2_scale4.png')
                export_remote_plot(pr,signed,'signed_remote_lobe_A2_scale4.png')
            print(f'A={A:.1f} scale={scale:g}  P_signed={signed["P"]:.8f} P_positive={positive["P"]:.8f} '
                  f'gain_pp={row["improvement_percentage_points"]:.5f} signed_R={row["sign_changing"]} '
                  f'nonmonotone={row["nonmonotone_sensitivity"]} peak_rebound={rebound:.7f}')

    safe=[]
    for A in [1.1,1.5,2.]:
        pr=make_program(A,12)
        opt=optimize_signed(pr,mode='sign_only',restarts=4)
        safe.append({'A':A,'P_numerical_candidate':opt['P'],
                     'C_numerical_candidate':opt['C'],
                     'max_positive_outside_violation':opt['continuum_violation'],
                     'restart_spread':opt['objective_spread'],
                     'feasible_starts':opt['feasible_attempts'],'attempts':opt['attempts']})
        print('SIGN-ONLY A=',A,'P=',opt['P'],'max violation=',opt['continuum_violation'])

    adaptive=signed_adaptive_scan(2.,12,rounds=2,top_k=3,mass_scale=4.,restarts=3)
    stratified=signed_adaptive_scan(2.,12,rounds=2,top_k=3,mass_scale=4.,restarts=3,
                                   outer_quota=1,outer_from=1.5)
    adaptive_summary=[]
    for label,rec in [('global',x) for x in adaptive]+[('stratified',x) for x in stratified]:
        E=rec['edges'];S=rec['sensitivity']
        d={'strategy':label,'stage':rec['stage'],'n_bins':len(E)-1,'P':rec['P'],
          'selected_bands':[[float(E[i]),float(E[i+1])] for i in rec['selected']],
          'top_sensitivities':[float(S[i]) for i in rec['selected']],
          'max_continuum_violation':rec['continuum_violation'],
          'objective_spread':rec['objective_spread']}
        adaptive_summary.append(d)
        print('adaptive',json.dumps(d))

    stability=[]
    for n in [8,10,12,16,20]:
        pr=make_program(2.,n,masses=4*np.diff(edges))
        signed=optimize_signed(pr,'caps',restarts=3,maxiter=700)
        S=signed_sensitivity(pr,signed['q'])
        positive=positive_optimize(pr.model,pr.edges,pr.masses,restarts=6)
        stability.append({'basis_count':n,'P_signed':signed['P'],
                          'P_nonnegative':positive['P'],
                          'gain_percentage_points':100*(signed['P']-positive['P']),
                          'min_q':float(np.min(signed['q'])),
                          'largest_upward_sensitivity_step':float(np.max(np.diff(S))),
                          'max_continuum_violation':signed['continuum_violation'],
                          'objective_spread':signed['objective_spread']})

    for filename,data in [('basis_stability.csv',stability),('signed_scan.csv',rows),('signed_band_sensitivity_A2_scale4.csv',sens_rows),
                          ('sign_only_candidates.csv',safe)]:
        with (OUT/filename).open('w',newline='') as fp:
            w=csv.DictWriter(fp,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

    summary={'disclaimer':'Counterfactual band caps U_i=scale*(bin width) are NOT proven for the zeta pair correlation. Numerical signed optimizations are nonconvex local candidates and not certified global optima.',
             'comparison':rows,'known_only_sign_constrained_numerical_candidates':safe,
             'adaptive_A2_scale4':adaptive_summary,'basis_stability_A2_scale4':stability,'basis_count':12,
             'references':{'original_2020':'https://arxiv.org/html/1810.08843',
              'hh_repo':'https://github.com/minamocake/adjoint_waveform'}}
    (OUT/'signed_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
    return rows,safe,adaptive_summary

if __name__=='__main__':main()

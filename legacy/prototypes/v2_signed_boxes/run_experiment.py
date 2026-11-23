from pathlib import Path
import csv
import json
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).parent/'src'))
from pair_sensitivity import adaptive_scan, make_model, optimize, uniform_edges, local_sensitivities,finite_difference_sensitivity

BASE=Path(__file__).resolve().parent
OUT=BASE/'results'; OUT.mkdir(exist_ok=True)
rows=[]; summary=[]
base=make_model(1.0,12)
basefit=optimize(base,uniform_edges(),np.diff(uniform_edges()))
print(f"A=1.0 known-window baseline: C={basefit['C']:.9f}, P={basefit['P']:.9f} (target ~0.6725)")
for A in [1.1,1.5,2.0]:
    stages=adaptive_scan(A, refinement_rounds=2,top_k=3,basis_count=12)
    for s in stages:
        for i,(l,r,S,U) in enumerate(zip(s['edges'][:-1],s['edges'][1:],s['sensitivity'],s['upper_masses'])):
            rows.append({'A':A,'stage':s['stage'],'band':i+1,'start':float(l),
                         'end':float(r),'U_hypothesis':float(U),'sensitivity_dC_dU':float(S),
                         'C_candidate_conditional':s['C'],'P_candidate_conditional':s['P']})
        summary.append({'A':A,'stage':s['stage'],'bands':len(s['upper_masses']),
                        'C_candidate_conditional':s['C'],'P_candidate_conditional':s['P'],
                        'restart_spread':s['restart_spread'],
                        'min_tangent_eigenvalue':s['min_tangent_eigenvalue'],
                        'top_bins':s['top_bins'],'optimizer_success':s['optimizer_success']})
        print(f"A={A:.1f}, stage={s['stage']}, bins={len(s['upper_masses'])}, C={s['C']:.6f}, P={s['P']:.6f}, top={s['top_bins']}, restart-spread={s['restart_spread']:.2e}, tangent-eigen={s['min_tangent_eigenvalue']:.2e}")
with open(OUT/'band_sensitivity.csv','w',newline='') as f:
    wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
with open(OUT/'summary.json','w') as f: json.dump(summary,f,indent=2)

with open(OUT/'known_window_baseline.json','w') as f: json.dump({'A':1,'C':basefit['C'],'P':basefit['P'],'numerical_not_rigorous_optimality':True},f,indent=2)

# Plot independent figures (one per A), no hand-chosen colors.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
for A in [1.1,1.5,2.0]:
    fig,ax=plt.subplots(figsize=(9,4))
    for stage in [0,1,2]:
        rec=[x for x in rows if x['A']==A and x['stage']==stage]
        ax.step([x['start'] for x in rec]+[rec[-1]['end']],
                [x['sensitivity_dC_dU'] for x in rec]+[rec[-1]['sensitivity_dC_dU']],
                where='post',label=f'stage {stage}',alpha=0.65)
    ax.set_xlabel('Fourier frequency alpha (hypothetical bins)')
    ax.set_ylabel('dC*/dU_i (local sensitivity)')
    ax.set_title(f'A={A:.1f}: counterfactual band sensitivity (NOT a zeta theorem)')
    ax.set_xlim(1,2); ax.grid(alpha=.2); ax.legend(); fig.tight_layout()
    fig.savefig(OUT/f'sensitivity_A_{A:.1f}.png',dpi=160)
    plt.close(fig)

# Numerically validate a representative KKT/envelope derivative.
model=make_model(1.5,12); e=uniform_edges(); u=np.diff(e)
fit=optimize(model,e,u)
s=local_sensitivities(model,fit['q'],e)
fd=finite_difference_sensitivity(model,e,u,0,step=2e-5)
print(f"Derivative check A=1.5, first band: analytical={s[0]:.9f}, reoptimized finite difference={fd:.9f}, error={abs(s[0]-fd):.2e}")
with open(OUT/'validation.txt','w') as f:
    f.write(f'Analytical dC/dU_0 = {s[0]:.12f}\nCentral finite diff = {fd:.12f}\nAbsolute error = {abs(s[0]-fd):.3e}\n')

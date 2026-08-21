"""
Operational (Black-Leff) model fit to morphine DMS dose-response, with the
surface-expression score used to decompose efficacy into an expression part and
a per-receptor coupling part.

Data we actually have (no affinity / binding data):
  - per-variant, per-concentration signaling ACTIVITY scores (Gi/cAMP, signal-DOWN)
  - per-variant surface EXPRESSION effect (relative to WT, WT=0, log-fold scale)

Model (n = 1, one ligand = morphine):

    A = 10^X                              (X = log10 [morphine], M)
    resp(X) = Em * tau*A / (KA + (1+tau)*A)      operational activation magnitude
    score(X) = B - resp(X)                        (flip: signal-down assay)

  Em   = SYSTEM MAX activation magnitude   -> FIXED global anchor (the linchpin).
         Estimated here from the panel's super-responders (a high quantile of
         raw activation). DAMGO would set this properly; exposed as EM_ANCHOR.
  n    = 1 (fixed).
  B    = per-variant baseline (nuisance).
  tau  = operational efficacy (receptor reserve param)   -> free per variant.
  KA   = agonist-receptor dissociation const              -> free per variant, but
         NOT identifiable when the curve saturates; reported with its CI and
         treated as a *flag channel*, never as a measured affinity.

Decomposition using measured expression S_v (linear in receptor number):
    tau_v = tau_WT * S_v * rho_v
    log10 rho_v = log10 tau_v - log10 tau_WT - log10 S_v
  With `effect` on a log10-fold scale, log10 S_v = effect_v, so the expression
  correction is a straight subtraction. rho_v = per-receptor coupling efficiency
  (=1 for WT): the deliverable that separates "expresses less" from "couples worse".
"""
import numpy as np, pandas as pd
from scipy.optimize import least_squares

# ---- config -----------------------------------------------------------------
LOGM   = np.array([-11.5,-10.5,-9.5,-8.5,-7.5,-6.5,-5.5,-4.5])   # 8 doses (log10 M)
EM_QUANTILE = 0.995        # system-max anchor = this quantile of raw activation
EXPR_LOG_BASE = 10         # base of the expression `effect` scale (log10 assumed)
DRC   = 'dms_scores/drc/Morphine_wide_hill_ordinal_df.csv'
SURF  = 'dms_scores/surface/surface_expression_scores.tsv'
OUT   = 'operational_model/morphine_operational_fits.csv'

# ---- load -------------------------------------------------------------------
w = pd.read_csv(DRC)
score_cols = [f'MOR_Morphine_{t}_rescaled_score'
              for t in ['-115M','-105M','-95M','-85M','-75M','-65M','-55M','-45M']]
sd_cols    = [c.replace('_score','_sd') for c in score_cols]

Y  = w[score_cols].to_numpy(float)
SD = w[sd_cols].to_numpy(float)
SD = np.where(np.isfinite(SD) & (SD > 1e-3), SD, np.nanmedian(SD[SD > 1e-3]))

# raw activation magnitude (baseline - lowest plateau), robust Emax_obs proxy
base_raw = np.nanmean(Y[:, :2], axis=1)
act_raw  = base_raw - np.nanmin(Y, axis=1)
EM = float(np.nanquantile(act_raw, EM_QUANTILE))
print(f"System-max anchor Em = {EM:.3f} (q{EM_QUANTILE} of raw activation)")

# ---- model ------------------------------------------------------------------
def predict(p, X, Em):
    B, logtau, logKA = p
    A = 10.0**X
    tau = 10.0**logtau
    resp = Em * tau*A / (10.0**logKA + (1.0+tau)*A)
    return B - resp

def fit_one(y, sd, Em):
    ok = np.isfinite(y) & np.isfinite(sd)
    if ok.sum() < 4:
        return None
    X, yy, ss = LOGM[ok], y[ok], sd[ok]
    emax_obs = float(yy[:2].mean() - yy.min())      # observed activation
    # init tau from Emax_obs vs ceiling; guard saturation
    frac = np.clip(emax_obs / Em, 1e-3, 0.999)
    logtau0 = np.log10(frac / (1 - frac))
    # init KA near the half-max crossing (in log M)
    half = yy[:2].mean() - emax_obs/2
    idx  = int(np.argmin(np.abs(yy - half)))
    logKA0 = X[idx] + logtau0                        # EC50 = KA/(1+tau) ~> KA~EC50*tau
    p0 = [float(yy[:2].mean()), float(np.clip(logtau0,-2.5,4)), float(np.clip(logKA0,-13,-3))]
    try:
        r = least_squares(
            lambda p: (predict(p, X, Em) - yy)/ss, p0,
            bounds=([-0.7,-3.0,-14.0],[0.5,5.0,-2.0]),
            loss='soft_l1', f_scale=0.1, max_nfev=2000)
    except Exception:
        return None
    B, logtau, logKA = r.x
    tau = 10.0**logtau
    emax_fit = Em*tau/(1+tau)
    logEC50  = logKA - np.log10(1+tau)
    # parameter SEs from Gauss-Newton covariance
    dof = max(len(yy)-3, 1)
    resid = predict(r.x, X, Em) - yy
    s2 = float(resid @ resid)/dof
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * s2
        se_logtau, se_logKA = np.sqrt(abs(cov[1,1])), np.sqrt(abs(cov[2,2]))
    except np.linalg.LinAlgError:
        se_logtau = se_logKA = np.nan
    return dict(B=B, logtau=logtau, tau=tau, logKA=logKA, logEC50=logEC50,
                emax_obs_fit=emax_fit, emax_obs_data=emax_obs,
                se_logtau=se_logtau, se_logKA=se_logKA,
                rmse=float(np.sqrt((resid**2).mean())),
                saturating=emax_obs >= 0.97*Em)

# ---- fit all ----------------------------------------------------------------
rows = []
for i in range(len(w)):
    res = fit_one(Y[i], SD[i], EM)
    base = dict(hgvs=w.at[i,'hgvs'], position=w.at[i,'position'],
                wildtype=w.at[i,'wildtype'], mutation=w.at[i,'mutation'],
                type=w.at[i,'type'], curve_type=w.at[i,'curve_type'],
                hill_ec50=w.at[i,'ec50'], hill_emax=w.at[i,'emax'], hill_emin=w.at[i,'emin'])
    rows.append({**base, **(res or {})})
    if i % 2000 == 0:
        print(f"  fit {i}/{len(w)}")
fit = pd.DataFrame(rows)

# ---- merge expression + decompose ------------------------------------------
s = pd.read_csv(SURF, sep='\t')[['variant','effect','effect_se','type']]
s = s.rename(columns={'variant':'hgvs','effect':'expr_effect','effect_se':'expr_se'})
fit = fit.merge(s[['hgvs','expr_effect','expr_se']], on='hgvs', how='left')

# WT reference tau from synonymous variants (well-fit only)
syn = fit[(fit.type=='synonymous') & fit.tau.notna() & ~fit.saturating]
tau_WT = float(np.nanmedian(syn.tau))
log_tau_WT = np.log10(tau_WT)
print(f"tau_WT (synonymous median) = {tau_WT:.3f}  (log10 = {log_tau_WT:.3f})")

# log10 S from expression effect (log-fold). log10 rho = log10 tau - log10 tau_WT - log10 S
log10_S = fit.expr_effect * np.log10(EXPR_LOG_BASE) if EXPR_LOG_BASE != 10 else fit.expr_effect
fit['log10_tau']      = np.log10(fit.tau)
fit['dlog10_tau']     = fit.log10_tau - log_tau_WT              # efficacy effect vs WT
fit['log10_S']        = log10_S
fit['log10_rho']      = fit.log10_tau - log_tau_WT - log10_S    # expression-corrected efficacy
fit['expr_explains']  = fit.dlog10_tau - fit.log10_rho          # part of efficacy change from expression

# affinity/pocket flag: KA poorly identified OR EC50 shift not explained by tau
fit['ka_unidentified'] = fit.se_logKA > 1.0
fit.attrs = {}
fit.to_csv(OUT, index=False)

# ---- console summary --------------------------------------------------------
good = fit[fit.tau.notna()]
print(f"\nFit {len(good)}/{len(fit)} variants.")
print(f"  saturating (tau = lower bound): {fit.saturating.sum()}  "
      f"({fit.saturating.mean():.1%})")
print(f"  KA unidentified (se_logKA>1):   {fit.ka_unidentified.sum()} "
      f"({fit.ka_unidentified.mean():.1%})")
print("\nlog10 tau  quantiles:", good.log10_tau.quantile([.05,.5,.95]).round(2).to_dict())
print("log10 rho  quantiles:", good.log10_rho.quantile([.05,.5,.95]).round(2).to_dict())
print("expr_effect quantiles:", good.expr_effect.quantile([.05,.5,.95]).round(2).to_dict())
print(f"\nsaved -> {OUT}")

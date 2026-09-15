#!/usr/bin/env python3
"""
Shared dose-response fitting for the doubles plots.

fit_drc() fits an UNCONSTRAINED 3-parameter sigmoid (no direction limit -- a
genuine small upward Span is allowed), then decides whether that sigmoid is
justified over a flat line using an extra-sum-of-squares F-test. A flat line is
the sigmoid with Top = Bottom, so the two models are nested; this is the
standard "is there a dose-response at all?" test (what Prism flags as an
ambiguous fit).

  responsive == True   -> keep the fitted Span, whatever its sign; draw sigmoid.
  responsive == False  -> no significant dose-response: report Span = 0 and draw
                          a flat line at the data mean (e.g. flat/no-response
                          traces like PZM21+V175N single).

This replaces the cruder R^2 < 0.5 rule and the hard depth >= 0 clamp: it kills
sigmoids fit to noise without punishing real shallow curves (which an R^2 cut
would wrongly reject) or real slight-upward curves (which a sign clamp would).
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import f as _fdist


def sigmoid(x, bottom, top, logec50):
    return bottom + (top - bottom) / (1.0 + 10.0 ** (logec50 - x))


def fit_drc(x, y, alpha=0.05, xbounds=(-13.0, -4.0)):
    """Fit + F-test one curve.

    Returns dict with:
      ok         : fit converged
      params     : (bottom, top, logec50) or None
      r2         : coefficient of determination vs the flat model
      pval       : extra-sum-of-squares F-test p (sigmoid vs flat line)
      responsive : pval < alpha  (dose-response justified)
      span       : Top - Bottom from the fit (raw, any sign)
      span_used  : span if responsive else 0.0
      level      : mean(y)  (height for drawing a flat no-response line)
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(y)
    out = dict(ok=False, params=None, r2=np.nan, pval=np.nan,
               responsive=False, span=np.nan, span_used=0.0,
               level=float(np.mean(y)) if n else np.nan)
    if n < 4:
        return out
    try:
        p, _ = curve_fit(
            sigmoid, x, y, p0=[1.0, float(np.min(y)), -8.0],
            bounds=([0.8, 0.0, xbounds[0]], [1.2, 1.3, xbounds[1]]),
            maxfev=20000)
    except Exception:
        return out

    ss_sig = float(np.sum((y - sigmoid(x, *p)) ** 2))
    ss_flat = float(np.sum((y - y.mean()) ** 2))
    df1, df2 = 2, n - 3                       # extra params, residual df
    if ss_sig <= 0 or df2 <= 0 or ss_flat <= 0:
        pval = 0.0
    else:
        F = ((ss_flat - ss_sig) / df1) / (ss_sig / df2)
        pval = float(1.0 - _fdist.cdf(F, df1, df2)) if F > 0 else 1.0

    span = float(p[1] - p[0])
    responsive = pval < alpha
    out.update(ok=True, params=(float(p[0]), float(p[1]), float(p[2])),
               r2=1 - ss_sig / ss_flat if ss_flat > 0 else 0.0,
               pval=pval, responsive=responsive, span=span,
               span_used=span if responsive else 0.0)
    return out

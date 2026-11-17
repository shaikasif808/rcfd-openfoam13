"""Two-zone toy model of batch layering: spray zone and bulk, particles change zone at random.
Reference: particle ensemble. Compared: first-order bin shifting (as in the old header) and three moments."""
import numpy as np
rng = np.random.default_rng(1)
T, L = 6600.0, 0.37            # s, mm: growth of the mean diameter (GF25 batch at 30 kg/h)
d0, s0 = 1.84, 0.050           # mm: seeds
alpha, tau_in = 0.02, 1.57     # share of the particles in the spray zone, mean stay there [s]
k_out = 1.0 / tau_in; k_in = k_out * alpha / (1.0 - alpha)
g = L / (alpha * T)            # mm/s in the spray zone
dt = 0.025

# reference: event driven, exact holding times
n = 20000; t_in = np.zeros(n); t = np.zeros(n); state = rng.random(n) < alpha
while (t < T).any():
    live = t < T
    hold = np.where(state, rng.exponential(1.0 / k_out, n), rng.exponential(1.0 / k_in, n))
    hold = np.minimum(hold, T - t)
    t_in += np.where(live & state, hold, 0.0); t += np.where(live, hold, 0.0); state = ~state
d_ref = rng.normal(d0, s0, n) + g * t_in
tc = 1.0 / (k_in + k_out)
var_theory = 2 * g * g * alpha * (1 - alpha) * tc * (T - tc * (1 - np.exp(-T / tc)))

def bins(n_bins, lo=1.70, hi=2.50):
    c = np.linspace(lo, hi, n_bins); dd = c[1] - c[0]
    edges = np.concatenate([[c[0] - dd / 2], c + dd / 2])
    from math import erf
    cdf = np.array([0.5 * (1 + erf((e - d0) / (s0 * 2 ** 0.5))) for e in edges]); p = np.diff(cdf); p /= p.sum()
    z = alpha * p; b = (1 - alpha) * p; f = g * dt / dd
    for _ in range(int(round(T / dt))):
        ex_out = k_out * dt * z; ex_in = k_in * dt * b
        z += ex_in - ex_out; b += ex_out - ex_in
        mv = f * z; mv[-1] = 0.0; z -= mv; z[1:] += mv[:-1]          # fraction f of every bin moves one bin up
    q = z + b
    m = (q * c).sum(); s = ((q * (c - m) ** 2).sum()) ** 0.5
    return c, q / dd, m, s, dd

def moments():
    P = np.array([alpha, 1 - alpha]); A1 = P * d0; A2 = P * (d0 * d0 + s0 * s0)
    for _ in range(int(round(T / dt))):
        for A in (P, A1, A2):
            o = k_out * dt * A[0]; i = k_in * dt * A[1]; A[0] += i - o; A[1] += o - i
        D = g * dt
        A2[0] += 2 * D * A1[0] + D * D * P[0]; A1[0] += D * P[0]
    m = A1.sum(); return m, (A2.sum() - m * m) ** 0.5

c24, q24, m24, s24, dd24 = bins(27)      # spacing of the old header: 0.7 mm / 23
c96, q96, m96, s96, dd96 = bins(106)     # four times finer
mm, sm = moments()
print('reference   mean %.4f sigma %.4f  (theory sigma %.4f)' % (d_ref.mean(), d_ref.std(), (s0 * s0 + var_theory) ** 0.5))
print('bins %.4f mm mean %.4f sigma %.4f' % (dd24, m24, s24))
print('bins %.4f mm mean %.4f sigma %.4f' % (dd96, m96, s96))
print('moments      mean %.4f sigma %.4f' % (mm, sm))
print('numerical spread of the bins: %.4f mm, estimate sqrt(dd*L) = %.4f' % ((s24 ** 2 - d_ref.std() ** 2) ** 0.5, (dd24 * L) ** 0.5))
np.savez('toy.npz', d_ref=d_ref, c24=c24, q24=q24, c96=c96, q96=q96, mm=mm, sm=sm, m24=m24, s24=s24, m96=m96, s96=s96, dd24=dd24, dd96=dd96)

"""One-parameter (w = Omega/T) families of Einstein kernels used by the cell certificates.

All in T = 1 units, b_k = 2 pi k, k >= 1.

MFamily (rule NM-G, docs/theory/03_certificates.md Sec. 2):
    l(k) = m0 / (w^2 + b_k^2)          (Einstein atom at w with lambda' w^2 = m0, lambda' = m0/w^2)
    diagonal shift s(w) = beta (lam - m0/w^2)   (Lagrange multiplier for the lambda constraint;
                                                 beta = 0 and lam = None for the lambda = inf rule)
    s is concave in w (s'' = -6 beta m0 / w^4 <= 0), so s(wc+t) <= s(wc) + t s'(wc).

AFamily (rule N1mu, Sec. 4): the wE1 Einstein kernel with A fixed,
    l(k) = 2 a w / (w^2 + b_k^2),   plus the constant Coulomb block -2 mu* on n, m < Nc.
"""
from flint import arb, arb_mat
from arbcore import PI, TWO_PI


class MFamily:
    def __init__(self, m0, lam, beta):
        self.m0 = arb(m0)
        self.lam = None if lam is None else arb(lam)
        self.beta = arb(beta) if lam is not None else arb(0)
        self.has_shift = lam is not None and float(beta) != 0.0

    def lvals(self, w, kmax):
        out = [None]
        w2 = w * w
        for k in range(1, kmax + 1):
            out.append(self.m0 / (w2 + (TWO_PI * k) ** 2))
        return out

    def dlvals(self, w, kmax):
        out = [None]
        w2 = w * w
        for k in range(1, kmax + 1):
            d = w2 + (TWO_PI * k) ** 2
            out.append(-2 * self.m0 * w / (d * d))
        return out

    def d2l_abs_upper(self, w, kmax):
        out = [None]
        w2 = w * w
        for k in range(1, kmax + 1):
            b2 = (TWO_PI * k) ** 2
            d = w2 + b2
            v = self.m0 * (6 * w2 - 2 * b2) / (d * d * d)
            out.append(arb(v.abs_upper()))
        return out

    def L_total(self, w):
        """L = sum_{k>=1} m0/(w^2+b_k^2) = (m0/2)(coth(w/2)/(2w) - 1/w^2) (Lemma 1); None if
        the ball touches 0."""
        if not w.lower() > 0:
            return None
        return (self.m0 / 2) * ((w / 2).coth() / (2 * w) - 1 / (w * w))

    def T_integral(self, w, j, lj):
        """upper bound for T(j) = sum_{k>=j} l(k) by integral comparison (l decreasing in k):
        l(j) + int_j^inf m0/(w^2 + 4 pi^2 t^2) dt <= l(j) + min(m0 atan(w/2pi j)/(2 pi w), m0/(4 pi^2 j))."""
        best = (lj + self.m0 / (4 * PI * PI * j)).upper()
        if w.lower() > 0:
            u = lj + self.m0 * (w / (TWO_PI * j)).atan() / (TWO_PI * w)
            best = min(best, u.upper())
        return arb(best)

    # diagonal shift
    def shift(self, w):
        if not self.has_shift:
            return arb(0)
        return self.beta * (self.lam - self.m0 / (w * w))

    def dshift(self, w):
        if not self.has_shift:
            return arb(0)
        return 2 * self.beta * self.m0 / (w * w * w)

    def shift_sup(self, w1, w2):
        """sup of s over [w1, w2] (s increasing in w)."""
        return arb(self.shift(arb(w2)).upper())

    def const_block(self, N):
        return None

    def describe(self):
        return dict(family="M", m0=str(self.m0), lam=None if self.lam is None else str(self.lam),
                    beta=str(self.beta))


class AFamily:
    def __init__(self, a, mu=0, Nc=0):
        self.a = arb(a)
        self.mu = arb(mu)
        self.Nc = int(Nc)
        self.has_shift = False

    def lvals(self, w, kmax):
        out = [None]
        for k in range(1, kmax + 1):
            out.append(2 * self.a * w / (w * w + (TWO_PI * k) ** 2))
        return out

    def dlvals(self, w, kmax):
        out = [None]
        for k in range(1, kmax + 1):
            b2 = (TWO_PI * k) ** 2
            out.append(2 * self.a * (b2 - w * w) / (w * w + b2) ** 2)
        return out

    def d2l_abs_upper(self, w, kmax):
        out = [None]
        for k in range(1, kmax + 1):
            b2 = (TWO_PI * k) ** 2
            v = -4 * self.a * w * (3 * b2 - w * w) / (w * w + b2) ** 3
            out.append(arb(v.abs_upper()))
        return out

    def L_total(self, w):
        if not w.lower() > 0:
            return None
        return (self.a / 2) * ((w / 2).coth() - 2 / w)

    def T_integral(self, w, j, lj):
        return arb((lj + (self.a / PI) * (w / (TWO_PI * j)).atan()).upper())

    def shift(self, w):
        return arb(0)

    def dshift(self, w):
        return arb(0)

    def shift_sup(self, w1, w2):
        return arb(0)

    def const_block(self, N):
        """Coulomb term -2 mu* on the Nc x Nc leading block (symmetric block convention)."""
        if self.Nc == 0 or not self.mu > 0:
            return None
        assert self.Nc <= N
        rows = [[(-2 * self.mu if (i < self.Nc and j < self.Nc) else arb(0)) for j in range(N)]
                for i in range(N)]
        return arb_mat(rows)

    def describe(self):
        return dict(family="A", a=str(self.a), mu=str(self.mu), Nc=self.Nc)

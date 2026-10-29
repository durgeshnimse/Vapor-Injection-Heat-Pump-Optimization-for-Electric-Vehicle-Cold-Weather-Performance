"""
co2_props.py  --  R744 (CO2) thermodynamic property backend
=============================================================
Units (all public functions):  P [MPa],  T [degC],  h [kJ/kg],  s [kJ/kg-K],  rho [kg/m3]
Reference state: IIR  (h = 200 kJ/kg, s = 1.000 kJ/kg-K for saturated liquid at 0 degC)

Backends
--------
1. CoolProp (Span-Wagner reference EOS)  -> used automatically if `import CoolProp` works.
   >>> This is the recommended backend for the FINAL report numbers. <<<
2. Built-in fallback: Peng-Robinson cubic EOS for the residual part + the Span-Wagner (1996)
   ideal-gas Helmholtz function for cp0(T).  Self-contained (numpy/scipy only).
   Saturation pressures are within ~1 % of Span-Wagner; enthalpy/entropy differences are
   typically within a few kJ/kg; liquid density is less accurate (PR is known to be ~5-10 % off).

Force a backend with the environment variable  CO2_BACKEND=PR  or  CO2_BACKEND=COOLPROP.
"""
import os
import math
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq

# ----------------------------------------------------------------------------------------
# Fluid constants (Span & Wagner, 1996)
# ----------------------------------------------------------------------------------------
M_CO2 = 44.0098e-3           # kg/mol
R_U = 8.314462618            # J/mol-K
R = R_U / M_CO2              # J/kg-K  (188.92)
TC = 304.1282                # K
PC = 7.3773e6                # Pa
RHOC = 467.6                 # kg/m3
OMEGA = 0.22394
T_TRIPLE = 216.592           # K

# Span-Wagner ideal-gas part: phi0 = ln(delta)+a1+a2*tau+a3*ln(tau)+sum a_i ln(1-exp(-th_i*tau))
_A1, _A2, _A3 = 8.37304456, -3.70454304, 2.5
_AI = np.array([1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678])
_TH = np.array([3.15163, 6.11190, 6.77708, 11.32384, 27.08792])

# Peng-Robinson parameters
_KAPPA = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA ** 2
_A0 = 0.45724 * R ** 2 * TC ** 2 / PC
_B = 0.07780 * R * TC / PC
_SQ2 = math.sqrt(2.0)


def _ideal(T):
    """Return (h0 [J/kg], s0 at 1 bar [J/kg-K], cp0 [J/kg-K]) of the ideal gas (arbitrary constants)."""
    tau = TC / T
    ex = np.exp(-_TH * tau)
    phi0_tau = _A2 + _A3 / tau + np.sum(_AI * _TH * ex / (1 - ex))
    h0 = R * T * (1.0 + tau * phi0_tau)
    delta = 1e5 / (R * T * RHOC)
    phi0 = math.log(delta) + _A1 + _A2 * tau + _A3 * math.log(tau) + np.sum(_AI * np.log(1 - ex))
    s0 = R * (tau * phi0_tau - phi0)
    cp0 = R * (1.0 + _A3 + np.sum(_AI * _TH ** 2 * tau ** 2 * ex / (1 - ex) ** 2))
    return h0, s0, cp0


def _pr_ab(T):
    Tr = T / TC
    alpha = (1 + _KAPPA * (1 - math.sqrt(Tr))) ** 2
    a = _A0 * alpha
    da = -_A0 * _KAPPA * math.sqrt(alpha) / math.sqrt(T * TC)
    return a, da


def _roots(T, P):
    a, _ = _pr_ab(T)
    A = a * P / (R * T) ** 2
    B = _B * P / (R * T)
    c = [1.0, -(1 - B), A - 3 * B ** 2 - 2 * B, -(A * B - B ** 2 - B ** 3)]
    r = np.roots(c)
    r = np.real(r[np.abs(np.imag(r)) < 1e-9])
    r = r[r > B]
    return np.sort(r), A, B


def _lnphi(Z, A, B):
    L = math.log((Z + (1 + _SQ2) * B) / (Z + (1 - _SQ2) * B))
    return Z - 1 - math.log(Z - B) - A / (2 * _SQ2 * B) * L


def _Z_select(T, P, phase="auto"):
    r, A, B = _roots(T, P)
    if len(r) == 1 or phase == "auto" and len(r) == 1:
        return r[0], A, B
    if phase == "liq":
        return r[0], A, B
    if phase == "vap":
        return r[-1], A, B
    lp = [_lnphi(z, A, B) for z in (r[0], r[-1])]
    return (r[0] if lp[0] < lp[1] else r[-1]), A, B


def _resid(T, P, Z, A, B):
    a, da = _pr_ab(T)
    L = math.log((Z + (1 + _SQ2) * B) / (Z + (1 - _SQ2) * B))
    h_res = R * T * (Z - 1) + (T * da - a) / (2 * _SQ2 * _B) * L
    s_res = R * math.log(Z - B) + da / (2 * _SQ2 * _B) * L
    return h_res, s_res


def _hs_raw(T, P, phase="auto"):
    """Absolute (un-referenced) h [J/kg], s [J/kg-K], rho for T [K], P [Pa]."""
    Z, A, B = _Z_select(T, P, phase)
    h_res, s_res = _resid(T, P, Z, A, B)
    h0, s0, _ = _ideal(T)
    s = s0 - R * math.log(P / 1e5) + s_res
    return h0 + h_res, s, P / (Z * R * T)


@lru_cache(maxsize=None)
def _psat_raw(T):
    """Saturation pressure [Pa] by successive substitution on fugacity."""
    Tr = T / TC
    P = PC * math.exp(5.373 * (1 + OMEGA) * (1 - 1 / Tr))
    for _ in range(300):
        r, A, B = _roots(T, P)
        if len(r) < 3:
            # outside the 3-root window: nudge pressure toward it
            P *= 0.97 if len(r) == 1 and r[0] > 0.3 else 1.03
            continue
        dl = _lnphi(r[0], A, B) - _lnphi(r[-1], A, B)
        P_new = P * math.exp(dl)
        if abs(P_new / P - 1) < 1e-11:
            return P_new
        P = P_new
    return P


@lru_cache(maxsize=None)
def _tsat_raw(P):
    f = lambda T: _psat_raw(T) - P
    return brentq(f, 235.0, TC - 0.2, xtol=1e-9)


_REF = {}


def _ref():
    if not _REF:
        T0 = 273.15
        P0 = _psat_raw(T0)
        h, s, _ = _hs_raw(T0, P0, "liq")
        _REF["h"] = 200e3 - h
        _REF["s"] = 1000.0 - s
    return _REF["h"], _REF["s"]


# ----------------------------------------------------------------------------------------
# Public API (PR backend)
# ----------------------------------------------------------------------------------------
class _PR:
    name = "Peng-Robinson + Span-Wagner ideal part (built-in fallback)"

    @staticmethod
    @lru_cache(maxsize=None)
    def _hsr(P, T, phase):
        dh, ds = _ref()
        h, s, rho = _hs_raw(T + 273.15, P * 1e6, phase)
        return (h + dh) / 1e3, (s + ds) / 1e3, rho

    def h(self, P, T, phase="auto"):   return self._hsr(round(P, 9), round(T, 9), phase)[0]
    def s(self, P, T, phase="auto"):   return self._hsr(round(P, 9), round(T, 9), phase)[1]
    def rho(self, P, T, phase="auto"): return self._hsr(round(P, 9), round(T, 9), phase)[2]

    def Psat(self, T):  return _psat_raw(T + 273.15) / 1e6
    def Tsat(self, P):  return _tsat_raw(P * 1e6) - 273.15

    def _sat(self, P):
        T = self.Tsat(P)
        return self._hsr(round(P, 9), round(T, 9), "liq"), self._hsr(round(P, 9), round(T, 9), "vap"), T

    def hf(self, P):   return self._sat(P)[0][0]
    def hg(self, P):   return self._sat(P)[1][0]
    def sf(self, P):   return self._sat(P)[0][1]
    def sg(self, P):   return self._sat(P)[1][1]
    def rhof(self, P): return self._sat(P)[0][2]
    def rhog(self, P): return self._sat(P)[1][2]

    def T_ps(self, P, s):
        lo = (self.Tsat(P) + 1e-4) if P < PC / 1e6 else -30.0
        f = lambda T: self.s(P, T) - s
        return brentq(f, lo, 400.0, xtol=1e-9)

    def h_ps(self, P, s):
        return self.h(P, self.T_ps(P, s))

    def T_ph(self, P, h):
        lo = (self.Tsat(P) + 1e-4) if P < PC / 1e6 else -30.0
        f = lambda T: self.h(P, T) - h
        return brentq(f, lo, 400.0, xtol=1e-9)


class _CP:
    name = "CoolProp (Span-Wagner EOS, IIR reference)"

    def __init__(self):
        import CoolProp.CoolProp as CPm
        self.C = CPm
        CPm.set_reference_state("CO2", "IIR")

    def h(self, P, T, phase="auto"):   return self.C.PropsSI("H", "P", P * 1e6, "T", T + 273.15, "CO2") / 1e3
    def s(self, P, T, phase="auto"):   return self.C.PropsSI("S", "P", P * 1e6, "T", T + 273.15, "CO2") / 1e3
    def rho(self, P, T, phase="auto"): return self.C.PropsSI("D", "P", P * 1e6, "T", T + 273.15, "CO2")
    def Psat(self, T):  return self.C.PropsSI("P", "T", T + 273.15, "Q", 0, "CO2") / 1e6
    def Tsat(self, P):  return self.C.PropsSI("T", "P", P * 1e6, "Q", 0, "CO2") - 273.15
    def hf(self, P):    return self.C.PropsSI("H", "P", P * 1e6, "Q", 0, "CO2") / 1e3
    def hg(self, P):    return self.C.PropsSI("H", "P", P * 1e6, "Q", 1, "CO2") / 1e3
    def sf(self, P):    return self.C.PropsSI("S", "P", P * 1e6, "Q", 0, "CO2") / 1e3
    def sg(self, P):    return self.C.PropsSI("S", "P", P * 1e6, "Q", 1, "CO2") / 1e3
    def rhof(self, P):  return self.C.PropsSI("D", "P", P * 1e6, "Q", 0, "CO2")
    def rhog(self, P):  return self.C.PropsSI("D", "P", P * 1e6, "Q", 1, "CO2")
    def T_ps(self, P, s): return self.C.PropsSI("T", "P", P * 1e6, "S", s * 1e3, "CO2") - 273.15
    def h_ps(self, P, s): return self.C.PropsSI("H", "P", P * 1e6, "S", s * 1e3, "CO2") / 1e3
    def T_ph(self, P, h): return self.C.PropsSI("T", "P", P * 1e6, "H", h * 1e3, "CO2") - 273.15


def get_props():
    want = os.environ.get("CO2_BACKEND", "").upper()
    if want != "PR":
        try:
            return _CP()
        except Exception:
            if want == "COOLPROP":
                raise
    return _PR()


if __name__ == "__main__":
    p = get_props()
    print("Backend:", p.name)
    print("\nSaturation pressure check [MPa]  (Span-Wagner / NIST literature values)")
    ref = {-20: 1.9696, -15: 2.2906, -10: 2.6492, 0: 3.4851, 10: 4.5020, 20: 5.7291, 30: 7.2142}
    for T, Pr in ref.items():
        print(f"  T={T:4d} C   model={p.Psat(T):7.4f}   ref={Pr:7.4f}   err={100*(p.Psat(T)/Pr-1):+5.2f} %")
    print("\nIIR reference: hf(0C)=%.2f  sf(0C)=%.4f" % (p.hf(p.Psat(0)), p.sf(p.Psat(0))))
    for T in (-15, 0, 20):
        P = p.Psat(T)
        print(f"  T={T:3d}C  hf={p.hf(P):7.2f} hg={p.hg(P):7.2f} hfg={p.hg(P)-p.hf(P):6.2f}  rhof={p.rhof(P):6.1f} rhog={p.rhog(P):6.2f}")
    print("  cp0(300 K)= %.1f J/kgK (lit. 846)" % _ideal(300.0)[2])
    print("  9.5 MPa, 35 C: h=%.2f s=%.4f rho=%.1f" % (p.h(9.5, 35), p.s(9.5, 35), p.rho(9.5, 35)))

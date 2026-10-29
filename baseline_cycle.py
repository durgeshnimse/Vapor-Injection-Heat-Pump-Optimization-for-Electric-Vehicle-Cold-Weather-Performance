"""
baseline_cycle.py  --  STUDENT 1 : conventional (non-VI) transcritical R744 heat pump, Case A
==============================================================================================
State points
  1  evaporator outlet / compressor suction   (P_evap , T_evap + superheat)
  2s isentropic compressor outlet             (P_dis , s1)
  2  actual compressor discharge              (P_dis , h2 = h1 + (h2s-h1)/eta_is)
  3  gas-cooler outlet                        (P_dis , T_gc_out)
  4  EEV outlet (isenthalpic)                 (P_evap, h4 = h3)   two-phase
"""
from dataclasses import dataclass, asdict
from config import BRIEF, ASSUME
from co2_props import get_props

props = get_props()


# ----------------------------------------------------------------------------------------
# helper models
# ----------------------------------------------------------------------------------------
def vol_eff(P_dis, P_suc):
    """Volumetric efficiency vs overall pressure ratio (assumption, see report)."""
    return ASSUME["EV0"] - ASSUME["EV1"] * (P_dis / P_suc)


def cabin_load(T_amb, T_cabin=BRIEF["T_cabin"]):
    """Cabin heating load [kW] = envelope loss + ventilation (illustrative steady-state model)."""
    dT = T_cabin - T_amb
    return (ASSUME["UA_cabin"] * dT / 1000.0) + ASSUME["m_vent"] * ASSUME["cp_air"] * dT


def evap_temperature(T_amb):
    return T_amb - ASSUME["dT_evap_approach"]


# ----------------------------------------------------------------------------------------
# baseline cycle
# ----------------------------------------------------------------------------------------
def baseline(T_evap=BRIEF["T_evap"], P_dis=BRIEF["P_dis"], T_gc_out=BRIEF["T_gc_out"],
             N=BRIEF["N"], V_disp=None, eta_is=BRIEF["eta_is"], SH=ASSUME["superheat"],
             m_dot=None):
    """Solve the steady-state baseline cycle.  V_disp [m3] fixes the compressor size.
    If `m_dot` is supplied (kg/s) it overrides the swept-volume flow (used for calibration)."""
    P_e = props.Psat(T_evap)
    T1 = T_evap + SH
    h1, s1, rho1 = props.h(P_e, T1), props.s(P_e, T1), props.rho(P_e, T1)
    h2s = props.h_ps(P_dis, s1)
    h2 = h1 + (h2s - h1) / eta_is
    T2 = props.T_ph(P_dis, h2)
    T2s = props.T_ps(P_dis, s1)
    h3, s3 = props.h(P_dis, T_gc_out), props.s(P_dis, T_gc_out)
    h4 = h3
    T4 = props.Tsat(P_e)
    x4 = (h4 - props.hf(P_e)) / (props.hg(P_e) - props.hf(P_e))
    s4 = props.sf(P_e) + x4 * (props.sg(P_e) - props.sf(P_e))
    PR = P_dis / P_e
    eta_v = vol_eff(P_dis, P_e)
    if m_dot is None:
        m_dot = rho1 * V_disp * (N / 60.0) * eta_v
    W = m_dot * (h2 - h1)
    Q_h = m_dot * (h2 - h3)
    Q_e = m_dot * (h1 - h4)
    res = dict(
        T_evap=T_evap, P_evap=P_e, P_dis=P_dis, PR=PR, N=N, eta_v=eta_v, eta_is=eta_is,
        m_dot=m_dot, rho1=rho1,
        h1=h1, s1=s1, T1=T1, h2s=h2s, T2s=T2s, h2=h2, T2=T2, h3=h3, s3=s3, T3=T_gc_out,
        h4=h4, T4=T4, x4=x4, s4=s4,
        q_h=h2 - h3, q_e=h1 - h4, w=h2 - h1,
        Q_h=Q_h, Q_e=Q_e, W=W, COP_h=Q_h / W, COP_c=Q_e / W,
        T_dis=T2, energy_balance_err=Q_h - Q_e - W,
        feasible=(T2 <= BRIEF["T_dis_max"]) and (P_dis <= BRIEF["P_max"]),
    )
    return res


def calibrate_displacement(Q_target=BRIEF["Q_load"], T_evap=BRIEF["T_evap"]):
    """Size the compressor so that Case A delivers exactly Q_target at the design point
    (baseline speed, baseline pressure).  Q_h is linear in V_disp, so one pass suffices."""
    ref = baseline(T_evap=T_evap, V_disp=1.0e-5)
    V = 1.0e-5 * Q_target / ref["Q_h"]
    return V


# ----------------------------------------------------------------------------------------
# gas-cooler air side consistency check
# ----------------------------------------------------------------------------------------
def gas_cooler_air_check(res):
    dT_air = res["Q_h"] / (ASSUME["m_air_gc"] * ASSUME["cp_air"])
    T_sup = ASSUME["T_air_gc_in"] + dT_air
    return dict(dT_air=dT_air, T_air_out=T_sup,
                approach_cold_end=res["T3"] - ASSUME["T_air_gc_in"],
                T_air_out_ok=T_sup < res["T2"])


if __name__ == "__main__":
    print("Backend:", props.name)
    V = calibrate_displacement()
    print("Calibrated displacement: %.2f cc/rev" % (V * 1e6))
    r = baseline(V_disp=V)
    for k in ("P_evap", "PR", "eta_v", "m_dot", "h1", "h2s", "h2", "h3", "h4", "x4", "T2s", "T2",
              "q_h", "w", "Q_h", "Q_e", "W", "COP_h", "energy_balance_err"):
        print(f"  {k:20s}{r[k]:12.4f}")
    print(gas_cooler_air_check(r))

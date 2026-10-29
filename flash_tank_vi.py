"""
flash_tank_vi.py  --  STUDENT 2 : flash-tank vapour-injection (FT-VI) R744 heat pump, Case B
============================================================================================
State points (numbering used in the report)
  1  evaporator outlet / suction            (P_evap , T_evap+SH)
  2s isentropic end of stage 1              (P_inj  , s1)
  2  actual end of stage 1                  (P_inj  , h2 = h1+(h2s-h1)/eta_is)
  3  after mixing with injected vapour      (P_inj  , h3 = [(1-x) h2 + x h7] )   [per kg total flow]
  4s isentropic end of stage 2              (P_dis  , s3)
  4  discharge                              (P_dis  , h4 = h3+(h4s-h3)/eta_is)
  5  gas-cooler outlet                      (P_dis  , T_gc_out)
  6  EEV-1 outlet = flash-tank inlet        (P_inj  , h6 = h5)  two-phase, quality x
  7  flash-tank vapour (sat. vapour)        (P_inj)   -> injected
  8  flash-tank liquid (sat. liquid)        (P_inj)
  9  EEV-2 outlet                           (P_evap , h9 = h8)

Mass balance on flash tank (per kg of gas-cooler flow, m_tot = 1):
    x = (h6 - h_f)/(h_g - h_f)      m_inj = x m_tot      m_suc = m_liq = (1-x) m_tot
    injection ratio  xi = m_inj / m_suc = x/(1-x)
Energy balance (adiabatic tank):  m_tot h6 = m_inj h7 + m_liq h8   (checked numerically)
"""
import math
from config import BRIEF, ASSUME
from co2_props import get_props
from baseline_cycle import vol_eff, baseline, calibrate_displacement

props = get_props()


def flash_tank_balance(P_inj, h_in, m_tot=1.0):
    """Mass & energy balance of an adiabatic flash tank fed by a two-phase stream of enthalpy h_in.
    Returns a dict with vapour quality, flows and a closure residual."""
    hf, hg = props.hf(P_inj), props.hg(P_inj)
    x = (h_in - hf) / (hg - hf)
    m_vap, m_liq = x * m_tot, (1 - x) * m_tot
    res_energy = m_tot * h_in - (m_vap * hg + m_liq * hf)
    return dict(x=x, m_vap=m_vap, m_liq=m_liq, hf=hf, hg=hg, xi=m_vap / m_liq if m_liq > 0 else float("inf"),
                energy_residual=res_energy,
                rho_f=props.rhof(P_inj), rho_g=props.rhog(P_inj), T_tank=props.Tsat(P_inj))


def flash_tank_size(m_vap, m_liq, rho_f, rho_g):
    """Indicative tank sizing: Souders-Brown vapour velocity limit + liquid hold-up time."""
    v_max = ASSUME["K_souders"] * math.sqrt((rho_f - rho_g) / rho_g)      # m/s
    A = m_vap / (rho_g * v_max)                                           # m2
    D = math.sqrt(4 * A / math.pi)
    V_liq = m_liq * ASSUME["t_hold"] / rho_f                              # m3
    return dict(v_max=v_max, D_mm=D * 1e3, V_liq_cc=V_liq * 1e6)


def ftvi(V_disp, T_evap=BRIEF["T_evap"], P_dis=BRIEF["P_dis"], P_inj=BRIEF["P_inj"],
         T_gc_out=BRIEF["T_gc_out"], N=BRIEF["N"], eta_is=BRIEF["eta_is"], SH=ASSUME["superheat"],
         f_vi=1.0):
    """f_vi (<=1) derates the suction flow of the VI compressor relative to the non-VI one
    (ideal thermodynamic balance: f_vi = 1).  COP is independent of f_vi; capacity scales with it."""
    P_e = props.Psat(T_evap)
    if not (P_e < P_inj < P_dis):
        raise ValueError("Constraint violated: need P_evap < P_inj < P_dis")
    T1 = T_evap + SH
    h1, s1, rho1 = props.h(P_e, T1), props.s(P_e, T1), props.rho(P_e, T1)
    # high side
    h5, s5 = props.h(P_dis, T_gc_out), props.s(P_dis, T_gc_out)
    h6 = h5
    tank = flash_tank_balance(P_inj, h6)
    x = tank["x"]
    h7, h8 = tank["hg"], tank["hf"]
    h9 = h8
    T9 = props.Tsat(P_e)
    x9 = (h9 - props.hf(P_e)) / (props.hg(P_e) - props.hf(P_e))
    # stage 1
    h2s = props.h_ps(P_inj, s1)
    h2 = h1 + (h2s - h1) / eta_is
    T2 = props.T_ph(P_inj, h2)
    # mixing (per kg of total flow): (1-x) kg from stage 1, x kg of saturated vapour
    h3 = (1 - x) * h2 + x * h7
    T3 = props.T_ph(P_inj, h3)
    s3 = props.s(P_inj, T3)
    # stage 2
    h4s = props.h_ps(P_dis, s3)
    h4 = h3 + (h4s - h3) / eta_is
    T4 = props.T_ph(P_dis, h4)
    T4s = props.T_ps(P_dis, s3)
    # flows (swept-volume model: suction flow fixed by the compressor displacement)
    eta_v = vol_eff(P_dis, P_e)
    m_suc = rho1 * V_disp * (N / 60.0) * eta_v * f_vi
    m_tot = m_suc / (1 - x)
    m_inj = x * m_tot
    # specific work (per kg total flow) and rates
    w1 = (1 - x) * (h2 - h1)
    w2 = (h4 - h3)
    w_tot = w1 + w2
    W = m_tot * w_tot
    q_h = h4 - h5
    Q_h = m_tot * q_h
    Q_e = m_suc * (h1 - h9)
    size = flash_tank_size(m_inj, m_suc, tank["rho_f"], tank["rho_g"])
    return dict(
        T_evap=T_evap, P_evap=P_e, P_inj=P_inj, P_dis=P_dis, PR=P_dis / P_e,
        PR1=P_inj / P_e, PR2=P_dis / P_inj, N=N, eta_v=eta_v, eta_is=eta_is,
        h1=h1, T1=T1, s1=s1, h2s=h2s, h2=h2, T2=T2, h3=h3, T3=T3, s3=s3, h4s=h4s, T4s=T4s, h4=h4, T4=T4,
        h5=h5, T5=T_gc_out, h6=h6, T6=tank["T_tank"], x=x, h7=h7, h8=h8, h9=h9, T9=T9, x9=x9,
        m_suc=m_suc, m_inj=m_inj, m_tot=m_tot, xi=tank["xi"], inj_frac=x,
        w1=w1, w2=w2, w_tot=w_tot, q_h=q_h,
        Q_h=Q_h, Q_e=Q_e, W=W, COP_h=Q_h / W, T_dis=T4,
        tank_energy_residual=tank["energy_residual"],
        energy_balance_err=Q_h - Q_e - W,
        tank_D_mm=size["D_mm"], tank_Vliq_cc=size["V_liq_cc"], tank_vmax=size["v_max"],
        feasible=(T4 <= BRIEF["T_dis_max"]) and (P_dis <= BRIEF["P_max"]),
    )


def sweep_Pinj(V_disp, P_list, **kw):
    out = []
    for p in P_list:
        try:
            out.append(ftvi(V_disp, P_inj=p, **kw))
        except ValueError:
            pass
    return out


def optimum_Pinj_geometric(P_evap, P_dis):
    """Rule-of-thumb intermediate pressure  sqrt(P_evap * P_dis)."""
    return math.sqrt(P_evap * P_dis)


if __name__ == "__main__":
    V = calibrate_displacement()
    A = baseline(V_disp=V)
    B = ftvi(V)
    print("Backend:", props.name)
    for k in ("x", "xi", "m_suc", "m_inj", "m_tot", "h2", "T2", "h3", "T3", "h4", "T4", "q_h", "w_tot",
              "Q_h", "W", "COP_h", "Q_e", "tank_energy_residual", "energy_balance_err", "tank_D_mm"):
        print(f"  {k:22s}{B[k]:12.4f}")
    print("Baseline:  Q=%.3f W=%.3f COP=%.3f Tdis=%.1f" % (A["Q_h"], A["W"], A["COP_h"], A["T2"]))
    print("Geometric-mean Pinj = %.2f MPa" % optimum_Pinj_geometric(B["P_evap"], B["P_dis"]))

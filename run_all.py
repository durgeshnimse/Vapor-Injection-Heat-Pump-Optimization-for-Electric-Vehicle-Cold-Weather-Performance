"""
run_all.py -- runs every calculation for Students 1 & 2, writes CSV/JSON results and figures.
Usage:  python run_all.py            (uses CoolProp if installed, otherwise the built-in PR backend)
"""
import json, csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

from config import BRIEF, ASSUME, SENS_AMBIENTS
from co2_props import get_props
import baseline_cycle as bc
import flash_tank_vi as fv

props = get_props()
OUT = os.path.dirname(os.path.abspath(__file__))
FIG, RES = os.path.join(OUT, "figures"), os.path.join(OUT, "results")
os.makedirs(FIG, exist_ok=True); os.makedirs(RES, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False})
CA, CB, CR = "#c0392b", "#1f6fb4", "#555555"


def save_csv(name, rows):
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


R = {}                       # everything the PDF builder needs
R["backend"] = props.name

# ======================================================================================
# 0. Compressor sizing + design point
# ======================================================================================
V = bc.calibrate_displacement()
R["V_disp_cc"] = V * 1e6
A = bc.baseline(V_disp=V)
B = fv.ftvi(V)
R["A"], R["B"] = A, B
R["air_check"] = {k: float(v) if not isinstance(v, (bool, np.bool_)) else bool(v)
                  for k, v in bc.gas_cooler_air_check(A).items()}

# cabin load vs ambient
R["cabin_load"] = {str(T): bc.cabin_load(T) for T in (-5, -10, -15)}

# ======================================================================================
# 1. Student 1 parametric studies (Case A)
# ======================================================================================
Pd = np.linspace(8.0, 12.0, 41)
S1_pd = [bc.baseline(P_dis=p, V_disp=V) for p in Pd]
save_csv("A_Pdis_sweep.csv", [dict(P_dis=r["P_dis"], Q_h=r["Q_h"], W=r["W"], COP=r["COP_h"], T_dis=r["T_dis"], m_dot=r["m_dot"]) for r in S1_pd])
cop_opt = max(S1_pd, key=lambda r: r["COP_h"])
feas = [r for r in S1_pd if r["T_dis"] <= BRIEF["T_dis_max"]]
R["A_cop_opt_unconstrained"] = dict(P=cop_opt["P_dis"], COP=cop_opt["COP_h"], Tdis=cop_opt["T_dis"], Q=cop_opt["Q_h"])
R["A_Pdis_max_Tlimit"] = max(r["P_dis"] for r in feas)
R["A_cop_opt_constrained"] = max(feas, key=lambda r: r["COP_h"])["P_dis"]
R["A_Pdis_table"] = [S1_pd[i] for i in range(0, 41, 5)]

Ns = np.linspace(BRIEF["N_min"], BRIEF["N_max"], 23)
S1_N = [bc.baseline(N=n, V_disp=V) for n in Ns]
save_csv("A_speed_sweep.csv", [dict(N=r["N"], Q_h=r["Q_h"], W=r["W"], COP=r["COP_h"], T_dis=r["T_dis"], m_dot=r["m_dot"]) for r in S1_N])
R["A_speed_table"] = [bc.baseline(N=n, V_disp=V) for n in (2500, 3000, 4000, 5000, 6000, 8000)]
R["A_N_for_load"] = BRIEF["Q_load"] / A["Q_h"] * BRIEF["N"]

Te = np.linspace(-25, -5, 21)
S1_Te = [bc.baseline(T_evap=t, V_disp=V) for t in Te]
save_csv("A_Tevap_sweep.csv", [dict(T_evap=r["T_evap"], Q_h=r["Q_h"], W=r["W"], COP=r["COP_h"], T_dis=r["T_dis"]) for r in S1_Te])

# ambient sensitivity, Case A and B (same compressor, same speed, Pinj = 4.5 MPa and geometric mean)
amb_rows = []
for Ta in SENS_AMBIENTS:
    Tev = bc.evap_temperature(Ta)
    a = bc.baseline(T_evap=Tev, V_disp=V)
    b = fv.ftvi(V, T_evap=Tev)
    Pgm = fv.optimum_Pinj_geometric(a["P_evap"], BRIEF["P_dis"])
    bg = fv.ftvi(V, T_evap=Tev, P_inj=Pgm)
    amb_rows.append(dict(T_amb=Ta, load=bc.cabin_load(Ta), A=a, B=b, Bg=bg, Pgm=Pgm))
R["amb"] = amb_rows
save_csv("A_vs_B_ambient.csv", [dict(T_amb=r["T_amb"], T_evap=r["A"]["T_evap"], cabin_load=r["load"],
        A_Q=r["A"]["Q_h"], A_W=r["A"]["W"], A_COP=r["A"]["COP_h"], A_Tdis=r["A"]["T_dis"],
        B_Q=r["B"]["Q_h"], B_W=r["B"]["W"], B_COP=r["B"]["COP_h"], B_Tdis=r["B"]["T_dis"],
        B_x=r["B"]["x"], B_xi=r["B"]["xi"]) for r in amb_rows])

# ======================================================================================
# 2. Student 2 parametric study:  P_inj sweep
# ======================================================================================
Pi = np.round(np.arange(3.5, 6.0001, 0.1), 3)
S2 = {}
for Ta in SENS_AMBIENTS:
    Tev = bc.evap_temperature(Ta)
    S2[Ta] = fv.sweep_Pinj(V, Pi, T_evap=Tev)
    save_csv(f"B_Pinj_sweep_amb{int(Ta)}.csv", [dict(P_inj=r["P_inj"], x=r["x"], xi=r["xi"], m_suc=r["m_suc"], m_inj=r["m_inj"],
             m_tot=r["m_tot"], Q_h=r["Q_h"], W=r["W"], COP=r["COP_h"], T_dis=r["T_dis"], T2=r["T2"], T3=r["T3"]) for r in S2[Ta]])
S2d = S2[-10.0]
R["S2_design"] = S2d
iC = max(range(len(S2d)), key=lambda i: S2d[i]["COP_h"])
iQ = max(range(len(S2d)), key=lambda i: S2d[i]["Q_h"])
R["S2_best_COP"] = S2d[iC]; R["S2_best_Q"] = S2d[iQ]
# continuous refinement of COP-optimal Pinj
f = lambda p: -fv.ftvi(V, P_inj=p)["COP_h"]
rr = minimize_scalar(f, bounds=(3.5, 6.0), method="bounded", options=dict(xatol=1e-4))
R["S2_Pinj_opt_COP"] = float(rr.x); R["S2_COP_opt"] = float(-rr.fun)
R["S2_Pgm"] = fv.optimum_Pinj_geometric(B["P_evap"], B["P_dis"])
R["S2_table"] = [r for r in S2d if abs(r["P_inj"] * 4 - round(r["P_inj"] * 4)) < 1e-6]   # every 0.25 MPa

# derating sensitivity
R["fvi_table"] = [dict(f=f_, **{k: v for k, v in fv.ftvi(V, f_vi=f_).items() if k in ("Q_h", "W", "COP_h", "T_dis", "m_tot")})
                  for f_ in (1.0, 0.95, 0.90, 0.85, 0.80)]

# P_inj x P_dis map
Pdm = np.linspace(8.0, 12.0, 17); Pim = np.linspace(3.5, 6.0, 26)
COPm = np.zeros((len(Pdm), len(Pim))); TDm = np.zeros_like(COPm); QM = np.zeros_like(COPm)
for i, pd_ in enumerate(Pdm):
    for j, pi_ in enumerate(Pim):
        r = fv.ftvi(V, P_dis=pd_, P_inj=pi_)
        COPm[i, j], TDm[i, j], QM[i, j] = r["COP_h"], r["T_dis"], r["Q_h"]

# ======================================================================================
# 3. Hand-off data for Students 3 & 4
# ======================================================================================
hand = dict(
    backend=props.name, V_disp_cc=V * 1e6,
    operating_conditions=dict(BRIEF=BRIEF, ASSUMPTIONS=ASSUME),
    CaseA={k: v for k, v in A.items()}, CaseB={k: v for k, v in B.items()},
)
json.dump(hand, open(os.path.join(RES, "handoff_students_3_4.json"), "w"), indent=2, default=float)
stateA = [("1", "Evap. outlet / suction", A["P_evap"], A["T1"], A["h1"], A["s1"], "-"),
          ("2s", "Isentropic discharge", A["P_dis"], A["T2s"], A["h2s"], A["s1"], "-"),
          ("2", "Actual discharge", A["P_dis"], A["T2"], A["h2"], "", "-"),
          ("3", "Gas-cooler outlet", A["P_dis"], A["T3"], A["h3"], A["s3"], "-"),
          ("4", "EEV outlet", A["P_evap"], A["T4"], A["h4"], A["s4"], f"{A['x4']:.4f}")]
stateB = [("1", "Evap. outlet / suction", B["P_evap"], B["T1"], B["h1"], B["s1"], "-"),
          ("2s", "Stage-1 isentropic", B["P_inj"], "", B["h2s"], B["s1"], "-"),
          ("2", "Stage-1 actual", B["P_inj"], B["T2"], B["h2"], "", "-"),
          ("3", "After mixing", B["P_inj"], B["T3"], B["h3"], B["s3"], "-"),
          ("4s", "Stage-2 isentropic", B["P_dis"], B["T4s"], B["h4s"], B["s3"], "-"),
          ("4", "Discharge", B["P_dis"], B["T4"], B["h4"], "", "-"),
          ("5", "Gas-cooler outlet", B["P_dis"], B["T5"], B["h5"], props.s(B["P_dis"], B["T5"]), "-"),
          ("6", "Flash-tank inlet", B["P_inj"], B["T6"], B["h6"], "", f"{B['x']:.4f}"),
          ("7", "Tank vapour (sat.)", B["P_inj"], B["T6"], B["h7"], props.sg(B["P_inj"]), "1"),
          ("8", "Tank liquid (sat.)", B["P_inj"], B["T6"], B["h8"], props.sf(B["P_inj"]), "0"),
          ("9", "EEV-2 outlet", B["P_evap"], B["T9"], B["h9"], "", f"{B['x9']:.4f}")]
R["stateA"], R["stateB"] = stateA, stateB
for nm, st in (("stateA", stateA), ("stateB", stateB)):
    save_csv(f"{nm}_table.csv", [dict(state=a, description=b, P_MPa=c, T_C=d, h_kJkg=e, s_kJkgK=f_, quality=g) for a, b, c, d, e, f_, g in st])

# ======================================================================================
# FIGURES
# ======================================================================================
def dome(ax, kind):
    Ts = np.linspace(-40, 30.5, 120)
    Ts = [t for t in Ts if props.Psat(t) < 7.34]
    if kind == "ph":
        ax.plot([props.hf(props.Psat(t)) for t in Ts], [props.Psat(t) for t in Ts], color="k", lw=1)
        ax.plot([props.hg(props.Psat(t)) for t in Ts], [props.Psat(t) for t in Ts], color="k", lw=1)
    else:
        ax.plot([props.sf(props.Psat(t)) for t in Ts], Ts, color="k", lw=1)
        ax.plot([props.sg(props.Psat(t)) for t in Ts], Ts, color="k", lw=1)


def path_ph(P, h_a, h_b, n=40):
    return np.linspace(h_a, h_b, n), np.full(n, P)


# --- Fig 1: P-h ---------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
a0, a1 = ax
dome(a0, "ph"); dome(a1, "ph")
a0.plot([A["h1"], A["h2"], A["h3"], A["h4"], A["h1"]], [A["P_evap"], A["P_dis"], A["P_dis"], A["P_evap"], A["P_evap"]], "-o", color=CA, ms=4)
for n, hh, pp in (("1", A["h1"], A["P_evap"]), ("2", A["h2"], A["P_dis"]), ("3", A["h3"], A["P_dis"]), ("4", A["h4"], A["P_evap"])):
    a0.annotate(n, (hh, pp), textcoords="offset points", xytext=(5, 5), fontsize=9, weight="bold")
a0.set_title("Case A: conventional heat pump (T_amb = -10 °C)")
pts = [(B["h1"], B["P_evap"]), (B["h2"], B["P_inj"]), (B["h3"], B["P_inj"]), (B["h4"], B["P_dis"]), (B["h5"], B["P_dis"]),
       (B["h6"], B["P_inj"]), (B["h8"], B["P_inj"]), (B["h9"], B["P_evap"]), (B["h1"], B["P_evap"])]
a1.plot([p[0] for p in pts], [p[1] for p in pts], "-o", color=CB, ms=4)
a1.plot([B["h6"], B["h7"]], [B["P_inj"], B["P_inj"]], "--", color=CB)
a1.plot([B["h7"], B["h3"]], [B["P_inj"], B["P_inj"]], ":", color=CB)
for n, hh, pp in (("1", B["h1"], B["P_evap"]), ("2", B["h2"], B["P_inj"]), ("3", B["h3"], B["P_inj"]), ("4", B["h4"], B["P_dis"]),
                  ("5", B["h5"], B["P_dis"]), ("6", B["h6"], B["P_inj"]), ("7", B["h7"], B["P_inj"]), ("8", B["h8"], B["P_inj"]), ("9", B["h9"], B["P_evap"])):
    a1.annotate(n, (hh, pp), textcoords="offset points", xytext=(5, 5), fontsize=9, weight="bold")
a1.set_title("Case B: flash-tank vapour injection (P_inj = 4.5 MPa)")
for a_ in ax:
    a_.set_yscale("log"); a_.set_xlabel("Specific enthalpy h [kJ/kg]"); a_.set_xlim(150, 600); a_.set_ylim(1.5, 14)
    a_.set_yticks([2, 3, 4, 5, 7, 10, 12]); a_.set_yticklabels(["2", "3", "4", "5", "7", "10", "12"])
a0.set_ylabel("Pressure P [MPa]")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_ph.png")); plt.close()

# --- Fig 2: T-s ---------------------------------------------------------------------
def ts_line(pts):
    """pts: list of (P1,h1,P2,h2) straight lines in (log P, h) → (s, T) arrays"""
    S, T = [], []
    eps = 1e-3          # kJ/kg – tolerance around the saturation enthalpies
    for (p1, h1, p2, h2) in pts:
        for f_ in np.linspace(0, 1, 30):
            p = np.exp(np.log(p1) + f_ * (np.log(p2) - np.log(p1)))
            h = h1 + f_ * (h2 - h1)

            if p < 7.3:
                hf_ = props.hf(p)
                hg_ = props.hg(p)
                if (hf_ - eps) <= h <= (hg_ + eps):
                    # treat as two-phase (or extremely close to it)
                    Tt = props.Tsat(p)
                    q  = np.clip((h - hf_) / (hg_ - hf_), 0.0, 1.0)
                    s  = props.sf(p) + q * (props.sg(p) - props.sf(p))
                    S.append(s)
                    T.append(Tt)
                    continue

            # single-phase (or supercritical)
            Tt = props.T_ph(p, h)
            try:
                s = props.s(p, Tt)
            except ValueError:
                # last-resort fallback if CoolProp still complains
                hf_ = props.hf(p)
                hg_ = props.hg(p)
                s   = props.sf(p) if abs(h - hf_) < abs(h - hg_) else props.sg(p)
            S.append(s)
            T.append(Tt)
    return S, T

fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
for a_ in ax: dome(a_, "ts")
Pe, Pd_, Pi_ = A["P_evap"], A["P_dis"], B["P_inj"]
sA, tA = ts_line([(Pe, A["h1"], Pd_, A["h2"]), (Pd_, A["h2"], Pd_, A["h3"]), (Pd_, A["h3"], Pe, A["h4"]), (Pe, A["h4"], Pe, A["h1"])])
ax[0].plot(sA, tA, color=CA); ax[0].set_title("Case A: T–s diagram")
pB = [(Pe, B["h1"], Pi_, B["h2"]), (Pi_, B["h3"], Pd_, B["h4"]), (Pd_, B["h4"], Pd_, B["h5"]), (Pd_, B["h5"], Pi_, B["h6"]),
      (Pi_, B["h8"], Pe, B["h9"]), (Pe, B["h9"], Pe, B["h1"])]
sB, tB = ts_line(pB)
ax[1].plot(sB, tB, ".", ms=2, color=CB)
for (p1, h1, p2, h2) in [(Pi_, B["h6"], Pi_, B["h7"]), (Pi_, B["h7"], Pi_, B["h3"]), (Pi_, B["h2"], Pi_, B["h3"])]:
    s_, t_ = ts_line([(p1, h1, p2, h2)]); ax[1].plot(s_, t_, "--", color=CB, lw=1)
ax[1].set_title("Case B: T–s diagram (dots = process paths, dashed = flash / mixing)")
for a_ in ax:
    a_.set_xlabel("Specific entropy s [kJ/kg·K]"); a_.set_xlim(0.7, 2.1); a_.set_ylim(-40, 140)
ax[0].set_ylabel("Temperature T [°C]")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_ts.png")); plt.close()

# --- Fig 3: schematic ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9.5, 4.6)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 5)
def box(x, y, w, h, t, c="#eaf2fb"):
    ax.add_patch(plt.Rectangle((x, y), w, h, fc=c, ec="k", lw=1.2)); ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=8.5)
def arrow(x1, y1, x2, y2, c="k", ls="-"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", color=c, lw=1.4, ls=ls))
def lab(x, y, t): ax.text(x, y, t, fontsize=9, weight="bold", color="#8e1b0e", bbox=dict(fc="w", ec="none", pad=0.5))
box(3.6, 3.6, 2.4, 1.0, "Two-stage compressor\n(vapour injection port)", "#fdebd0")
box(7.6, 3.4, 1.9, 1.4, "Gas cooler\n(cabin heater)", "#fadbd8")
box(0.2, 3.4, 1.9, 1.4, "Evaporator\n(ambient air)", "#d6eaf8")
box(4.1, 1.6, 1.4, 1.2, "Flash tank\n(P_inj)", "#d5f5e3")
box(7.4, 0.5, 1.2, 0.7, "EEV-1")
box(1.6, 0.5, 1.2, 0.7, "EEV-2")
arrow(2.1, 4.1, 3.6, 4.1); lab(2.55, 4.25, "1")
arrow(6.0, 4.1, 7.6, 4.1); lab(6.6, 4.25, "4")
arrow(8.55, 3.4, 8.0, 1.2); lab(8.55, 2.3, "5")
arrow(7.4, 0.85, 5.5, 1.9); lab(6.5, 1.5, "6")
arrow(4.8, 2.8, 4.8, 3.6, c=CB); lab(4.95, 3.1, "7 (vapour)")
arrow(4.3, 1.6, 2.8, 0.9); lab(3.1, 1.6, "8 (liquid)")
arrow(1.6, 0.85, 1.15, 3.4); lab(0.85, 2.0, "9")
ax.text(5.0, 0.1, "R744 flash-tank vapour-injection (FT-VI) cycle — Case B", ha="center", fontsize=9, style="italic")
plt.savefig(os.path.join(FIG, "fig_schematic.png"), bbox_inches="tight"); plt.close()

# --- Fig 4: baseline P_dis sweep ----------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
x = [r["P_dis"] for r in S1_pd]
ax[0].plot(x, [r["COP_h"] for r in S1_pd], color=CA); ax[0].set_ylabel("Heating COP [-]")
ax[1].plot(x, [r["Q_h"] for r in S1_pd], color=CA); ax[1].set_ylabel("Heating capacity [kW]")
ax[1].axhline(5.0, color=CR, ls="--", lw=1); ax[1].text(8.1, 5.03, "5.0 kW load", fontsize=8)
ax[2].plot(x, [r["T_dis"] for r in S1_pd], color=CA); ax[2].axhline(130, color="k", ls="--", lw=1)
ax[2].text(8.1, 131, "130 °C limit", fontsize=8); ax[2].set_ylabel("Discharge temperature [°C]")
for a_ in ax:
    a_.axvline(9.5, color=CR, ls=":", lw=1); a_.set_xlabel("High-side pressure [MPa]")
plt.suptitle("Case A — effect of high-side pressure (T_evap = −15 °C, 4000 rpm, T_gc,out = 35 °C)", fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_A_Pdis.png")); plt.close()

# --- Fig 5: baseline speed + evap temp ----------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
ax[0].plot([r["N"] for r in S1_N], [r["Q_h"] for r in S1_N], color=CA, label="Q_h")
ax[0].plot([r["N"] for r in S1_N], [r["W"] for r in S1_N], color=CR, label="W")
ax[0].axhline(5, ls="--", color="k", lw=1); ax[0].set_xlabel("Compressor speed [rpm]"); ax[0].set_ylabel("Q_h, W [kW]"); ax[0].legend()
ax[0].set_title("Capacity and power vs speed")
ax[1].plot([r["T_evap"] for r in S1_Te], [r["COP_h"] for r in S1_Te], color=CA)
ax[1].set_xlabel("Evaporating temperature [°C]"); ax[1].set_ylabel("Heating COP [-]"); ax[1].set_title("COP vs T_evap")
ax[2].plot([r["T_evap"] for r in S1_Te], [r["T_dis"] for r in S1_Te], color=CA); ax[2].axhline(130, ls="--", color="k", lw=1)
ax[2].set_xlabel("Evaporating temperature [°C]"); ax[2].set_ylabel("T_dis [°C]"); ax[2].set_title("Discharge T vs T_evap")
ax2 = ax[1].twinx(); ax2.plot([r["T_evap"] for r in S1_Te], [r["Q_h"] for r in S1_Te], color=CB, ls="--"); ax2.set_ylabel("Q_h [kW] (dashed)")
ax2.grid(False)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_A_speed_Te.png")); plt.close()

# --- Fig 6: P_inj sweep, design ambient ---------------------------------------------
fig, ax = plt.subplots(2, 3, figsize=(10, 6))
xp = [r["P_inj"] for r in S2d]
def mark(a_):
    a_.axvline(4.5, color=CR, ls=":", lw=1); a_.axvline(R["S2_Pgm"], color="g", ls=":", lw=1); a_.set_xlabel("Intermediate pressure P_inj [MPa]")
ax[0, 0].plot(xp, [r["x"] for r in S2d], color=CB); ax[0, 0].set_ylabel("Flash-tank vapour fraction x [-]"); ax[0, 0].set_title("(a) Vapour quality after EEV-1")
ax[0, 1].plot(xp, [r["xi"] for r in S2d], color=CB); ax[0, 1].set_ylabel("Injection ratio ξ = m_inj/m_suc [-]"); ax[0, 1].set_title("(b) Injection ratio")
ax[0, 2].plot(xp, [r["m_inj"] * 1000 for r in S2d], color=CB, label="m_inj"); ax[0, 2].plot(xp, [r["m_tot"] * 1000 for r in S2d], color=CR, label="m_tot (gas cooler)")
ax[0, 2].plot(xp, [r["m_suc"] * 1000 for r in S2d], color="g", label="m_suc"); ax[0, 2].legend(fontsize=7); ax[0, 2].set_ylabel("Mass flow [g/s]"); ax[0, 2].set_title("(c) Flow rates")
ax[1, 0].plot(xp, [r["Q_h"] for r in S2d], color=CB); ax[1, 0].axhline(A["Q_h"], color=CA, ls="--"); ax[1, 0].text(3.55, A["Q_h"] + 0.05, "Case A", color=CA, fontsize=8)
ax[1, 0].set_ylabel("Heating capacity [kW]"); ax[1, 0].set_title("(d) Heating capacity")
ax[1, 1].plot(xp, [r["COP_h"] for r in S2d], color=CB); ax[1, 1].axhline(A["COP_h"], color=CA, ls="--"); ax[1, 1].text(3.55, A["COP_h"] + 0.01, "Case A", color=CA, fontsize=8)
ax[1, 1].plot([R["S2_Pinj_opt_COP"]], [R["S2_COP_opt"]], "k*", ms=9); ax[1, 1].set_ylabel("Heating COP [-]"); ax[1, 1].set_title("(e) Heating COP (★ = optimum)")
ax[1, 2].plot(xp, [r["T_dis"] for r in S2d], color=CB, label="T_dis (4)"); ax[1, 2].plot(xp, [r["T3"] for r in S2d], color="g", label="T after mixing (3)")
ax[1, 2].axhline(A["T_dis"], color=CA, ls="--"); ax[1, 2].text(3.55, A["T_dis"] + 1, "Case A T_dis", color=CA, fontsize=8)
ax[1, 2].axhline(130, color="k", ls="--", lw=1); ax[1, 2].legend(fontsize=7, loc="center right"); ax[1, 2].set_ylabel("Temperature [°C]"); ax[1, 2].set_title("(f) Discharge temperature")
for a_ in ax.flat: mark(a_)
plt.suptitle("Case B — effect of flash-tank pressure (T_amb = −10 °C, 4000 rpm, P_dis = 9.5 MPa).  Dotted: 4.5 MPa (red), √(P_e·P_d) (green)", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_B_Pinj.png")); plt.close()

# --- Fig 7: P_inj sweep for three ambients -------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
cols = {-5.0: "#27ae60", -10.0: CB, -15.0: "#8e44ad"}
for Ta, rows in S2.items():
    ax[0].plot([r["P_inj"] for r in rows], [r["Q_h"] for r in rows], color=cols[Ta], label=f"{Ta:.0f} °C")
    ax[1].plot([r["P_inj"] for r in rows], [r["COP_h"] for r in rows], color=cols[Ta])
    ax[2].plot([r["P_inj"] for r in rows], [r["T_dis"] for r in rows], color=cols[Ta])
for Ta, rr_ in zip(SENS_AMBIENTS, amb_rows):
    ax[0].axhline(rr_["A"]["Q_h"], color=cols[Ta], ls=":", lw=1); ax[1].axhline(rr_["A"]["COP_h"], color=cols[Ta], ls=":", lw=1)
    ax[2].axhline(rr_["A"]["T_dis"], color=cols[Ta], ls=":", lw=1)
ax[2].axhline(130, color="k", ls="--", lw=1)
for a_, yl in zip(ax, ("Heating capacity [kW]", "Heating COP [-]", "Discharge temperature [°C]")):
    a_.set_xlabel("P_inj [MPa]"); a_.set_ylabel(yl)
ax[0].legend(title="T_amb (dotted = Case A)", fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_B_Pinj_amb.png")); plt.close()

# --- Fig 8: A vs B bars -------------------------------------------------------------
fig, ax = plt.subplots(1, 4, figsize=(11, 3.2))
lab_ = [f"{r['T_amb']:.0f} °C" for r in amb_rows]; xx = np.arange(3); w = 0.38
for a_, key, yl in zip(ax, ("Q_h", "W", "COP_h", "T_dis"), ("Heating capacity [kW]", "Compressor power [kW]", "Heating COP [-]", "Discharge T [°C]")):
    va = [r["A"][key] for r in amb_rows]; vb = [r["B"][key] for r in amb_rows]
    a_.bar(xx - w / 2, va, w, color=CA, label="A: baseline"); a_.bar(xx + w / 2, vb, w, color=CB, label="B: FT-VI")
    a_.set_xticks(xx); a_.set_xticklabels(lab_); a_.set_ylabel(yl); a_.set_xlabel("Ambient"); a_.grid(axis="x", visible=False)
    if key == "T_dis": a_.axhline(130, color="k", ls="--", lw=1)
    if key == "Q_h":
        a_.plot(xx, [r["load"] for r in amb_rows], "k_", ms=22, mew=2, label="cabin load")
ax[0].legend(fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_A_vs_B.png")); plt.close()

# --- Fig 9: maps --------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 3.9))
c1 = ax[0].contourf(Pim, Pdm, COPm, 18, cmap="viridis"); plt.colorbar(c1, ax=ax[0], label="Heating COP [-]")
c2 = ax[0].contour(Pim, Pdm, TDm, levels=[130], colors="r", linewidths=2)
ax[0].clabel(c2, fmt="T_dis = 130 °C", fontsize=8)
ax[0].set_xlabel("P_inj [MPa]"); ax[0].set_ylabel("P_dis [MPa]"); ax[0].set_title("FT-VI COP map (red = T_dis limit; to the right/below is feasible)")
c3 = ax[1].contourf(Pim, Pdm, QM, 18, cmap="plasma"); plt.colorbar(c3, ax=ax[1], label="Heating capacity [kW]")
ax[1].set_xlabel("P_inj [MPa]"); ax[1].set_ylabel("P_dis [MPa]"); ax[1].set_title("FT-VI heating capacity map (4000 rpm)")
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_B_map.png")); plt.close()

# --- Fig 10: property validation ----------------------------------------------------
fig, ax = plt.subplots(figsize=(4.8, 3.4))
refP = {-20: 1.9696, -15: 2.2906, -10: 2.6492, 0: 3.4851, 10: 4.5020, 20: 5.7291, 30: 7.2142}
tt = np.linspace(-25, 30, 50)
ax.plot(tt, [props.Psat(t) for t in tt], color=CB, label="model"); ax.plot(list(refP), list(refP.values()), "ko", ms=4, label="literature (Span–Wagner)")
ax.set_xlabel("Saturation temperature [°C]"); ax.set_ylabel("Saturation pressure [MPa]"); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG, "fig_props_check.png")); plt.close()
R["psat_check"] = [(T, props.Psat(T), P) for T, P in refP.items()]

# store sweep data used in text
R["S2_by_amb_best"] = {str(Ta): max(rows, key=lambda r: r["COP_h"]) for Ta, rows in S2.items()}
import pickle
pickle.dump(R, open(os.path.join(RES, "all_results.pkl"), "wb"))

# console summary
print("Backend :", props.name)
print("V_disp  : %.2f cc" % R["V_disp_cc"])
print("Case A  : Q=%.3f W=%.3f COP=%.3f Tdis=%.1f" % (A["Q_h"], A["W"], A["COP_h"], A["T_dis"]))
print("Case B  : Q=%.3f W=%.3f COP=%.3f Tdis=%.1f x=%.3f xi=%.3f" % (B["Q_h"], B["W"], B["COP_h"], B["T_dis"], B["x"], B["xi"]))
print("A: COP-opt Pdis unconstrained", R["A_cop_opt_unconstrained"], " Pdis max for Tdis<=130:", R["A_Pdis_max_Tlimit"])
print("B: Pinj opt (COP) = %.3f MPa, COP = %.4f ; geometric mean = %.3f" % (R["S2_Pinj_opt_COP"], R["S2_COP_opt"], R["S2_Pgm"]))
for r in amb_rows:
    print("Tamb %.0f: load %.2f | A Q %.2f COP %.2f Td %.1f | B Q %.2f COP %.2f Td %.1f | Bgm(%.2f) COP %.2f"
          % (r["T_amb"], r["load"], r["A"]["Q_h"], r["A"]["COP_h"], r["A"]["T_dis"], r["B"]["Q_h"], r["B"]["COP_h"], r["B"]["T_dis"], r["Pgm"], r["Bg"]["COP_h"]))

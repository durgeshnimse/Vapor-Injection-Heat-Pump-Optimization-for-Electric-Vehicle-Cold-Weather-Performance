"""Common input data (project brief) and modelling assumptions."""

# ----------------------------- Given in the project brief -------------------------------
BRIEF = dict(
    T_amb=-10.0,          # degC   design ambient
    T_cabin=20.0,         # degC   cabin target
    Q_load=5.0,           # kW     cabin heating load
    T_evap=-15.0,         # degC   baseline evaporating temperature
    T_gc_out=35.0,        # degC   gas cooler outlet
    P_dis=9.5,            # MPa    high-side pressure
    P_max=12.0,           # MPa    hard limit
    N=4000.0,             # rpm    baseline speed
    N_min=2500.0, N_max=8000.0,
    eta_is=0.70,          # -      isentropic efficiency
    P_inj=4.5,            # MPa    intermediate pressure
    P_inj_range=(3.5, 6.0),
    eev_open=50.0,        # %      (used by Students 3/4)
    T_dis_max=130.0,      # degC   hard limit
)

# ----------------------------- Assumptions added by us (stated in report) ---------------
ASSUME = dict(
    dT_evap_approach=5.0,   # K   evaporating temp = ambient - 5 K  (so -10 degC amb -> -15 degC evap)
    superheat=5.0,          # K   evaporator-outlet (suction) superheat
    no_pressure_drops=True, # -   no line / heat-exchanger pressure drops
    adiabatic_compressor=True,
    # volumetric efficiency  eta_v = EV0 - EV1 * (P_dis/P_suc)   (generic CO2 scroll/rotary trend)
    EV0=0.97, EV1=0.03,
    # cabin envelope used only to show that a 5.0 kW design load is physically consistent
    UA_cabin=125.0,         # W/K   (glazing + body shell)
    m_vent=0.040,           # kg/s  fresh-air ventilation flow
    cp_air=1.006,           # kJ/kg-K
    m_air_gc=0.15,          # kg/s  air flow through gas cooler
    T_air_gc_in=20.0,       # degC  air entering gas cooler (recirculated cabin air)
    # flash-tank sizing (indicative)
    K_souders=0.06,         # m/s   Souders-Brown coefficient for vapour/liquid separation
    t_hold=10.0,            # s     liquid hold-up time
)

SENS_AMBIENTS = [-5.0, -10.0, -15.0]

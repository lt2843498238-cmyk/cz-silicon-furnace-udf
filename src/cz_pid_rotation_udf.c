/*===========================================================================
 *  CZ Single Crystal Silicon Furnace -- PID + Crystal Rotation UDF  v3.0
 *  Fluent 24.2 / Double Precision (3ddp) / Windows / Parallel
 *
 *  ┌─────────────────────────────────────────────────────────────────┐
 *  │  DESIGN PRINCIPLES                                              │
 *  │                                                                 │
 *  │  1. Power state: static g_P, NOT from UDM read-back.            │
 *  │  2. No cross-process broadcast; PRF_GRSUM1 keeps all identical. │
 *  │  3. No host_to_node_real macros (broken in 3ddp).               │
 *  │  4. UDM[0] = heat source density, written by write_udm().       │
 *  │  5. pid_state.pid persists ALL state (PID + rotation).          │
 *  │  6. Crystal rotation via Carman-Kozeny momentum source on       │
 *  │     X and Z directions. Enabled only after |e| < 10K for 50    │
 *  │     consecutive PID updates. Suspended if |e| > 50K.           │
 *  └─────────────────────────────────────────────────────────────────┘
 *
 *  Incremental PID:
 *    dP(n) = Kp*[e(n)-e(n-1)] + Ki*e(n) + Kd*[e(n)-2e(n-1)+e(n-2)]
 *    P(n)  = P(n-1) + dP(n)
 *
 *  Crystal Rotation Source (X-momentum example):
 *    S_x = -C_ROT * (1-f_l)^2 / (f_l^3 + eps) * (v_x - v_rot_x)
 *    v_rot_x = -omega * z,   v_rot_z = +omega * x
 *    Only active when g_rot_on == 1.
 *
 *  [SETUP]
 *  1. Define -> User-Defined -> Memory -> 1 -> OK
 *  2. Compile & Load UDF
 *  3. Function Hooks: Execute at End -> heater_pid_control
 *  4. Zone 4508 -> Source Terms -> Energy  -> heater_main_source
 *     Zone 4505 -> Source Terms -> Energy  -> heater_bottom_source
 *     Zone ???? -> Source Terms -> X-Mom   -> crystal_xmom_source
 *     Zone ???? -> Source Terms -> Z-Mom   -> crystal_zmom_source
 *     (???? = silicon melt+crystal zone ID, set ID_SILICON below)
 *  5. Solution Initialize -> Calculate
 *
 *  Resume : open case+dat.h5, Load UDF, Calculate
 *  Fresh  : delete pid_state.pid, Load UDF, Calculate
 *  Reset  : Execute On Demand -> reset_pid
 *===========================================================================*/

#include "udf.h"
#include <math.h>
#include <stdio.h>

/*---------------------------------------------------------------------------
 *  Zone / Mesh parameters
 *---------------------------------------------------------------------------*/
#define ID_INTERFACE      4
#define ID_HEATER_MAIN    4508
#define ID_HEATER_BOTTOM  4505
#define ID_SILICON        4592    /* silicon melt + crystal zone */
#define V_HEATER_MAIN     0.0235083163
#define V_HEATER_BOTTOM   0.0132345022
#define R_MIN             0.1255
#define R_MAX             0.1380

/*---------------------------------------------------------------------------
 *  UDM
 *---------------------------------------------------------------------------*/
#define UDM_QDOT   0

/*---------------------------------------------------------------------------
 *  PID control parameters
 *---------------------------------------------------------------------------*/
#define T_SET          1685.0
#define P_START        250000.0
#define P_LIMIT_STEP   5000.0
#define P_TOTAL_MIN    50000.0
#define RATIO_MAIN     0.8

/*---------------------------------------------------------------------------
 *  PID zone thresholds [K] with hysteresis
 *---------------------------------------------------------------------------*/
#define E_Z1   100.0
#define E_Z2    20.0
#define E_Z3     5.0
#define E_HYST   1.0

/*---------------------------------------------------------------------------
 *  D-term low-pass filter
 *---------------------------------------------------------------------------*/
#define ALPHA_D  0.3

/*---------------------------------------------------------------------------
 *  Crystal rotation parameters
 *---------------------------------------------------------------------------*/
#define ROT_OMEGA      0.9948    /* 9.5 rpm = 9.5*2*pi/60 rad/s          */
#define ROT_C          1.0e5    /* Carman-Kozeny constant (=Mushy Zone)  */
#define ROT_EPS        1.0e-3   /* prevent division by zero              */
#define ROT_ENABLE_AE  10.0     /* |e| threshold to start counting      */
#define ROT_CONSEC     50       /* consecutive stable PID updates needed */
#define ROT_SUSPEND_AE 50.0     /* emergency suspend threshold           */

/*---------------------------------------------------------------------------
 *  Probe points for rotation velocity monitoring
 *  Probe A: near crystal center (r ~ 0.04m)
 *  Probe B: near crystal edge   (r ~ 0.10m)
 *  Crystal: Y range 1.2~2.5m, diameter 251mm (radius 0.1255m)
 *  probe_y = 1.8m (mid-crystal, well above interface, pure solid)
 *---------------------------------------------------------------------------*/
#define PROBE_A_R   0.04
#define PROBE_B_R   0.10
#define PROBE_Y     1.8
#define PROBE_TOL_Y 0.05
#define PROBE_TOL_R 1.0

/*---------------------------------------------------------------------------
 *  Persistence file & console
 *---------------------------------------------------------------------------*/
#define STATE_FILE  "pid_state.pid"
#define BOX_W       66

/*===========================================================================
 *  Static globals  (each process maintains its own copy; all identical)
 *===========================================================================*/
/* --- PID state --- */
static real g_P      = P_START;
static real g_ep1    = 0.0;
static real g_ep2    = 0.0;
static real g_dDf    = 0.0;
static int  g_ucnt   = 0;
static int  g_pzone  = -1;
static int  g_first  = 1;

static real g_Th[5]  = {T_SET, T_SET, T_SET, T_SET, T_SET};
static int  g_Tcnt   = 0;

/* --- Rotation state --- */
static real g_ae       = 999.0;  /* current |error|, read by DEFINE_SOURCE */
static int  g_rot_on   = 0;     /* 0=off, 1=on, 2=suspended              */
static int  g_rot_cnt  = 0;     /* consecutive-stable counter             */


/*===========================================================================
 *  PID zone / gain helpers  (unchanged from v2.1)
 *===========================================================================*/
static int zone_of(real ae, int prev_zone)
{
    real enter_z0 = (real)E_Z1;
    real enter_z1 = (real)E_Z2;
    real enter_z2 = (real)E_Z3;
    real leave_z1 = enter_z0 + (real)E_HYST;
    real leave_z2 = enter_z1 + (real)E_HYST;
    real leave_z3 = enter_z2 + (real)E_HYST;

    if (prev_zone < 0 || prev_zone > 3)
    {
        if      (ae > enter_z0) return 0;
        else if (ae > enter_z1) return 1;
        else if (ae > enter_z2) return 2;
        else                    return 3;
    }

    switch (prev_zone)
    {
        case 0:  return (ae <= enter_z0) ? 1 : 0;
        case 1:  if (ae > leave_z1) return 0;
                 if (ae <= enter_z1) return 2;
                 return 1;
        case 2:  if (ae > leave_z2) return 1;
                 if (ae <= enter_z2) return 3;
                 return 2;
        case 3:  if (ae > leave_z3) return 2;
                 return 3;
        default: return 3;
    }
}

static void gains(int z, real *kp, real *ki, real *kd, int *fr)
{
    switch(z) {
        case 0: *kp=400;  *ki=15;  *kd=30;  *fr=10; break;
        case 1: *kp=500;  *ki=20;  *kd=60;  *fr=10; break;
        case 2: *kp=400;  *ki=12;  *kd=80;  *fr=10; break;
        default:*kp=150;  *ki=5;   *kd=40;  *fr=5;  break;
    }
}

static const char* zname(int z)
{
    switch(z) {
        case 0: return "Coarse  |e|>100K";
        case 1: return "Trans   20~100K ";
        case 2: return "Fine     5~20K  ";
        default:return "Prec      <5K   ";
    }
}

static const char* zstatus(real ae)
{
    if      (ae < 0.5)  return "** CONVERGED **";
    else if (ae < 5.0)  return "   Precision   ";
    else if (ae < 20.0) return "   Fine Tune   ";
    else                return "   Heating Up  ";
}


/*===========================================================================
 *  Console helpers  (unchanged)
 *===========================================================================*/
static void pb(int t)
{
    int i;
    char c = t ? '-' : '=';
    Message0("  +");
    for (i = 0; i < BOX_W; i++) Message0("%c", c);
    Message0("+\n");
}

static void pl(const char *s)
{
    int n = 0;
    const char *p = s;
    int pad, i;
    while (*p++) n++;
    Message0("  | %s", s);
    pad = BOX_W - 1 - n;
    for (i = 0; i < pad; i++) Message0(" ");
    Message0(" |\n");
}


/*===========================================================================
 *  Persistence  (v3.0: 15 fields = v2.1's 12 + g_ae + g_rot_on + g_rot_cnt)
 *
 *  Backward compatible: reads v3.0 (15), v2.1 (12), v1.0 (11) formats.
 *  All processes read and write (same strategy as v2.1, proven safe).
 *===========================================================================*/
static void pid_save(void)
{
    FILE *fp = fopen(STATE_FILE, "w");
    if (!fp) { Message0("  [PID] WARNING: cannot write %s\n", STATE_FILE); return; }
    fprintf(fp, "%.15le %.15le %.15le %.15le %d %d "
                "%.6le %.6le %.6le %.6le %.6le %d "
                "%.15le %d %d\n",
            g_P, g_ep1, g_ep2, g_dDf, g_ucnt, g_pzone,
            g_Th[0], g_Th[1], g_Th[2], g_Th[3], g_Th[4], g_Tcnt,
            g_ae, g_rot_on, g_rot_cnt);
    fclose(fp);
}

static int pid_load(void)
{
    FILE *fp;
    real p, e1, e2, ddf, ae_val, t0, t1, t2, t3, t4;
    int uc, pz, tc, ron, rcnt, n;

    fp = fopen(STATE_FILE, "r");
    if (!fp) return 0;

    /* Try v3.0: 15 fields */
    n = fscanf(fp, "%le %le %le %le %d %d %le %le %le %le %le %d %le %d %d",
               &p, &e1, &e2, &ddf, &uc, &pz,
               &t0, &t1, &t2, &t3, &t4, &tc,
               &ae_val, &ron, &rcnt);
    if (n == 15) {
        fclose(fp);
        g_P = p;  g_ep1 = e1;  g_ep2 = e2;  g_dDf = ddf;
        g_ucnt = uc;  g_pzone = pz;
        g_Th[0] = t0;  g_Th[1] = t1;  g_Th[2] = t2;
        g_Th[3] = t3;  g_Th[4] = t4;  g_Tcnt = tc;
        g_ae = ae_val;  g_rot_on = ron;  g_rot_cnt = rcnt;
        return 1;
    }

    /* Try v2.1: 12 fields */
    rewind(fp);
    n = fscanf(fp, "%le %le %le %le %d %d %le %le %le %le %le %d",
               &p, &e1, &e2, &ddf, &uc, &pz,
               &t0, &t1, &t2, &t3, &t4, &tc);
    if (n == 12) {
        fclose(fp);
        g_P = p;  g_ep1 = e1;  g_ep2 = e2;  g_dDf = ddf;
        g_ucnt = uc;  g_pzone = pz;
        g_Th[0] = t0;  g_Th[1] = t1;  g_Th[2] = t2;
        g_Th[3] = t3;  g_Th[4] = t4;  g_Tcnt = tc;
        g_ae = 999.0;  g_rot_on = 0;  g_rot_cnt = 0;
        Message0("  [PID] Loaded v2.1 state (rotation state reset)\n");
        return 1;
    }

    /* Try v1.0: 11 fields */
    rewind(fp);
    n = fscanf(fp, "%le %le %le %d %d %le %le %le %le %le %d",
               &p, &e1, &e2, &uc, &pz,
               &t0, &t1, &t2, &t3, &t4, &tc);
    fclose(fp);
    if (n == 11) {
        g_P = p;  g_ep1 = e1;  g_ep2 = e2;  g_dDf = 0.0;
        g_ucnt = uc;  g_pzone = pz;
        g_Th[0] = t0;  g_Th[1] = t1;  g_Th[2] = t2;
        g_Th[3] = t3;  g_Th[4] = t4;  g_Tcnt = tc;
        g_ae = 999.0;  g_rot_on = 0;  g_rot_cnt = 0;
        Message0("  [PID] Loaded v1.0 state (rotation state reset)\n");
        return 1;
    }

    return 0;
}


/*===========================================================================
 *  write_udm  (unchanged from v2.1)
 *===========================================================================*/
static void write_udm(Thread *tm, Thread *tb, real power)
{
    cell_t c;
    real dm = (power * (real)RATIO_MAIN)         / (real)V_HEATER_MAIN;
    real db = (power * (1.0 - (real)RATIO_MAIN)) / (real)V_HEATER_BOTTOM;

    if (NNULLP(tm))
    {
        begin_c_loop(c, tm)
        {
            C_UDMI(c, tm, UDM_QDOT) = dm;
        }
        end_c_loop(c, tm)
    }

    if (NNULLP(tb))
    {
        begin_c_loop(c, tb)
        {
            C_UDMI(c, tb, UDM_QDOT) = db;
        }
        end_c_loop(c, tb)
    }
}


/*===========================================================================
 *  DEFINE_ON_DEMAND: reset_pid  (v3.0: also resets rotation)
 *===========================================================================*/
DEFINE_ON_DEMAND(reset_pid)
{
    Domain *d  = Get_Domain(1);
    Thread *tm = Lookup_Thread(d, ID_HEATER_MAIN);
    Thread *tb = Lookup_Thread(d, ID_HEATER_BOTTOM);
    int i;

    g_P = P_START;
    g_ep1 = 0.0;  g_ep2 = 0.0;  g_dDf = 0.0;
    g_ucnt = 0;   g_pzone = -1;  g_first = 0;
    g_Tcnt = 0;
    for (i = 0; i < 5; i++) g_Th[i] = (real)T_SET;

    g_ae = 999.0;  g_rot_on = 0;  g_rot_cnt = 0;

    write_udm(tm, tb, g_P);
    Message0("\n[PID-RESET] P=%.0fW  Rotation=OFF  All state reset.\n\n",
             (real)P_START);
}


/*===========================================================================
 *  DEFINE_EXECUTE_AT_END: heater_pid_control
 *===========================================================================*/
DEFINE_EXECUTE_AT_END(heater_pid_control)
{
    Domain *domain  = Get_Domain(1);
    Thread *t_iface = Lookup_Thread(domain, ID_INTERFACE);
    Thread *t_main  = Lookup_Thread(domain, ID_HEATER_MAIN);
    Thread *t_bot   = Lookup_Thread(domain, ID_HEATER_BOTTOM);

    face_t f;
    real   x[ND_ND], T_sum, A_sum, r, area, T_avg;
    real   e_now, de1, de2, ae;
    real   kp, ki, kd, tP, tI, tD_raw, tD, dPr, dPo, newP;
    int    fcnt, i, clipped, freq, zone;
    char   buf[128];

    /* ==============================================================
     *  FIRST CALL: restore state
     * ============================================================== */
    if (g_first)
    {
        int ok = pid_load();
        write_udm(t_main, t_bot, g_P);
        g_first = 0;
        pb(0);
        if (ok) {
            pl(" CZ-PID v3.0  Restored  (Resume Mode)");
            sprintf(buf, " P=%.0fW  e_prev=%+.4fK  Upd=%d  Rot=%s",
                    g_P, g_ep1, g_ucnt,
                    g_rot_on == 1 ? "ON" :
                    (g_rot_on == 2 ? "SUSPENDED" : "OFF"));
        } else {
            pl(" CZ-PID v3.0  Fresh Start");
            sprintf(buf, " T_set=%.1fK  P_start=%.0fW  r:[%.4f,%.4f]m",
                    (real)T_SET, (real)P_START, (real)R_MIN, (real)R_MAX);
        }
        pl(buf);
        pb(2);
        Message0("\n");
    }

    /* ==============================================================
     *  TEMPERATURE SAMPLING  (unchanged)
     * ============================================================== */
    T_sum = 0.0;  A_sum = 0.0;  fcnt = 0;

    if (NNULLP(t_iface))
    {
        begin_f_loop(f, t_iface)
        {
            if (PRINCIPAL_FACE_P(f, t_iface))
            {
                F_CENTROID(x, f, t_iface);
                r = sqrt(x[0]*x[0] + x[2]*x[2]);
                if (r >= (real)R_MIN && r <= (real)R_MAX)
                {
                    area   = NV_MAG(F_AREA_CACHE(f, t_iface));
                    T_sum += F_T(f, t_iface) * area;
                    A_sum += area;
                    fcnt++;
                }
            }
        }
        end_f_loop(f, t_iface)
    }

    T_sum = PRF_GRSUM1(T_sum);
    A_sum = PRF_GRSUM1(A_sum);
    fcnt  = (int)PRF_GRSUM1((real)fcnt);
    T_avg = (A_sum > 1.0e-20) ? (T_sum / A_sum) : (real)T_SET;

    /* ==============================================================
     *  FREQUENCY GATE
     * ============================================================== */
    if (N_ITER == 0 || N_ITER % 5 != 0) return;

    e_now = (real)T_SET - T_avg;
    ae    = fabs(e_now);
    g_ae  = ae;   /* update for DEFINE_SOURCE to read */

    zone  = zone_of(ae, g_pzone);
    gains(zone, &kp, &ki, &kd, &freq);
    if (N_ITER % freq != 0) return;

    if (A_sum < 1.0e-20) {
        Message0("  [PID-WARN] face area=0! Check R_MIN/R_MAX.\n");
        return;
    }

    /* ==============================================================
     *  FIRST PID STEP
     * ============================================================== */
    if (g_ucnt == 0) {
        g_ep1 = e_now;
        g_ep2 = e_now;
    }

    /* ==============================================================
     *  INCREMENTAL PID  (unchanged from v2.1)
     * ============================================================== */
    de1 = e_now - g_ep1;
    de2 = g_ep1 - g_ep2;

    tP     = kp * de1;
    tI     = ki * e_now;
    tD_raw = kd * (de1 - de2);
    tD     = (real)ALPHA_D * tD_raw + (1.0 - (real)ALPHA_D) * g_dDf;
    dPr    = tP + tI + tD;

    dPo     = dPr;
    clipped = 0;
    if (dPo >  (real)P_LIMIT_STEP) { dPo =  (real)P_LIMIT_STEP; clipped = 1; }
    if (dPo < -(real)P_LIMIT_STEP) { dPo = -(real)P_LIMIT_STEP; clipped = 1; }

    newP = g_P + dPo;
    if (newP < (real)P_TOTAL_MIN) newP = (real)P_TOTAL_MIN;

    g_P     = newP;
    g_ep2   = g_ep1;
    g_ep1   = e_now;
    g_dDf   = tD;
    g_ucnt++;
    g_pzone = zone;

    for (i = 4; i > 0; i--) g_Th[i] = g_Th[i-1];
    g_Th[0] = T_avg;
    if (g_Tcnt < 5) g_Tcnt++;

    /* ==============================================================
     *  ROTATION ENABLE / SUSPEND LOGIC
     *
     *  State machine:
     *    OFF (0) --[|e|<10K x50]--> ON (1)
     *    ON  (1) --[|e|>50K]------> SUSPENDED (2)
     *    SUSPENDED (2) --[|e|<10K x50]--> ON (1)
     *    OFF/SUSPENDED: counter resets if |e| >= 10K
     * ============================================================== */
    if (g_rot_on == 1)
    {
        if (ae > (real)ROT_SUSPEND_AE) {
            g_rot_on  = 2;
            g_rot_cnt = 0;
        }
    }
    else
    {
        if (ae < (real)ROT_ENABLE_AE) {
            g_rot_cnt++;
            if (g_rot_cnt >= ROT_CONSEC) {
                g_rot_on  = 1;
                g_rot_cnt = ROT_CONSEC;
            }
        } else {
            g_rot_cnt = 0;
        }
    }

    /* ==============================================================
     *  SAVE & WRITE HEAT SOURCE
     * ============================================================== */
    pid_save();
    write_udm(t_main, t_bot, newP);

    /* ==============================================================
     *  PROBE ROTATION VELOCITY
     *
     *  When rotation is ON, scan silicon zone for two probe points
     *  to verify solid cells are actually rotating.
     *
     *  Tangential velocity: v_theta = (-v_x*z + v_z*x) / r
     *  Target:              v_theta_target = omega * r
     *
     *  Runs only at PID update frequency, negligible cost.
     *  Message0 prints from node_zero only -- we get node_zero's
     *  closest cells, which is approximate but sufficient.
     * ============================================================== */
    {
        Thread *t_si;
        cell_t  cc;
        real    xc[ND_ND], rc, dist, vt, fl_val;
        real    bestA_d = 1.0e10, bestA_vt = 0.0, bestA_r = 0.0, bestA_fl = 1.0;
        real    bestB_d = 1.0e10, bestB_vt = 0.0, bestB_r = 0.0, bestB_fl = 1.0;

        t_si = Lookup_Thread(domain, ID_SILICON);

        if (g_rot_on == 1 && NNULLP(t_si))
        {
            begin_c_loop(cc, t_si)
            {
                C_CENTROID(xc, cc, t_si);

                /* Only look near probe height */
                if (fabs(xc[1] - (real)PROBE_Y) < (real)PROBE_TOL_Y)
                {
                    rc = sqrt(xc[0]*xc[0] + xc[2]*xc[2]);
                    fl_val = C_VOF(cc, t_si);

                    /* Tangential velocity */
                    if (rc > 1.0e-6)
                        vt = (-C_U(cc, t_si) * xc[2]
                              + C_W(cc, t_si) * xc[0]) / rc;
                    else
                        vt = 0.0;

                    /* Probe A: closest to r = PROBE_A_R */
                    dist = fabs(rc - (real)PROBE_A_R);
                    if (dist < bestA_d) {
                        bestA_d = dist;  bestA_r = rc;
                        bestA_vt = vt;   bestA_fl = fl_val;
                    }

                    /* Probe B: closest to r = PROBE_B_R */
                    dist = fabs(rc - (real)PROBE_B_R);
                    if (dist < bestB_d) {
                        bestB_d = dist;  bestB_r = rc;
                        bestB_vt = vt;   bestB_fl = fl_val;
                    }
                }
            }
            end_c_loop(cc, t_si)
        }

        /* ==============================================================
         *  CONSOLE OUTPUT
         * ============================================================== */
        Message0("\n");
        pb(0);
        sprintf(buf, " [CZ-PID]  Iter:%-6d  Upd#:%-4d  %-16s  %-15s",
                N_ITER, g_ucnt, zname(zone), zstatus(ae));
        pl(buf);
        pb(1);

        sprintf(buf, " [TEMP]   Faces:%-3d  T_meas:%10.4fK  Target:%7.1fK",
                fcnt, T_avg, (real)T_SET);
        pl(buf);
        sprintf(buf, "          Error:%+10.4fK", e_now);
        pl(buf);
        pb(1);

        sprintf(buf, " [POWER]  Total:%9.0fW  Main(80%%):%9.0fW",
                newP, newP * (real)RATIO_MAIN);
        pl(buf);
        sprintf(buf, "          Bot(20%%):%9.0fW",
                newP * (1.0 - (real)RATIO_MAIN));
        pl(buf);
        pb(1);

        sprintf(buf, " [PID]    Kp:%-5.0f Ki:%-5.0f Kd:%-5.0f  Zone:%d",
                kp, ki, kd, zone);
        pl(buf);
        sprintf(buf, "          P-term:%+10.2fW  I-term:%+10.2fW", tP, tI);
        pl(buf);
        sprintf(buf, "          D-raw: %+10.2fW  D-filt:%+10.2fW", tD_raw, tD);
        pl(buf);
        if (clipped)
            sprintf(buf,
                "          dP_raw:%+10.2fW  dP_out:%+10.2fW [CLIPPED]",
                dPr, dPo);
        else
            sprintf(buf,
                "          dP_raw:%+10.2fW  dP_out:%+10.2fW", dPr, dPo);
        pl(buf);
        pb(1);

        /* --- Rotation status --- */
        if (g_rot_on == 1) {
            sprintf(buf, " [ROT]    ON   omega=%.4f rad/s  (%.1f rpm)",
                    (real)ROT_OMEGA, 9.5);
            pl(buf);
            if (bestA_d < (real)PROBE_TOL_R) {
                sprintf(buf,
                    "          ProbeA r=%.3fm v_th=%+.5f tgt=%+.5f fl=%.2f",
                    bestA_r, bestA_vt,
                    (real)ROT_OMEGA * bestA_r, bestA_fl);
                pl(buf);
            }
            if (bestB_d < (real)PROBE_TOL_R) {
                sprintf(buf,
                    "          ProbeB r=%.3fm v_th=%+.5f tgt=%+.5f fl=%.2f",
                    bestB_r, bestB_vt,
                    (real)ROT_OMEGA * bestB_r, bestB_fl);
                pl(buf);
            }
        } else if (g_rot_on == 2) {
            sprintf(buf, " [ROT]    ** SUSPENDED **  |e|=%.1fK > %.0fK",
                    ae, (real)ROT_SUSPEND_AE);
            pl(buf);
            sprintf(buf, "          Count: %d/%d  (re-stabilizing...)",
                    g_rot_cnt, ROT_CONSEC);
            pl(buf);
        } else {
            sprintf(buf,
                " [ROT]    OFF  Count: %d/%d  |e|=%.1fK  (need <%.0fK)",
                g_rot_cnt, ROT_CONSEC, ae, (real)ROT_ENABLE_AE);
            pl(buf);
        }
        pb(1);

        /* --- Temperature history --- */
        {
            char h[120];
            int p2 = 0;
            p2 += sprintf(h + p2, " [T-HIST] ");
            for (i = 0; i < 5; i++)
                p2 += sprintf(h + p2,
                    i < g_Tcnt ? (i == 0 ? "*%8.3f " : " %8.3f ")
                               : "  ------- ",
                    i < g_Tcnt ? g_Th[i] : 0.0);
            p2 += sprintf(h + p2, "K");
            pl(h);
        }
        {
            char d[120];
            int p2 = 0;
            p2 += sprintf(d + p2, " [dT/stp] ");
            if (g_Tcnt >= 2) {
                p2 += sprintf(d + p2, "          ");
                for (i = 0; i < g_Tcnt - 1; i++)
                    p2 += sprintf(d + p2, " %+8.4f", g_Th[i] - g_Th[i+1]);
            } else {
                p2 += sprintf(d + p2, " (accumulating...)");
            }
            pl(d);
        }
        pb(2);
        Message0("\n");
    }
}


/*===========================================================================
 *  DEFINE_SOURCE: heater_main_source / heater_bottom_source
 *  (unchanged from v2.1)
 *===========================================================================*/
DEFINE_SOURCE(heater_main_source, c, t, dS, eqn)
{
    dS[eqn] = 0.0;
    return C_UDMI(c, t, UDM_QDOT);
}

DEFINE_SOURCE(heater_bottom_source, c, t, dS, eqn)
{
    dS[eqn] = 0.0;
    return C_UDMI(c, t, UDM_QDOT);
}


/*===========================================================================
 *  DEFINE_SOURCE: crystal_xmom_source / crystal_zmom_source
 *
 *  Crystal rotation via Carman-Kozeny momentum source.
 *  Applied to the silicon cell zone (melt + crystal in one zone).
 *
 *  Physics:
 *    S = -C * (1-fl)^2 / (fl^3+eps) * (v - v_rot)
 *    v_rot_x = -omega * z    (rotation about Y axis)
 *    v_rot_z = +omega * x
 *
 *  Behavior by liquid fraction:
 *    fl = 1 (liquid):  (1-fl)^2 = 0  -> S = 0      (free flow)
 *    fl = 0 (solid):   huge force    -> v = v_rot   (locked to rotation)
 *    0 < fl < 1:       smooth transition (mushy zone)
 *
 *  When g_rot_on != 1: returns 0 (no rotation applied).
 *
 *  dS[eqn] provides implicit coupling for numerical stability.
 *
 *  SETUP: Assign to X-Momentum / Z-Momentum source of silicon zone.
 *===========================================================================*/
DEFINE_SOURCE(crystal_xmom_source, c, t, dS, eqn)
{
    real fl, sw, xc[ND_ND], v_rot_x;

    if (g_rot_on != 1) {
        dS[eqn] = 0.0;
        return 0.0;
    }

    fl = C_VOF(c, t);
    sw = (1.0 - fl) * (1.0 - fl) / (fl * fl * fl + (real)ROT_EPS);

    C_CENTROID(xc, c, t);
    v_rot_x = -(real)ROT_OMEGA * xc[2];

    dS[eqn] = -(real)ROT_C * sw;
    return -(real)ROT_C * sw * (C_U(c, t) - v_rot_x);
}

DEFINE_SOURCE(crystal_zmom_source, c, t, dS, eqn)
{
    real fl, sw, xc[ND_ND], v_rot_z;

    if (g_rot_on != 1) {
        dS[eqn] = 0.0;
        return 0.0;
    }

    fl = C_VOF(c, t);
    sw = (1.0 - fl) * (1.0 - fl) / (fl * fl * fl + (real)ROT_EPS);

    C_CENTROID(xc, c, t);
    v_rot_z = (real)ROT_OMEGA * xc[0];

    dS[eqn] = -(real)ROT_C * sw;
    return -(real)ROT_C * sw * (C_W(c, t) - v_rot_z);
}

/* ========================== END OF FILE ================================== */

import sys
sys.path.insert(0, './builddir/dynmod/apps/libs/pymodule')

import pyBioCMAMCST
import numpy as np
import scipy.optimize 

class SimpleModel:
    tauPTS = 25.0
    tau_f = 5.0
    tau_d = 5.0
    NPermease_max = 5e4
    tauAu = 5.0
    tauAd = 5.0
    tau_metabolisme = 1.0
    phi_pts_max = 4.454e-12
    kpts = 1e-3
    kppermease = 1e-2
    psi_permease = 1.25e-13
    YXS = 0.5
    YXO = 1e-4
    YSA = 1e-4
    psi_o_meta = 20e-3 / 3600 * 32e-3
    critcal_division_mass = 1.7e-5


    class Xi:
        def __init__(self):
            self.mass = 0.0
            self.mu_eff = 0.0
            self.a_pts = 0.0
            self.a_permease = 0.0
            self.n_permease = 0.0

        def __add__(self, rhs):
            result = SimpleModel.Xi()
            result.mass = self.mass + rhs.mass
            result.mu_eff = self.mu_eff + rhs.mu_eff
            result.a_pts = self.a_pts + rhs.a_pts
            result.a_permease = self.a_permease + rhs.a_permease
            result.n_permease = self.n_permease + rhs.n_permease
            return result

        def __mul__(self, scalar):
            result = SimpleModel.Xi()
            result.mass = self.mass * scalar
            result.mu_eff = self.mu_eff * scalar
            result.a_pts = self.a_pts * scalar
            result.a_permease = self.a_permease * scalar
            result.n_permease = self.n_permease * scalar
            return result

    def __init__(self):
        self.phi_s_in = 0.0
        self.phi_o_in = 0.0
        self.mu_p = 0.0
        self.xi = self.Xi()
        self.xi_dot = self.Xi()

    # def step(self, d,xi_dot):
        

    def __str__(self):
        return f'SimpleModel: phi_s_in={self.phi_s_in}, phi_o_in={self.phi_o_in}, mu_p={self.mu_p},mass={self.xi.mass}'

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (isinstance(other, SimpleModel) and
                self.phi_s_in == other.phi_s_in and
                self.phi_o_in == other.phi_o_in and
                self.mu_p == other.mu_p and
                self.xi == other.xi and
                self.xi_dot == other.xi_dot)

def phi_pts(model, S):
    return model.xi.a_pts * (SimpleModel.phi_pts_max * S / (SimpleModel.kpts + S))

def phi_permease(model, n_permease, S):
    return n_permease * SimpleModel.psi_permease * model.xi.a_permease * \
           (SimpleModel.phi_pts_max * S / (SimpleModel.kppermease + S))

def uptake_glucose(model, n_permease, S):
    if S == 0:
        model.phi_s_in = 0.0
        return 0.0

    tau_m = 1e-3

    def get_phi(Si):
        phi_s_pts = phi_pts(model, Si)
        phi_s_permease = phi_permease(model, n_permease, Si)
        rhs = phi_s_pts + phi_s_permease
        lhs = (S - Si) / tau_m
        return abs(rhs - lhs)

    Si = scipy.optimize.newton(get_phi, S, tol=1e-5)
    if Si < 0 :
        return 0.0

    phi_s_pts = phi_pts(model, Si)
    phi_s_in = phi_s_pts + phi_permease(model, n_permease, Si)
    phi_s_in = 0 if np.isnan(phi_s_in) else phi_s_in

    model.phi_s_in = phi_s_in

    gamma_PTS_S = phi_s_pts / SimpleModel.phi_pts_max
    return gamma_PTS_S



def uptake_o2(model, O):
    tau_m = 1e-3
    growth = 0

    def get_phi(Oi):
        phi_o_growth = SimpleModel.YXO * model.xi.mu_eff
        phi_o_in = model.xi.mass*SimpleModel.psi_o_meta + growth * phi_o_growth
        rhs = phi_o_growth + phi_o_in
        lhs = (O - Oi) / tau_m
        return abs(rhs - lhs)

    if O == 0:
        return 0.0

    try:
        Oi = scipy.optimize.newton(get_phi, O, tol=1e-5)
    except RuntimeError:
        return 0.0

    if Oi < 0:
        return 0.0

    growth = 1
    try:
        Oi2 = scipy.optimize.newton(get_phi, O, tol=1e-5)
    except RuntimeError:
        return model.xi.mass*SimpleModel.psi_o_meta

    phi_o_growth = SimpleModel.YXO * model.xi.mu_eff
    return SimpleModel.psi_o_meta + growth * phi_o_growth

def update_xi_dot(model, gamma_PTS_S, n_permease, S):
    model.xi_dot.mass = model.xi.mu_eff

    model.xi_dot.a_pts = (1.0 / SimpleModel.tauPTS) * \
                         ((S / SimpleModel.kpts + S) - model.xi.a_pts)

    model.xi_dot.a_permease = \
        ((1.0 / SimpleModel.tauAu) * gamma_PTS_S +
         (1.0 / SimpleModel.tauAd) * (1.0 - gamma_PTS_S)) * \
        (1.0 - gamma_PTS_S - model.xi.a_permease)

    model.xi_dot.n_permease = \
        (1.0 / SimpleModel.tau_f) * \
        (SimpleModel.kppermease / (SimpleModel.kppermease + S)) + \
        (1.0 / SimpleModel.tau_d) * \
        (S / (SimpleModel.kpts + S)) * \
        (SimpleModel.NPermease_max * (1.0 - gamma_PTS_S) - n_permease)

    model.xi_dot.mu_eff = \
        (1.0 / SimpleModel.tau_metabolisme) * \
        (model.mu_p - model.xi.mu_eff)



def init_kernel(particle):
    # print("init")
    ptr = particle.getOpaque()
    model = SimpleModel()
    model.xi.mass = SimpleModel.critcal_division_mass-1e-6
    model.xi.a_permease = 0.5
    model.xi.a_pts = 0.5
    model.xi.mu_eff = 1e-5
    model.xi.n_permease = 1e3

    model.mu_p = 0.0
    model.phi_o_in = 0.0
    model.phi_s_in = 0.0
    ptr.init(model)
   
   
def update_kernel(d_t, p, concentrations):
    opaque = p.getOpaque()
    model = opaque.cast()
    S = concentrations[0]
    n_permease = int(model.xi.n_permease)

    gamma_PTS_S = uptake_glucose(model, n_permease, S)
    phi_o_in = uptake_o2(model, concentrations[1])

    mu_p = min(SimpleModel.YXS * model.phi_s_in, SimpleModel.YXO * phi_o_in)

    model.mu_p = mu_p
    model.phi_o_in = phi_o_in

    update_xi_dot(model, gamma_PTS_S, n_permease, S)

    model.xi = model.xi + model.xi_dot*d_t
   
    if model.xi.mass >= SimpleModel.critcal_division_mass:
        p.status = pyBioCMAMCST.CellStatus.CYTOKINESIS



def contribution_kernel(p, contribution):
    opaque = p.getOpaque()
    model = opaque.cast()
    ic = int(p.current_container)
    contribution[0, ic] += -model.phi_s_in * p.weight
    contribution[1, ic] += -model.phi_o_in * p.weight
    contribution[2, ic] += SimpleModel.YSA * max(0.0, model.mu_p - model.xi.mu_eff) * p.weight

    

   
def division_kernel(p,child):
    opaque = p.getOpaque()
    model = opaque.cast()
    model.xi.mass = model.xi.mass / 2
    p.status = pyBioCMAMCST.CellStatus.IDLE
    child.status = pyBioCMAMCST.CellStatus.IDLE

    child_opaque = child.getOpaque()

    child_opaque.init(model)
    

    return 
  
def __debug(p):
    opaque = p.getOpaque()
    model = opaque.cast()
    print(model)



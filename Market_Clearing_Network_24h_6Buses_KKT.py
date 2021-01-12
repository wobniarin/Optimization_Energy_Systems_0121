# Optimization Course - 24h Market Clearing with network KKT Conditions - IEEE 6 Bus network
# Author: Íngrid Munné-Collado
# Date: 12/01/2021

# Requirements: Install pyomo, and ipopt version ipopt=3.11.1
from pyomo.environ import *
from pyomo.mpec import *
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

# Creating the model
model = AbstractModel()

# Defining Sets
model.G = Set() # generators
model.D = Set() # demand
model.N = Set() # buses in the network
model.L = Set() # Lines in the network
model.T = Set() # Time periods

# Defining Parameters
model.Pgmax = Param(model.G, model.T)
model.Pdmax = Param(model.D, model.T)
model.costs_g = Param(model.G)
model.costs_d = Param(model.D)
model.Fmaxnn = Param(model.N, model.N, mutable=True)
model.Bnn = Param(model.N, model.N)
model.location_generators = Param(model.G, model.N)
model.location_demands = Param(model.D, model.N)

# Defining Variables
model.pd = Var(model.D, model.T, within=NonNegativeReals)
model.pg = Var(model.G, model.T, within=NonNegativeReals)
model.thetan = Var(model.N, model.T,  within=Reals)
model.flownm = Var(model.L, model.T, within=Reals)

# Defining Dual Variables
model.lambda_nt= Var(model.N, model.T, within=Reals)
model.mu_dt_MIN = Var(model.D, model.T, within=NonNegativeReals)
model.mu_dt_MAX = Var(model.D, model.T, within=NonNegativeReals)
model.mu_gt_MIN = Var(model.G, model.T, within=NonNegativeReals)
model.mu_gt_MAX= Var(model.G, model.T, within=NonNegativeReals)
model.eta_nm_t_up = Var(model.N, model.N, model.T, within=NonNegativeReals)
model.eta_nm_t_down = Var(model.N, model.N, model.T, within=NonNegativeReals)
model.gamma_t = Var(model.N, model.T, within=Reals)

# Defining PRIMAL constraints
# C1 demand  max constraint
def pd_MAX_limit(model,d, t):
    return model.pd[d,t] <= model.Pdmax[d,t]
model.pd_max_limit = Constraint(model.D,model.T, rule=pd_MAX_limit)

# C2 generators max constraint
def pg_MAX_limit(model,g, t):
    return model.pg[g,t] <= model.Pgmax[g,t]
model.pgmax_limit = Constraint(model.G, model.T, rule=pg_MAX_limit)

# # C4 Power flow Upper bound
# def Powerflownm(model, n, m, t):
#     pf_nm = model.Bnn[n,m] * (model.thetan[n,t] - model.thetan[m,t])
#     return pf_nm <= model.Fmaxnn[n,m]
# model.powerflow = Constraint(model.N, model.N, model.T, rule=Powerflownm)
#
# # C5 SLACK BUS
# def slack(model,t):
#     return model.thetan[0,t] == 0
# model.slackbus = Constraint(model.T, rule=slack)
#
# # C6 NODAL POWER BALANCE
# def nodalPowerBalancen(model, n,t):
#     gen_node_n_t = sum(model.pg[g,t] * model.location_generators[g,n] for g in model.G)
#     dem_node_n_t = sum(model.pd[d,t] * model.location_demands[d,n] for d in model.D)
#     powerflow_n_t = sum(model.Bnn[n,m] * (model.thetan[n,t] - model.thetan[m,t]) for m in model.N)
#     return dem_node_n_t + powerflow_n_t - gen_node_n_t == 0
# model.nodalPFB = Constraint(model.N, model.T, rule=nodalPowerBalancen)

# DUAL CONSTRAINTS
# C1 DUAL DEMAND
def dual_c_demand(model,d, t, n):
    return model.mu_dt_MAX[d,t] + model.costs_d[d] - model.mu_dt_MIN[d,t] +  model.location_demands[d,n] * model.lambda_nt[n,t] == 0
model.dualcdemand = Constraint(model.D, model.T, model.N, rule = dual_c_demand)

# C2 DUAL GENERATION
def dual_c_gen(model,g, t,n):
    return model.mu_gt_MAX[g,t] + model.costs_g[g] - model.mu_gt_MIN[g,t] - model.location_generators[g,n]*model.lambda_nt[n,t] == 0
model.dualcgen = Constraint(model.G, model.T, model.N, rule = dual_c_gen)

# C3 DUAL THETA ANGLE
def dual_angle(model, n, m, t ):
    return sum( model.Bnn[n,m] * (model.lambda_nt[n,t] - model.lambda_nt[m,t] + model.eta_nm_t_up[n,m,t] \
                                  - model.eta_nm_t_up[m,n,t] -  model.eta_nm_t_down[n,m,t] \
                                  + model.eta_nm_t_down[m,n,t] ) for m in model.N) == 0
model.dual_angle = Constraint(model.N, model.N, model.T, rule= dual_angle)

# C4 DUAL THETA SLACK BUS
def dual_angle_slack(model, n, m, t ):
    return sum( model.Bnn[n,m] * (model.lambda_nt[n,t] - model.lambda_nt[m,t] + model.eta_nm_t_up[n,m,t] \
                                  - model.eta_nm_t_up[m,n,t] -  model.eta_nm_t_down[n,m,t] \
                                  + model.eta_nm_t_down[m,n,t] ) for m in model.N)  + model.gamma_t[0,t] == 0
model.dual_angle_SLACK = Constraint(model.N, model.N, model.T, rule= dual_angle_slack)

# Complementarity conditions
# CC1
def cc1(model, d, t):
    return complements(pd_MAX_limit(model,d, t), model.mu_dt_MAX[d,t] >= 0)
model.comp_1 = Complementarity(model.D, model.T, rule=cc1)

# CC2
def cc2(model, d, t):
    return complements(0 <= model.pd[d,t], model.mu_dt_MIN[d,t] >= 0)
model.comp_2= Complementarity(model.D, model.T, rule=cc2)

# CC3
def cc3(model, g, t):
    return complements(pg_MAX_limit(model,g, t), model.mu_gt_MAX[g,t] >=0)
model.comp_3 = Complementarity(model.G, model.T, rule=cc3)

# CC4
def cc4(model, g, t):
    return complements(0 <= model.pg[g,t], model.mu_gt_MIN[g,t] >= 0)
model.comp_4 = Complementarity(model.G, model.T, rule=cc4)

# CC5
def cc5(model,n,m,t):
    return complements(0 <= model.Bnn[n,m] * (model.thetan[n,t] - model.thetan[m,t]) + model.Fmaxnn[n,m], model.eta_nm_t_down[n,m,t] >=0)
model.comp_5 = Complementarity(model.N, model.N, model.T, rule=cc5)

# CC6
def cc6(model,n,m,t):
    return complements(0 <= - model.Bnn[n,m] * (model.thetan[n,t] - model.thetan[m,t]) + model.Fmaxnn[n,m], model.eta_nm_t_down[n,m,t] >=0)
model.comp_6 = Complementarity(model.N, model.N, model.T, rule=cc6)

# choose the solver
opt = SolverFactory('mpec_nlp')

# Create a model instance and optimize
instance = model.create_instance('data_Market_Clearing_Network_24_6bus.dat')

# so the solver plugin will know which suffixes to collect
instance.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)

# Solve the optimization problem
results = opt.solve(instance)

# Display results of the code
# instance.display()
results.write()
print("\nDisplaying Solution\n" + '-' * 60)
model.display()

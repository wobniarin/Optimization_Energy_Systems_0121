# Optimization Course - 24h Market Clearing with network DUAL - IEEE 6 Bus network
# Author: Íngrid Munné-Collado
# Date: 11/01/2021

# Requirements: Install pyomo, glpk and gurobi. You should apply for an academic license in gurobi
from pyomo.environ import *
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np

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

# Defining Dual Variables
model.lambda_nt= Var(model.N, model.T, within=Reals)
model.mu_dt = Var(model.D, model.T, within=NonNegativeReals)
model.mu_gt = Var(model.G, model.T, within=NonNegativeReals)
model.eta_nm_t_up = Var(model.N, model.N, model.T, within=NonNegativeReals)
model.eta_nm_t_down = Var(model.N, model.N, model.T, within=NonNegativeReals)
model.gamma_t = Var(model.N, model.T, within=Reals)

# Defining DUAL Objective Function
def SW(model):
    return sum(sum(model.costs_d[d] * model.mu_dt[d,t] for d in model.D) \
               + sum(model.costs_g[g] * model.mu_gt[g,t] for g in model.G) \
               + sum(sum(model.Fmaxnn[n,m] *(model.eta_nm_t_down[n,m,t] + model.eta_nm_t_up[n,m,t]) for m in model.N) for n in model.N) for t in model.T)
model.social_welfare = Objective(rule=SW, sense=minimize)

# Defining constraints
# C1 DUAL DEMAND
def dual_c_demand(model,d, t, n):
    return - model.costs_d[d] + model.mu_dt[d,t] +  model.location_demands[d,n] * model.lambda_nt[n,t] >= 0
model.dualcdemand = Constraint(model.D, model.T, model.N, rule = dual_c_demand)

# C2 DUAL GENERATION
def dual_c_gen(model,g, t,n):
    return model.costs_g[g] + model.mu_gt[g,t] +  model.location_generators[g,n]*model.lambda_nt[n,t] >= 0
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

# choose the solver
opt = pyo.SolverFactory('gurobi')

## in order to solve the problem we have to run this command in a terminal prompt
## pyomo solve --solver=glpk Transport_problem_example_pyomo.py datos.dat

# Create a model instance and optimize
instance = model.create_instance('data_Market_Clearing_Network_24_6bus.dat')

# Create a "dual" suffic component on the instance
# so the solver plugin will know which suffixes to collect
instance.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)

# Solve the optimization problem
results = opt.solve(instance)
# Display results of the code.
instance.display()

# Display all dual variables
print("Duals")
for c in instance.component_objects(pyo.Constraint, active = True):
    print("    Constraint, c")
    for index in c:
        print("      ", index, instance.dual[c[index]])


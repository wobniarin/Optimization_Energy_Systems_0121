# Optimization Course - Lecture 1 Exercise 1 using PYOMO
# Author: Íngrid Munné-Collado
# Date: 09/01/2021

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

# Defining Variables
model.pd = Var(model.D, model.T, within=NonNegativeReals)
model.pg = Var(model.G, model.T, within=NonNegativeReals)
model.thetan = Var(model.N, model.T,  within=Reals)
model.flownm = Var(model.L, model.T, within=Reals)

# Defining Objective Function
def SW(model):
    return sum(sum(model.costs_d[d] * model.pd[d,t] for d in model.D) - sum(model.costs_g[g] * model.pg[g,t] for g in model.G) for t in model.T)
model.social_welfare = Objective(rule=SW, sense=maximize)

# Defining constraints
# C1 demand  max constraint
def pd_MAX_limit(model,d, t):
    return model.pd[d,t] <= model.Pdmax[d,t]
model.pd_max_limit = Constraint(model.D,model.T, rule=pd_MAX_limit)

# C2 generators max constraint
def pg_MAX_limit(model,g, t):
    return model.pg[g,t] <= model.Pgmax[g,t]
model.pgmax_limit = Constraint(model.G, model.T, rule=pg_MAX_limit)

# C4 Power flow Upper bound
def Powerflownm(model, n, m, t):
    pf_nm = model.Bnn[n,m] * (model.thetan[n,t] - model.thetan[m,t])
    return pf_nm <= model.Fmaxnn[n,m]
model.powerflow = Constraint(model.N, model.N, model.T, rule=Powerflownm)

# C5 SLACK BUS
def slack(model,t):
    return model.thetan[0,t] == 0
model.slackbus = Constraint(model.T, rule=slack)

# C6 NODAL POWER BALANCE
def nodalPowerBalancen(model, n,t):
    gen_node_n_t = sum(model.pg[g,t] * model.location_generators[g,n] for g in model.G)
    dem_node_n_t = sum(model.pd[d,t] * model.location_demands[d,n] for d in model.D)
    powerflow_n_t = sum(model.Bnn[n,m] * (model.thetan[n,t] - model.thetan[m,t]) for m in model.N)
    return dem_node_n_t + powerflow_n_t - gen_node_n_t == 0
model.nodalPFB = Constraint(model.N, model.T, rule=nodalPowerBalancen)

# choose the solver
opt = pyo.SolverFactory('gurobi')

## in order to solve the problem we have to run this command in a terminal prompt
## pyomo solve --solver=glpk Transport_problem_example_pyomo.py datos.dat

# Create a model instance and optimize
instance = model.create_instance('data_L1_E1_NETWORK_24_6BUS.dat')

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




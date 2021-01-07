# Optimization Course - Lecture 1 Exercise 1 using PYOMO
# Author: Íngrid Munné-Collado
# Date: 07/01/2021

# Requirements: Install pyomo, glpk and gurobi. You should apply for an academic license in gurobi
from pyomo.environ import *
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np

#Import data
file = 'data_L1_E1.xlsx'
xl = pd.ExcelFile(file)

# Read each sheet
generators = xl.parse('Generator')
demand = xl.parse('Demand')
capacity = xl.parse('Capacity_lines')
susceptance = xl.parse('Susceptance_lines')

# Creating the model
model=AbstractModel()

# Defining Sets
model.G = Set() #generators
model.D = Set() #demand
model.N = Set() #buses in the network
model.L = Set(initialize=model.N*model.N) # lines in the network

# Parameters
init_pg_max = generators['Pmax'].to_dict()
init_pd_max = demand['Pmax'].to_dict()
init_cost_gen = generators['Offer'].to_dict()
init_cost_dem = demand['Bid'].to_dict()
# capacity =
# susceptance =

# Defining Parameters
model.Pgmax = Param(model.G, initialize=init_pg_max, within=NonNegativeReals)
model.Pdmax = Param(model.D, initialize=init_pd_max, within=NonNegativeReals)
model.costs_g = Param(model.G, initialize=init_cost_gen, within=NonNegativeReals)
model.costs_d = Param(model.G, initialize=init_cost_dem, within=NonNegativeReals)
# model.Fmax = Param(model.N, model.N, initialize={line:capacity[line] for line in model.L}, within=NonNegativeReals)
# model.B = Param(model.N, model.N, initialize={line:susceptance[line] for line in model.L}, within=NonNegativeReals)


# Defining Variables
model.pd = Var(model.D,initialize=50, within=NonNegativeReals)
model.pg = Var(model.G,initialize=50, within=NonNegativeReals)
model.theta = Var(model.N,initialize=0, within=Reals)
model.fnm = Var(model.N, initialize=0, within=Reals)

# Defining Objective Function
def SW(model):
    return sum(model.costs_d[d] * model.pd[d] for d in model.D) - sum(model.costs_g[g] * model.pg[g] for g in model.G)
model.social_welfare = Objective(rule=SW, sense=maximize)

# Defining constraints
# C1 demand  max constraint
def pd_MAX_limit(model,d):
    return model.pd[d] <= model.Pdmax[d]
model.pd_max_limit = Constraint(model.D, rule=pd_MAX_limit)

# C2 generatos max constraint
def pg_MAX_limit(model,g):
    return model.pg[g] <= model.Pgmax[g]
model.pgmax_limit = Constraint(model.G, rule=pg_MAX_limit)

# # C3 PowerFlow Definition
# def fnm(model, n,m):
#     return model.fnm = B[n][m] * (model.theta[n] - model.theta[m])
# model.fnm_flow = Constraint(N, rule=fnm)

# # C4 Powerflow Limits
# def
#
# # C5 SLACK
# def theta_ref(model):
#     return theta[1] = 0
# model.theta_slack = Constraint(N, rule= theta_ref)

# C6 PowerFlow Balance at node

# choose the solver
opt = pyo.SolverFactory('gurobi')

# create the model instance, solve and display results.
instance = model.create_instance()
results = opt.solve(instance)
instance.display()

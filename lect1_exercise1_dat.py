# Optimization Course - Lecture 1 Exercise 1 using PYOMO
# Author: Íngrid Munné-Collado
# Date: 07/01/2021

# Requirements: Install pyomo, glpk and gurobi. You should apply for an academic license in gurobi
from pyomo.environ import *
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np

# Creating the model
model = AbstractModel()

# Defining Sets
model.G = Set() #generators
model.D = Set() #demand
model.N = Set() #buses in the network


# Defining Parameters
model.Pgmax = Param(model.G)
model.Pdmax = Param(model.D)
model.costs_g = Param(model.G)
model.costs_d = Param(model.D)

# Defining Variables
model.pd = Var(model.D, within=NonNegativeReals)
model.pg = Var(model.G, within=NonNegativeReals)


# Defining Objective Function
def SW(model):
    return sum(model.costs_d[d] * model.pd[d] for d in model.D) - sum(model.costs_g[g] * model.pg[g] for g in model.G)
model.social_welfare = Objective(rule=SW, sense=maximize)

# Defining constraints
# C1 demand  max constraint
def pd_MAX_limit(model,d):
    return model.pd[d] <= model.Pdmax[d]
model.pd_max_limit = Constraint(model.D, rule=pd_MAX_limit)

# C2 generators max constraint
def pg_MAX_limit(model,g):
    return model.pg[g] <= model.Pgmax[g]
model.pgmax_limit = Constraint(model.G, rule=pg_MAX_limit)

def power_balance(model, g, d):
    return sum(model.pd[d] for d in model.D) - sum(model.pg[g] for g in model.G) == 0
model.powerbalance = Constraint(model.G, model.D, rule=power_balance)

# choose the solver
opt = pyo.SolverFactory('gurobi')

## in order to solve the problem we have to run this command in a terminal prompt
## pyomo solve --solver=glpk Transport_problem_example_pyomo.py datos.dat

# Create a model instance and optimize
instance = model.create_instance('data_L1_E1.dat')
results = opt.solve(instance)
instance.display()

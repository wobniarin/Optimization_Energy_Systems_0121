# Optimization Course - Lecture 1 Exercise 1 using PYOMO
# Author: Íngrid Munné-Collado
# Date: 06/01/2021

# Requirements: Install pyomo, glpk and gurobi. You should apply for an academic license in gurobi
from pyomo.environ import *
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np

# Importing excel sheet data into data frames
Unit = pd.read_excel(r'C:\Users\ingri\PycharmProjects\Optimization_Energy_Systems_0121\Lect1Ex1.xlsx',sheet_name='Unit')
Bus = pd.read_excel(r'C:\Users\ingri\PycharmProjects\Optimization_Energy_Systems_0121\Lect1Ex1.xlsx',sheet_name='Bus')
Branch = pd.read_excel(r'C:\Users\ingri\PycharmProjects\Optimization_Energy_Systems_0121\Lect1Ex1.xlsx',sheet_name='Branch')


#split and list the essential variables that will be used to construct the model
Offer = list(Unit["Offer"])
pmax = list(Unit["Pmax"])
pmin = list(Unit["Pmin"])
Pdmax = list(Bus["Pmax"])
Bid = list(Bus["Bid"])
X = list(Branch["X"])
Lim = list(Branch["Lim"])

# This portion overcomplicates the code and the nodal balance constraint, can be made simpler by using nxm matrix for branch parameters.
Branch2 = Branch[['To','From','R','X','Lim']]
Branch2=Branch2.rename(columns={"To": "From", "From": "To"})
FullBranch= pd.concat([Branch,Branch2], axis=0, ignore_index=True,sort=0)
From = list(FullBranch["From"])
To =list(FullBranch["To"])

#Define Sets of Indexes, this will not work on a problem where multiple generators are connected to the same bus.
N = M = list(Bus.index.map(int))
G = list(Unit.index.map(int))
D = list(Bus.index.map(int))


# Start the model
model = ConcreteModel()
model.pd = Var(D,initialize=50, within=Reals)
model.pg = Var(G,initialize=50, within=Reals)
model.d = Var(N,initialize=0, within=Reals)

# Objective function
def obj_rule(mod):
    return sum(Bid[d]*mod.pd[d] for d in D)-sum(Offer[g]*mod.pg[g] for g in G)
model.obj= Objective(rule=obj_rule,sense=maximize)

# Constraints
def PDmax_rule(mod,d):
    return mod.pd[d]<=Pdmax[d]
model.PDmax = Constraint(D,rule=PDmax_rule)

#This constraint and all other non-negativity constraints can be removed by adding a zero to inf bound to the
# variables (e.g. model.pd=Var(D,bound(0,inf)))

def PDmin_rule(mod,d):
    return mod.pd[d]>=0
model.PDmin = Constraint(D,rule=PDmin_rule)

def PGmax_rule(mod,g):
    return mod.pg[g]<=pmax[g]
model.PGmax = Constraint(G,rule=PGmax_rule)

def PGmin_rule(mod,g):
    return mod.pg[g]>=pmin[g]
model.PGmin = Constraint(G,rule=PGmin_rule)

#There must be a better way to do this one
def NodBal_rule(mod,n):
    return mod.pd[n]-mod.pg[n] == sum((mod.d[i-1]-mod.d[k-1])/FullBranch.loc[x,"X"] for x,(i,k) in
                                      enumerate(zip(FullBranch["From"], FullBranch["To"])) if i == n+1)
model.NodBal = Constraint(N, rule=NodBal_rule)

def Cap_rule(mod,g):
    return Lim[g] >= (mod.d[From[g]-1]-mod.d[To[g]-1])/X[g]
model.Cap = Constraint(G, rule=Cap_rule)

def Cap2_rule(mod,g):
    return Lim[g] >= (mod.d[To[g]-1]-mod.d[From[g]-1])/X[g]
model.Cap2 = Constraint(G, rule=Cap2_rule)

#choose the solver
opt = pyo.SolverFactory('gurobi')

#create the model instance, solve and display results.

instance = model.create_instance()
results = opt.solve(instance)
instance.display()

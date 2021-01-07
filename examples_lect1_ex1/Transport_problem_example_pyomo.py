# Optimization Course - Lecture 1 Exercise 1 using PYOMO
# Author: Íngrid Munné-Collado
# Date: 07/01/2021
# code based on https://www.youtube.com/watch?v=hTKa7Qc38b4&t=184s


# Requirements: Install pyomo, glpk and gurobi. You should apply for an academic license in gurobi
from pyomo.environ import *
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import pandas as pd
import numpy as np

# Read import data

# Creating the model
model = AbstractModel()
""" This is useful because everything we define as a model will be the equations. 
All data used to solve the model will come in a different file"""

# Defining Sets
model.plantas = Set()
model.mercados = Set()

# Defining Parameters for each set
model.produccion = Param(model.plantas)
model.demanda = Param(model.mercados)
model.costes = Param(model.plantas, model.mercados)

# Defining variables
model.unidades = Var(model.plantas, model.mercados, within = NonNegativeReals)

# Defining ObjectiveFunction in an abstract way
def CosteTotal(model):
    return sum(model.costes[n,i] * model.unidades[n,i]
               for n in model.plantas
               for i in model.mercados)

model.costefinal = Objective(rule=CosteTotal)
# here we are saying that we will save the value of the OF defined previously in
# this variable of the model.

# Defining constraints
def MinDemanda(model, mercado): #it depends on model and each of the markets
    return sum(model.unidades[i,mercado] for i in model.plantas) >= model.demanda[mercado]

model.DemandaConstraint = Constraint(model.mercados, rule=MinDemanda)

def MaxProduccion(model, planta):
        return sum(model.unidades[planta, j] for j in model.mercados) <= model.produccion[planta]
model.ProdConstraint = Constraint(model.plantas, rule=MaxProduccion)

## in order to solve the problem we have to run this command in a terminal prompt
## pyomo solve --solver=glpk Transport_problem_example_pyomo.py datos.dat

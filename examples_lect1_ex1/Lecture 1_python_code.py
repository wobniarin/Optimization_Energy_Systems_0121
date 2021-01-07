# -*- coding: utf-8 -*-
"""
Created on Mon Jan 08 20:53:03 2018

@author: trva
"""

# REMEMBER ALL LISTS / VECTORS START FROM 0 - NEVER FORGET !!!

import numpy as np
import gurobipy as gb
import matplotlib.pyplot as plt
import pickle
from itertools import izip
from collections import defaultdict
import time
import pprint
import scipy.stats as sp
import scipy.integrate as spint
import json, ast
from tempfile import TemporaryFile
from datetime import datetime

print ("--- --- --- ---")

#Conventional Model

def Prim_prob(data):
    
    conv_gen = data['generators']		#index for conventional generators    
    wind_gen = data['wind_gen']   
    nodes = data['nodes']                 #index for nodes     
    demand = data['demand']              # index for demand  
    S = np.size(data['w1']['p_hat']) 		#index for scenarios,  bids must be an array      
   
    # Model alias
    model_primal = gb.Model('primal')

    # Scheduled production of generator gen in DA
    P_G = {}
    for g in conv_gen:
        P_G[g] = model_primal.addVar(lb=0.0, ub=data[g]['capacity'], name='Pg {}'.format(g))
    
    # Wind, dispatched and bid per node
    P_W= {}
    for w in wind_gen:
        P_W[w] = model_primal.addVar(lb=0.0, ub=data[w]['capacity'],vtype=gb.GRB.CONTINUOUS, name='P_W {}'.format(w))

    # Initial angles of the buses
    delta_DA = {}
    for n in nodes: # Bus angles
        delta_DA[n] = model_primal.addVar(lb=-gb.GRB.INFINITY, ub=gb.GRB.INFINITY, vtype=gb.GRB.CONTINUOUS, name='delta_DA {}:'.format(n))
             

    # Real time adjustment for conventional generations
    r = {}
    for g in conv_gen:
        for s in xrange(S):
            r[g,s] = model_primal.addVar(lb=-data[g]['Rmax'], ub=data[g]['Rmax'], vtype=gb.GRB.CONTINUOUS, name='r {0} {1}'.format(g,s))

    # Wind spillage per scenario
    P_W_sp = {}
    for k in wind_gen:
        for s in xrange(S):
            P_W_sp[k,s] = model_primal.addVar(lb=0.0, ub=data[k]['W_hat'][s], name='W spill {0} {1}'.format(k,s))
        
    # Load shed per scenario # GROUP VARIABLES
    l_sh = {}
    for d in demand:
        for s in xrange(S):
            l_sh[d,s] = model_primal.addVar(lb=0.0, ub=data[d]['demand'], name='load shed {0} {1}'.format(d,s))
    
    # Angles
    delta_RT={}
    for n in nodes:
        for s in xrange(S):
            delta_RT[n,s] = model_primal.addVar(lb=-gb.GRB.INFINITY, ub=gb.GRB.INFINITY,vtype=gb.GRB.CONTINUOUS, name='delta_RT {0} {1}:'.format(n,s))
    
    #Flow DA
    f_DA={}
    for n in nodes:
        for m in data[n]['toNode']:
            f_DA[n,m] = model_primal.addVar(lb=-gb.GRB.INFINITY, ub=gb.GRB.INFINITY,vtype=gb.GRB.CONTINUOUS, name='f_DA {0} {1}:'.format(n,m)) 

    #Flow in RT
    f_RT={}
    for n in nodes:
        for m in data[n]['toNode']:
            for s in xrange(S):
                f_RT[n,m,s] = model_primal.addVar(lb=-gb.GRB.INFINITY, ub=gb.GRB.INFINITY,vtype=gb.GRB.CONTINUOUS, name='f_RT {0} {1} {2}:'.format(n,m,s)) 

    
    model_primal.update()  # update of the model with the new data

    # Set the objective in the DA market and balancing market
    model_primal.setObjective(
            gb.quicksum(data[g]['cost']*P_G[g] for g in conv_gen) # DA cost
            +gb.quicksum(data[k]['cost']*P_W[k] for k in wind_gen) # wind cost is zero
            + gb.quicksum( data['w1']['p_hat'][s]*(gb.quicksum(data[g]['cost']*r[g,s] for g in conv_gen)
                    + gb.quicksum(data['VOLL']*l_sh[d,s] for d in demand)) for s in xrange(S)),gb.GRB.MINIMIZE)
    

            
    model_primal.update()  # update of the model with the objective
    
    ##### DA_Constraints ######
    
    #DA power flow definition
    DA_flow={}
    for n in nodes:
        for m in data[n]['toNode']:
            DA_flow[n,m]=model_primal.addConstr(f_DA[n,m],gb.GRB.EQUAL,data['B'][n][m]*(delta_DA[n] - delta_DA[m]))
    
    # Power balance constraints
    pb_DA = {}
    for n in nodes:
        pb_DA[n] = model_primal.addConstr(gb.quicksum(P_G[g] for g in data[n]['conv_gen']) + gb.quicksum(P_W[k] for k in data[n]['wind_gen']) - data[n]['load'] - gb.quicksum(f_DA[n,m] for m in data[n]['toNode']),gb.GRB.EQUAL,0, "pb_DA")

    # Line capacity DA
    f_DA_max={}
    for n in nodes:
	    for m in data[n]['toNode']:
	        f_DA_max[n,m]=model_primal.addConstr(f_DA[n,m],gb.GRB.LESS_EQUAL,data['lineCapacity'][n][m])
#	        m.addConstr(f_DA[n,m],gb.GRB.GREATER_EQUAL,-data['lineCapacity'][n][m])        
            
    #RT constraints

    #RT power flow definition
    RT_Flow={}
    for n in nodes:
        for m in data[n]['toNode']:
            for s in xrange(S):
                RT_Flow[n,m,s] = model_primal.addConstr(f_RT[n,m,s],gb.GRB.EQUAL,data['B'][n][m]*(delta_RT[n,s] - delta_RT[m,s]))
    
    #Operation constraints    
    pb_RT = {}
    for n in nodes:
        for s in xrange(S):
            pb_RT[n,s] = model_primal.addConstr( gb.quicksum(r[g,s] for g in data[n]['conv_gen']) # RT adj
            + gb.quicksum(l_sh[d,s] for d in data[n]['demand'])  # load shedding
            + gb.quicksum(data[k]['W_hat'][s] - P_W[k] - P_W_sp[k,s] for k in data[n]['wind_gen']) # wind scenario
            -gb.quicksum(f_RT[n,m,s] - f_DA[n,m] for m in data[n]['toNode']),gb.GRB.EQUAL,0.0 #power flow
            ,name="pb_RT in node{0} of scenario{1}".format(n,s)
            )
        
    # Power flow limits
    f_RT_max={}
    for n in nodes:
        for m in data[n]['toNode']:
            for s in xrange(S):
               f_RT_max[n,m,s]=model_primal.addConstr(f_RT[n,m,s],gb.GRB.LESS_EQUAL,data['lineCapacity'][n][m])                           
#               f_RT_max[n,m,k]=m.addConstr(f_RT[n,m,k],gb.GRB.GREATER_EQUAL,-data['lineCapacity'][n][m])    
    
    # Constraint for the initialisation of angles
    model_primal.addConstr(delta_DA['n1'],gb.GRB.EQUAL,0)

    for s in xrange(S):
        model_primal.addConstr(delta_RT['n1',s],gb.GRB.EQUAL,0)    
        
    #RT plus DA generation limit
    upgenlimits = {}
    downgenlimits = {}
    for g in conv_gen:
        for s in xrange(S):
            upgenlimits[g,s] = model_primal.addConstr(r[g,s] + P_G[g] ,gb.GRB.LESS_EQUAL,data[g]['capacity'] )
            downgenlimits[g,s] = model_primal.addConstr(r[g,s] + P_G[g] ,gb.GRB.GREATER_EQUAL,0 )

    model_primal.update()         
    
    model_primal.optimize()
    
    # Cost recovery
    
#    # 1. Per scenario DA +RT
    profit_g={}
    for g in conv_gen:
        for s in xrange(S):    
            if data[g]['bus']==1: 
                profit_g[g,s] = P_G[g].x*pb_DA['n1'].pi  + r[g,s].x*((pb_RT['n1',s]).pi/data['w1']['p_hat'][s]) - data[g]['cost']*(P_G[g].x+r[g,s].x)
            else:
                profit_g[g,s] = P_G[g].x*pb_DA['n2'].pi  + r[g,s].x*((pb_RT['n2',s]).pi/data['w1']['p_hat'][s]) - data[g]['cost']*(P_G[g].x+r[g,s].x)
    
    
    profit_w={}
    for k in wind_gen:
        for s in xrange(S):
            if data[k]['bus']==1:
                profit_w[k,s] = P_W[k].x*pb_DA['n1'].pi + (data[k]['W_hat'][s] - P_W[k].x - P_W_sp[k,s].x)*((pb_RT['n1',s]).pi/data['w1']['p_hat'][s])
            else:           
                profit_w[k,s] = P_W[k].x*pb_DA['n2'].pi + (data[k]['W_hat'][s] - P_W[k].x - P_W_sp[k,s].x)*((pb_RT['n2',s]).pi/data['w1']['p_hat'][s])

    payment_load={}
    for d in demand:
        for s in xrange(S):
            payment_load[d,s] = data[d]['demand']*pb_DA['n1'].pi + data['VOLL']*l_sh[d,s].x
    
    #Expected DA + RT
    exp_prof_g = {}
    for g in conv_gen:
        exp_prof_g[g] = sum(data['w1']['p_hat'][s]*profit_g[g,s] for s in xrange(S) )

    exp_prof_w = {}
    for k in wind_gen:
        exp_prof_w[k] = sum(data['w1']['p_hat'][s]*profit_w[k,s] for s in xrange(S) )  
        
    # Revenue Adequacy
    Rev_ad={}
    for s in xrange(S):
        Rev_ad[s] =( 
                P_G['g1'].x*pb_DA['n1'].pi + P_G['g2'].x*pb_DA['n1'].pi + P_G['g3'].x*pb_DA['n2'].pi + P_W['w2'].x*pb_DA['n2'].pi
                + r['g1',s].x*((pb_RT['n1',s]).pi/data['w1']['p_hat'][s]) + r['g2',s].x*((pb_RT['n1',s]).pi/data['w1']['p_hat'][s]) + r['g3',s].x*((pb_RT['n2',s]).pi/data['w1']['p_hat'][s])
                + (data['w2']['W_hat'][s] - P_W['w2'].x - P_W_sp['w2',s].x)*((pb_RT['n2',s]).pi/data['w2']['p_hat'][s])
                - data['d1']['demand']*pb_DA['n1'].pi + l_sh['d1',s].x*((pb_RT['n1',s]).pi/data['w2']['p_hat'][s])
                ) 
 
    # Expected
    Rev_ad_exp ={}
    Rev_ad_exp= (
            P_G['g1'].x*pb_DA['n1'].pi + P_G['g2'].x*pb_DA['n1'].pi +P_G['g3'].x*pb_DA['n2'].pi + P_W['w2'].x*pb_DA['n2'].pi - data['d1']['demand']*pb_DA['n1'].pi # DA
            + sum( data['w1']['p_hat'][s]*r['g1',s].x*((pb_RT['n1',s]).pi/data['w1']['p_hat'][s]) for s in xrange(S) )
            + sum( data['w1']['p_hat'][s]*r['g2',s].x*((pb_RT['n1',s]).pi/data['w1']['p_hat'][s]) for s in xrange(S) )
            + sum( data['w1']['p_hat'][s]*r['g3',s].x*((pb_RT['n2',s]).pi/data['w1']['p_hat'][s]) for s in xrange(S) )
            + data['w2']['p_hat'][s]*l_sh['d1',s].x*((pb_RT['n1',s]).pi/data['w2']['p_hat'][s])
            + sum( data['w1']['p_hat'][s]*(data['w2']['W_hat'][s] - P_W['w2'].x - P_W_sp['w2',s].x)*((pb_RT['n2',s]).pi/data['w1']['p_hat'][s]) for s in xrange(S) )
            )
                
    
    solution_primal = {
        'Op. Cost' : model_primal.ObjVal,
        'Wind generators DA': [P_W[k].x for k in wind_gen],
        'Conventional generators DA': [P_G[g].x for g in conv_gen],
        'DA_Price' : [pb_DA[n].pi for n in nodes],
        'Balancing Prices' : [[pb_RT[n,s].pi for n in nodes]/data['w1']['p_hat'][s] for s in xrange(S)],
        'Balancing Duals' : [[pb_RT[n,s].pi for n in nodes] for s in xrange(S)],
        'angles_DA': [ [delta_DA[n].x] for n in nodes],
        'angles': [ [[delta_RT[n,s].x] for n in nodes] for s in xrange(S)],
        'load shed' : [ [l_sh[d,s].x for d in demand] for s in xrange(S)],
        'wind spilled' : [ [P_W_sp[k,s].x for k in wind_gen] for s in xrange(S)],
        'Adjustment volumes' : [[r[g,s].x for g in conv_gen] for s in xrange(S) ],
        'f_DA': [f_DA[n,m].x for n in nodes for m in data[n]['toNode'] ],
        'f_RT': [[f_RT[n,m,s].x for n in nodes for m in data[n]['toNode']] for s in xrange(S) ],
        'RT_FLow': [[RT_Flow[n,m,s].pi for n in nodes for m in data[n]['toNode']] for s in xrange(S) ], 
        'FLow_max': [[f_RT_max[n,m,s].pi for n in nodes for m in data[n]['toNode']] for s in xrange(S) ],
        'Profits of generators' :  [[profit_g[g,s] for g in conv_gen] for s in xrange(S) ],
        'Profits of wind generators' :  [[profit_w[k,s] for k in wind_gen] for s in xrange(S) ],
        'Charge to consumers' : [[payment_load[d,s] for d in demand] for s in xrange(S) ],
        'Expected profits RT + DA of conv. gen.' : [exp_prof_g[g] for g in conv_gen],
        'Expected profits RT + DA of wind. gen.' : [exp_prof_w[k] for k in wind_gen],
        'Revenue adequacy per scenario.' : [Rev_ad[s] for s in xrange(S)],
        'Expected Revenue adequacy.' : [Rev_ad_exp],
        
    }
   
    model_primal.write('primal.lp') 
    return solution_primal

def system_input():
    return {
    'generators' : ['g1', 'g2','g3'],
    'wind_gen': ['w1', 'w2'],
    'demand' : ['d1','d2'],
    'nodes' : ['n1','n2'],
    
    'g1': {'bus':1, 'capacity': 50, 'cost': 10, 'Rmax' : 0},
    'g2': {'bus':1, 'capacity': 110,'cost': 25, 'Rmax' : 0},
    'g3': {'bus':2, 'capacity': 100,'cost': 35, 'Rmax' : 45},
    
    'w1': {'bus':1,'capacity': 0, 'min':0, 'cost': 0, 'W_hat' : bids1, 'p_hat' : probability},
    'w2': {'bus':2,'capacity': 50, 'min':0, 'cost': 0, 'W_hat' : bids2, 'p_hat' : probability},
    
    'd1' : {'demand' : 200},
    'd2' : {'demand' : 0},
    
    'n1' : {'conv_gen' : ['g1','g2'], 'wind_gen' : ['w1'], 'demand':['d1'], 'load' : 200, 'toNode' : ['n2'], 'wind_min':0, 'wind_max':0 },
    'n2' : {'conv_gen' : ['g3'], 'wind_gen' : ['w2'], 'demand':['d2'], 'load' : 0, 'toNode' : ['n1'], 'wind_min':0, 'wind_max':70},	

    'lineCapacity':{'n1': {'n2':200}, 'n2':{'n1':200} },
    'B': {'n1':{'n2':1000}, 'n2': {'n1':1000} },
    
    'BigM' : 1000,
    'round' : 5,
    'VOLL' : 200
    }


# SCENARIOS
probability = np.array([0.2, 0.5, 0.3])
bids1 = np.array([0, 0, 0])
bids2 = np.array([50, 22, 10])

#t=wind_input(bid,probability)

pp1 = pprint.PrettyPrinter(indent = 2, width = 50, depth = None)
pp2 = pprint.PrettyPrinter(indent = 2, width = 250, depth = None)

# Generating results


primal_s= Prim_prob(system_input())
pp2.pprint(primal_s)






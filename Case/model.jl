using JuMP
using Gurobi
using DelimitedFiles
using CSV
using DataFrames

# import, define and pretreat data
node_matrix = CSV.read("data//node_matrix.csv",DataFrame)
prod_matrix = CSV.read("data//product_matrix.csv",DataFrame)
demand_matrix = CSV.read("data//demand_matrix.csv",DataFrame)
supply_matrix = CSV.read("data//supply_matrix.csv",DataFrame)
alpha_matrix = readdlm("data//alpha_matrix.csv",',');
site_matrix = CSV.read("data//site_matrix.csv",DataFrame)
time_matrix = CSV.read("data//time_matrix.csv",DataFrame)
flow_matrix = CSV.read("data//flow_matrix.csv",DataFrame)
storage_matrix = CSV.read("data//storage_matrix.csv",DataFrame)

# nodes, products, customers, suppliers and technologies
NODES = node_matrix[:,1]; # all nodes
PRODS = prod_matrix[:,1] # all products
DEMS  = demand_matrix[:,1] # all demands
SUPS  = supply_matrix[:,1] # all supply
TECH_PRVD = site_matrix[:,1] # all technology providers
TIME  = time_matrix[:,1] # all period
FLOW  = flow_matrix[:,1] # all flows
STOR  = storage_matrix[:,1] # all storage

# node properties
node_alia = Dict(zip(NODES, node_matrix[:,2])); # node alias

# product properties
prod_name = Dict(zip(PRODS, prod_matrix[:,2])); # product names

# customer properties
dem_node  = Dict(zip(DEMS, demand_matrix[:,2])); # demand locations
dem_prod  = Dict(zip(DEMS, demand_matrix[:,3])); # demand products
dem_bid   = Dict(zip(DEMS, demand_matrix[:,4])); # demand bids
dem_cap   = Dict(zip(DEMS, demand_matrix[:,5])); # demand capacities
dem_time  = Dict(zip(DEMS, demand_matrix[:,6])); # demand time

# supplier properties
sup_node  = Dict(zip(SUPS, supply_matrix[:,2])); # supply locations
sup_prod  = Dict(zip(SUPS, supply_matrix[:,3])); # supply products
sup_bid   = Dict(zip(SUPS, supply_matrix[:,4])); # supply bids
sup_cap   = Dict(zip(SUPS, supply_matrix[:,5])); # supply capacities
sup_time  = Dict(zip(SUPS, supply_matrix[:,6])); # supply period

# technology providers properties
tp_node = Dict(zip(TECH_PRVD, site_matrix[:,2])); # node location of the technology provider
tp_cap  = Dict(zip(TECH_PRVD, site_matrix[:,3])); # technology provider capacities
tp_bid  = Dict(zip(TECH_PRVD, site_matrix[:,4])); # technology provider bid
tp_time = Dict(zip(TECH_PRVD, site_matrix[:,5])); # technology provider time

# flow providers properties
flow_send = Dict(zip(FLOW, flow_matrix[:,2])) # flow sending node
flow_recv = Dict(zip(FLOW, flow_matrix[:,3])) # flow receving node
flow_prod = Dict(zip(FLOW, flow_matrix[:,4])) # product of flow
flow_cap  = Dict(zip(FLOW, flow_matrix[:,5])) # flow capacities
flow_bid  = Dict(zip(FLOW, flow_matrix[:,6])) # flow bid
flow_send_time = Dict(zip(FLOW, flow_matrix[:,7]))  # flow sending time
flow_recv_time = Dict(zip(FLOW, flow_matrix[:,8]))  # flow receving time

# storage properties
stor_node  = Dict(zip(STOR, storage_matrix[:,2])); # storage locations
stor_prod  = Dict(zip(STOR, storage_matrix[:,3])); # storage products
stor_bid   = Dict(zip(STOR, storage_matrix[:,4])); # storage bids
stor_cap   = Dict(zip(STOR, storage_matrix[:,5])); # storage capacities
stor_time  = Dict(zip(STOR, storage_matrix[:,6])); # storage period

# time subsets for balance constraint definition
Z = TIME
Z1 = [Z[1]] # the first time index needs its own balance; create as a 1D array
ZZ = Z[2:end] # the "interior" time points share a common structure

# transformation factors
transfer = Dict((TECH_PRVD[1],PRODS[1]) => 0.5);
for i in 1:length(TECH_PRVD)
    for k in 1: length(PRODS)
        transfer[(i,k)] = alpha_matrix[i,k];
    end
end

######################################################################
m = Model(Gurobi.Optimizer)

# flow
@variable(m, flow[FLOW] >= 0);
@variable(m, fin[j in NODES, i in NODES, pr in PRODS, z in TIME]>= 0);
@variable(m, fout[i in NODES, j in NODES, pr in PRODS, z in TIME]>= 0);

# demand and supply
@variable(m, dem[DEMS] >= 0);
@variable(m, d[NODES,PRODS,TIME] >= 0);
@variable(m, sup[SUPS] >= 0);
@variable(m, s[NODES,PRODS,TIME] >= 0);

# tech provider
@variable(m, x[TECH_PRVD] >= 0);

# storage
@variable(m, vstor[STOR] >= 0);
@variable(m, v[NODES, PRODS, TIME]>= 0);

# cost
@variable(m,opcost);
@variable(m,storcost);
@variable(m,flowcost);
@variable(m,demrevn);
@variable(m,profit);

# demand, supply, storage and flow
@constraint(m, demeq[n in NODES, pr in PRODS, z in TIME], d[n,pr,z] == sum(dem[dd] for dd in DEMS if dem_prod[dd]==pr
                && dem_node[dd]==n && dem_time[dd]==z));

@constraint(m, supeq[n in NODES, pr in PRODS, z in TIME], s[n,pr,z] == sum(sup[ss] for ss in SUPS if sup_prod[ss]==pr
                && sup_node[ss]==n && sup_time[ss]==z));

@constraint(m, storeq[n in NODES, pr in PRODS, z in TIME], v[n,pr,z] == sum(vstor[vv] for vv in STOR if stor_prod[vv]==pr
                && stor_node[vv]==n && stor_time[vv]==z));

@constraint(m, fineq[j in NODES, i in NODES, pr in PRODS, z in TIME], fin[j,i,pr,z] == sum(flow[ff] for ff in FLOW if flow_prod[ff]==pr
                && flow_recv[ff]==i && flow_send[ff]==j && flow_recv_time[ff]==z));

@constraint(m, fouteq[i in NODES, j in NODES, pr in PRODS, z in TIME], fout[i,j,pr,z] == sum(flow[FF] for FF in FLOW if flow_prod[FF]==pr
                && flow_send[FF]==i && flow_recv[FF]==j && flow_send_time[FF]==z));
                

# balance
@constraint(m, Balance1[i in NODES, pr in PRODS, z in Z1], s[i,pr,z] + sum(fin[j,i,pr,z] for j in NODES) 
                 + sum(x[t]*transfer[t,pr] for t in TECH_PRVD if tp_node[t]==i && tp_time[t]==z)
                == d[i,pr,z] + v[i,pr,z] + sum(fout[i,j,pr,z] for j in NODES));

@constraint(m, BalanceZ[i in NODES, pr in PRODS, z in ZZ], s[i,pr,z] + v[i,pr,z-1] + sum(fin[j,i,pr,z] for j in NODES)
                 + sum(x[t]*transfer[t,pr] for t in TECH_PRVD if tp_node[t]==i && tp_time[t]==z)
                == d[i,pr,z] + v[i,pr,z] + sum(fout[i,j,pr,z] for j in NODES));


                                    
# demand capacity constraints
@constraint(m, demand_capacity[i in DEMS],  dem[i] <= dem_cap[i]);

# # demand fixed constraints
# @constraint(m, demand_fixed[i in NODES, pr in PRODS, z in TIME],  d[5,2,z] == 2500);

# supply capacity constraints
@constraint(m, supply_capacity[i in SUPS],  sup[i] <= sup_cap[i]);

# flow capacity constraints
@constraint(m, flow_capacity[i in FLOW], flow[i] <= flow_cap[i]);

# processing capacity constraints
@constraint(m, x_capacity[t in TECH_PRVD], x[t] <= tp_cap[t] );

# storage capacity constraints
@constraint(m, storage_capacity[i in STOR],  vstor[i] <= stor_cap[i]);

##Objective
@constraint(m, opcost == sum(x[t]*tp_bid[t] for t in TECH_PRVD));
@constraint(m, storcost == sum(vstor[i]*stor_bid[i] for i in STOR));
@constraint(m, flowcost == sum(flow[i]*flow_bid[i] for i in FLOW));
@constraint(m, demrevn == sum(dem[i]*dem_bid[i] for i in DEMS));
@constraint(m, profit == - opcost - storcost - flowcost + demrevn );
@objective(m, Max, profit)
optimize!(m)

Record_Num_Continuous_Var=num_variables(m)
Record_Num_Integer_Var=num_constraints(m, VariableRef, MOI.Integer)
Record_Num_Ineq_con=num_constraints(m,AffExpr, MOI.GreaterThan{Float64})+num_constraints(m,AffExpr, MOI.LessThan{Float64})
Record_Num_Eq_con=num_constraints(m,AffExpr, MOI.EqualTo{Float64})



####INPUT PARAMETERS
param nL;	 			#number of jobs
param nT default 2;  	#number of control options
param nS;	 			#number of scenarios

set L = 1 .. nL; 		#set of jobs -- index l
set T = 1 .. nT; 		#set of control options -- t
set S = 1 .. nS; 		#set of scenarios -- index (i,j)

param p{i in S}:=1/nS; 				#scenario probabilities
param w{l in L};					#job weights
param h{t in T, l in L}; 			#cost coefficient for the control options
param xi{l in L, i in S, t in  T}; 	#processing times

param B default 0; 					#available budget 
param kappa; 						#radius
param alpha; 						#CVAR confidence level


####upper bounds
param M{i in S, j in S, l in L} default 2*max{t in T} abs(xi[l,i,t]-xi[l,j,t]);
param U_tau default 0; 				#Updated in the data file
param kappa_max default 0;

###For the SCHEDULING PROBLEM WITH THE COMONOTONE STRUCTURE
param a{t in T, l in L} default 0; 						#speedup factors a[t,l]\in [0,1]
param xi_base{l in L, i in S} default 0;				#baseline processing times


################################################################
#####Decision Variables

var eta >= 0;
var v{i in S} >= 0;
var tau >= 0, <= U_tau;
var z{t in T, k in L, l in L} >=0, <=1;
var u{t in T, k in L}  binary, >=0, <=1;
var theta{k in L, l in L}  binary, >=0, <=1; 
var lambda{i in S, j in S, l in L: i<>j} binary;	
var y{t in T, k in L}  >=0, <=U_tau;
var r{i in S, j in S, l in L: i<>j} >=0, <=U_tau;
var sigma{i in S, j in S, l in L} >= 0;		


################################################################
##############OBJECTIVE FUNCTIONS
minimize OBJ: sum{t in T, l in L} h[t,l]*u[t,l] + eta + (1/(1-alpha))*sum{i in S}p[i]*v[i] + (1/(1-alpha))*kappa*tau;

param RCVAR_weight default 0.5; 					#relative significance weight coefficient for Pareto analysis
minimize OBJ_Pareto: (1-RCVAR_weight)*(sum{t in T, l in L} h[t,l]*u[t,l]) + RCVAR_weight*(eta + (1/(1-alpha))*sum{i in S}p[i]*v[i] + (1/(1-alpha))*kappa*tau);


minimize RobustifiedCVaR: eta + (1/(1-alpha))*sum{i in S}p[i]*v[i] + (1/(1-alpha))*kappa*tau;


################################################################
##############CONSTRAINTS



subject to C_RobustificationLinear{i in S, j in S}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi[k,j,t]*z[t,k,l]) - eta - sum{l in L} sigma[i,j,l];

subject to C_Linearization_z_1{l in L, k in L, t in T}: 	
z[t,k,l] <= u[t,k];

subject to C_Linearization_z_2{l in L, k in L, t in T}: 
z[t,k,l] <= theta[k,l];

subject to C_Linearization_z_3{l in L, k in L, t in T}: 
z[t,k,l] >= u[t,k] + theta[k,l] - 1;


##McCormick inequalities associated with C_RobustificationLinear

#subject to C_varsigma_1{i in S, j in S, l in L}:  
#sigma[i,j,l] <=  sum{t in T} (y[t,l] *(xi[l,i,t]- xi[l,j,t])) + M[i,j,l]*r[i,j,l];

#subject to C_varsigma_2{i in S, j in S, l in L}:  
#sigma[i,j,l] <= sum{t in T} (y[t,l] * (-xi[l,i,t]+ xi[l,j,t])) + M[i,j,l]*(tau-r[i,j,l]);

##for the first compact formulation; later we improve it by using the symmetry: sigma[i,j,l]=sigma[j,i,l], and dropping the variables sigma[i,i,l] and sigma[i,j,l] for i>j
subject to C_varsigma_1_compact{i in S, j in S, l in L: i<>j}: 
sigma[i,j,l] <=  sum{t in T} (y[t,l] *(xi[l,i,t]- xi[l,j,t])) + M[i,j,l]*r[i,j,l];
 
subject to C_varsigma_2_compact{i in S, j in S, l in L: i<>j}: 
sigma[i,j,l] <= sum{t in T} (y[t,l] * (-xi[l,i,t]+ xi[l,j,t])) + M[i,j,l]*(tau-r[i,j,l]);

subject to C_sigma_fix{i in S, l in L}:
sigma[i,i,l] = 0;



###Type 1: linearizing the terms of the form u[t,l] * tau 
subject to C_LinearizationTypeI_1{l in L, t in T}:  
y[t,l] <= U_tau * u[t,l];

subject to C_LinearizationTypeI_2{l in L, t in T}:  
y[t,l] <= tau;

subject to C_LinearizationTypeI_3{l in L, t in T}:  
y[t,l] >= U_tau * u[t,l] + tau - U_tau;


###Type 2: linearizing the terms of the form lambda[i,j,l]*tau 
#subject to C_LinearizationTypeII_1{i in S, j in S, l in L}:  
#r[i,j,l] <= tau;

#subject to C_LinearizationTypeII_2{i in S, j in S, l in L}:  
#r[i,j,l] <= U_tau * lambda[i,j,l];

#subject to C_LinearizationTypeII_3{i in S, j in S, l in L}:  
#r[i,j,l] >= U_tau * lambda[i,j,l] + tau - U_tau;

###Type 2: linearizing the terms of the form lambda[i,j,l]*tau 
subject to C_LinearizationTypeII_1_compact{i in S, j in S, l in L: i<>j}:
r[i,j,l] <= tau;

subject to C_LinearizationTypeII_2_compact{i in S, j in S, l in L: i<>j}:
r[i,j,l] <= U_tau * lambda[i,j,l];

subject to C_LinearizationTypeII_3_compact{i in S, j in S, l in L: i<>j}:
r[i,j,l] >= U_tau * lambda[i,j,l] + tau - U_tau;


################################################################


###feasible schedules
subject to C_LinearOrdering_1{l in L}: 
theta[l,l] = 1;

subject to C_LinearOrdering_2{l in L, k in L: k < l}: 
theta[k,l] + theta[l,k] = 1;

subject to C_LinearOrdering_3{l in L, k in L, h1 in L: k < l < h1}: 
theta[k,l] + theta[l,h1] + theta[h1,k] <=2;

###feasible control decisions
subject to C_Control{l in L}: 
sum{t in T} u[t,l] = 1;

###this one is to get results for the decision-independent case
subject to C_Control_fixing{l in L}: 
u[1,l] = 1;


subject to C_Budget:
sum{t in T, l in L} h[t,l]*u[t,l] <= B;


################################################################################
###SOME ENHANCEMENTS FOR THE NON-COMONOTONE MIP FORMULATION

###drop the variables sigma[i,i,l] (also drop constraint C_sigma_fix)
###use the symmetry: sigma[i,j,l]=sigma[j,i,l], drop the variables sigma[i,j,l] for i>j (and exclude the corresponding absolute value constraints)


var lambda_sym{i in S, j in S, l in L: i<j} binary;	
var r_sym{i in S, j in S, l in L: i<j} >=0, <=U_tau;
var sigma_sym{i in S, j in S, l in L:i<j} >= 0;

subject to C_RobustificationLinear_compact{i in S}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi[k,i,t]*z[t,k,l]) - eta;

subject to C_RobustificationLinear_sym1{i in S, j in S: i<j}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi[k,j,t]*z[t,k,l]) - eta - sum{l in L} sigma_sym[i,j,l];

subject to C_RobustificationLinear_sym2{i in S, j in S: i>j}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi[k,j,t]*z[t,k,l]) - eta - sum{l in L} sigma_sym[j,i,l];

subject to C_varsigma_1_sym{i in S, j in S, l in L: i<j}: 
sigma_sym[i,j,l] <=  sum{t in T} (y[t,l] *(xi[l,i,t]- xi[l,j,t])) + M[i,j,l]*r_sym[i,j,l];
 
subject to C_varsigma_2_sym{i in S, j in S, l in L: i<j}: 
sigma_sym[i,j,l] <= sum{t in T} (y[t,l] * (-xi[l,i,t]+ xi[l,j,t])) + M[i,j,l]*(tau-r_sym[i,j,l]);


subject to C_LinearizationTypeII_1_sym{i in S, j in S, l in L: i<j}:
r_sym[i,j,l] <= tau;

subject to C_LinearizationTypeII_2_sym{i in S, j in S, l in L: i<j}:
r_sym[i,j,l] <= U_tau * lambda_sym[i,j,l];

subject to C_LinearizationTypeII_3_sym{i in S, j in S, l in L: i<j}:
r_sym[i,j,l] >= U_tau * lambda_sym[i,j,l] + tau - U_tau;

##Problem Description:
##Model_Basic_MIP_Sym: 
#problem Model_Basic_MIP_Sym: eta,v,tau,z,u,theta,lambda_sym,y,r_sym,sigma_sym,OBJ, C_RobustificationLinear_compact,C_RobustificationLinear_sym1,C_RobustificationLinear_sym2,
#C_Linearization_z_1,C_Linearization_z_2,C_Linearization_z_3, C_varsigma_1_sym, C_varsigma_2_sym,
#C_LinearizationTypeI_1,C_LinearizationTypeI_2,C_LinearizationTypeI_3,C_LinearizationTypeII_1_sym,C_LinearizationTypeII_2_sym,C_LinearizationTypeII_3_sym,
#C_Control, C_LinearOrdering_1,C_LinearOrdering_2, C_LinearOrdering_3;

###
#############################
###############HYBRID VERSION For NON-COMONOTONE Model
####benefit from the comonotone structure whenever possible; eliminate the indices for (i,j,l) in I_+ and I_-

set Iplus default {i in S, j in S, l in L: xi[l,i,1] >= xi[l,j,1] and xi[l,i,2] >= xi[l,j,2]};
set Iminus default {i in S, j in S, l in L: xi[l,i,1] <= xi[l,j,1] and xi[l,i,2] <= xi[l,j,2]};
set I0 default {i in S, j in S, l in L} diff {Iplus union Iminus};

var sigma_hybrid{(i,j,l) in I0: i<j} >= 0;		
var lambda_hybrid{(i,j,l) in I0: i<j} binary;	
var r_hybrid{(i,j,l) in I0: i<j} >=0, <=U_tau;


subject to C_RobustificationLinear_hybrid1{i in S, j in S: i<j}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi[k,j,t]*z[t,k,l]) - eta - (sum{l in L: (i,j,l) in Iplus} (sum{t in T} (y[t,l]*(xi[l,i,t]- xi[l,j,t]))) +
sum{l in L: (i,j,l) in Iminus} (sum{t in T} (y[t,l]*(-xi[l,i,t]+ xi[l,j,t]))) + sum{l in L: (i,j,l) in I0} sigma_hybrid[i,j,l]);

subject to C_RobustificationLinear_hybrid2{i in S, j in S: i>j}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi[k,j,t]*z[t,k,l]) - eta - (sum{l in L: (i,j,l) in Iplus} (sum{t in T} (y[t,l]*(xi[l,i,t]- xi[l,j,t]))) +
sum{l in L: (i,j,l) in Iminus} (sum{t in T} (y[t,l]*(-xi[l,i,t]+ xi[l,j,t]))) + sum{l in L: (i,j,l) in I0}sigma_hybrid[j,i,l]);

#subject to C_varsigma_compact_hybrid_a{(i,j,l) in Iplus}: 
#sigma[i,j,l] =  sum{t in T} (y[t,l] *(xi[l,i,t]- xi[l,j,t])); #plugged the values in the above constraints

#subject to C_varsigma_compact_hybrid_b{i in S, j in S, l in L: xi[l,i,1] <= xi[l,j,1] and xi[l,i,2] <= xi[l,j,2] and i<>j}: 
#sigma[i,j,l] = sum{t in T} (y[t,l] * (-xi[l,i,t]+ xi[l,j,t])); #plugged the values in the above constraints


subject to C_varsigma_1_compact_hybrid{(i,j,l) in I0: i<j}: 
sigma_hybrid[i,j,l] <=  sum{t in T} (y[t,l] *(xi[l,i,t]- xi[l,j,t])) + M[i,j,l]*r_hybrid[i,j,l];

subject to C_varsigma_2_compact_hybrid{(i,j,l) in I0: i<j}: 
sigma_hybrid[i,j,l] <= sum{t in T} (y[t,l] * (-xi[l,i,t]+ xi[l,j,t])) + M[i,j,l]*(tau-r_hybrid[i,j,l]);


subject to C_LinearizationTypeII_1_hybrid{(i,j,l) in I0: i<j}:
r_hybrid[i,j,l] <= tau;

subject to C_LinearizationTypeII_2_hybrid{(i,j,l) in I0: i<j}:
r_hybrid[i,j,l] <= U_tau * lambda_hybrid[i,j,l];

subject to C_LinearizationTypeII_3_hybrid{(i,j,l) in I0: i<j}:
r_hybrid[i,j,l] >= U_tau * lambda_hybrid[i,j,l] + tau - U_tau;


##Problem Description:
###Sym + Hybrid
#problem Model_Basic_MIP_Hybrid: eta,v,tau,z,u,theta,lambda_hybrid,y,r_hybrid,sigma_hybrid,OBJ, C_RobustificationLinear_compact,C_RobustificationLinear_hybrid1,C_RobustificationLinear_hybrid2,
#C_Linearization_z_1,C_Linearization_z_2,C_Linearization_z_3, C_varsigma_1_compact_hybrid, C_varsigma_2_compact_hybrid,
#C_LinearizationTypeI_1,C_LinearizationTypeI_2,C_LinearizationTypeI_3,C_LinearizationTypeII_1_hybrid,C_LinearizationTypeII_2_hybrid,C_LinearizationTypeII_3_hybrid,
#C_Control, C_LinearOrdering_1,C_LinearOrdering_2, C_LinearOrdering_3;


################################################################################
###SCHEDULING PROBLEM WITH THE COMONOTONE STRUCTURE

subject to C_Robustification_Comonotone{i in S, j in S}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi_base[k,j]*a[t,k]*z[t,k,l]) - eta - sum{l in L,t in T} abs(xi_base[l,i]- xi_base[l,j])*a[t,l]*y[t,l];

#subject to C_Linearization_z_1{l in L, k in L, t in T}: 	
#z[t,k,l] <= u[t,k];

#subject to C_Linearization_z_2{l in L, k in L, t in T}: 
#z[t,k,l] <= theta[k,l];

#subject to C_Linearization_z_3{l in L, k in L, t in T}: 
#z[t,k,l] >= u[t,k] + theta[k,l] - 1;

###Type 1: linearizing the terms of the form u[t,l] * tau 
#subject to C_LinearizationTypeI_1{l in L, t in T}:  
#y[t,l] <= U_tau * u[t,l];

#subject to C_LinearizationTypeI_2{l in L, t in T}:  
#y[t,l] <= tau;

#subject to C_LinearizationTypeI_3{l in L, t in T}:  
#y[t,l] >= U_tau * u[t,l] + tau - U_tau;

###feasible schedules
#subject to C_LinearOrdering_1{l in L}: 
#theta[l,l] = 1;

#subject to C_LinearOrdering_2{l in L, k in L: k < l}: 
#theta[k,l] + theta[l,k] = 1;

#subject to C_LinearOrdering_3{l in L, k in L, h1 in L: k < l < h1}: 
#theta[k,l] + theta[l,h1] + theta[h1,k] <=2;

###feasible control decisions
#subject to C_Control{l in L}: 
#sum{t in T} u[t,l] = 1;

#subject to C_Budget:
#sum{t in T, l in L} h[t,l]*u[t,l] <= B;

subject to C_ValidIneq1{l in L}:
sum{t in T} y[t,l] =tau; 

subject to C_ValidIneq2{k in L,l in L}:
sum{t in T} z[t,k,l]=theta[k,l];



######################################################
#######additional parameters
param temp_ampl_user default 0;
param temp_ampl_elapsed default 0;
param temp_total_solve_user default 0;
param temp_total_solve_elapsed default 0;
param com_measure1 default 0; #quantifies the comonotonicity level (Kendall's tau)
param com_measure2 default 0; #quantifies the comonotonicity level 

set SelectedModels default {1..12};
param M_Cost{i in SelectedModels, 1..2} default 0;
param M_INDICATE{i in SelectedModels} default 0;
param M_tau{i in SelectedModels} default 0;
param M_OBFV{i in SelectedModels} default 0;
param M_Vectheta{i in SelectedModels, k in L, l in L} default 0;
param M_Sequence{i in SelectedModels,l in L} default 0;
param M_Vecu{i in SelectedModels, t in T, k in L} default 0;
param M_Time{i in SelectedModels,1..4} default 0;
param M_ResNum{i in SelectedModels} default 0;
param M_Bestnode{i in SelectedModels} default 0;
param M_Message{i in SelectedModels} symbolic;
param Name{i in SelectedModels} symbolic;

#######################################################
 
##SOME QUADRATIC formulations; not used [can be ignored]; they can be solved by some nonlinear optimization solvers [tried a few, not promising]
var nu{i in S, j in S, l in L} >= 0;		#introduced for the quadratic version (robustification constraint involves the terms of the form nu[i,j,l]*tau)

subject to C_nu_fix{i in S, l in L}:
nu[i,i,l] = 0;

###Robustification constraint with quadratic terms of the form nu[i,j,l]*tau
subject to C_RobustificationQuad_1{i in S, j in S}:
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi[k,j,t]*z[t,k,l]) - eta - sum{l in L} nu[i,j,l]*tau;

###Robustification constraint with quadratic terms of the form u[t,k]*theta[k,l] and nu[i,j,l]*tau
subject to C_RobustificationFullQuad{i in S, j in S}:
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi[k,j,t]*u[t,k]*theta[k,l]) - eta - sum{l in L} nu[i,j,l]*tau;

subject to C_RobustificationQuad_2{i in S, j in S, l in L: i<>j}:
nu[i,j,l] <= sum{t in T} u[t,l]*(xi[l,i,t]-xi[l,j,t]) + M[i,j,l]*lambda[i,j,l];

subject to C_RobustificationQuad_3{i in S, j in S, l in L: i<>j}:
nu[i,j,l]<= -sum{t in T} u[t,l]*(xi[l,i,t]-xi[l,j,t]) + M[i,j,l]*(1-lambda[i,j,l]);

###SCHEDULING PROBLEM WITH THE COMONOTONE STRUCTURE

####Robustification constraint with quadratic terms of the form u[t,l]*tau
subject to C_Robustification_Comonotone_Quad{i in S, j in S}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi_base[k,j]*a[t,k]*z[t,k,l]) - eta - sum{l in L,t in T} abs(xi_base[l,i]- xi_base[l,j])*a[t,l]*u[t,l]*tau;

####Robustification constraint with quadratic terms of the form u[t,k]*theta[k,l] and u[t,l]*tau
subject to C_Robustification_ComonotoneFullQuad{i in S, j in S}: 
v[i] >= sum{l in L, k in L, t in T}(w[l]*xi_base[k,j]*a[t,k]*u[t,k]*theta[k,l]) - eta - sum{l in L,t in T} abs(xi_base[l,i]- xi_base[l,j])*a[t,l]*u[t,l]*tau;

 

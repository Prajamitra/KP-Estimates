


################ MSM in Manzini and Mbabane with Jeffrey's Prior for Ns in MSM  ####################


rm(list=ls())

################### R packages ##############

library(ggplot2)
library(MCMCpack)
library(HDInterval)


#############################################

tot=500000

th=20
lw=tot*0.5
up=tot


########## Initial parameter specification ###############


####### For Manzini ###############

N_2=450

beta01_2=0.2
beta02_2=0.4
beta03_2=0.3

p_1_2=0.6
p_2_2=0.45
p_3_2=0.3

####### For Mbabane ###############

N_3=475

beta01_3=0.25
beta02_3=0.15
beta03_3=0.45

p_1_3=0.65
p_2_3=0.15
p_3_3=0.50


#### Hyper parameters of the priors

## For Manzini 

a1_2=0.5
b1_2=0.5

a2_2=0.5
b2_2=0.5

a3_2=0.5
b3_2=0.5

beta_prior_2=array(0.5,4) # beta parameters are 0.5 refer Jeffrey's Dirichlet prior

## For Mbabane 

a1_3=0.5
b1_3=0.5

a2_3=0.5
b2_3=0.5

a3_3=0.5
b3_3=0.5

beta_prior_3=array(0.5,4) # beta parameters are 0.5 refer Jeffrey's Dirichlet prior


NN_2=0.8*N_2


##########################  Data ################################


###### TRS (MANZINI) ########
x_111_2=16
x_101_2=60
x_011_2=8
x_001_2=93
x_1.._2=194

x_1.0_2=x_1.._2-x_111_2-x_101_2

x_.11_2=x_111_2+x_011_2

x_.01_2=x_101_2+x_001_2

x_..1_2=x_.11_2+x_.01_2

x0_star_2=x_111_2+x_101_2+x_011_2+x_001_2

##### TRS (MBABANE) ##########
x_111_3=21
x_101_3=81
x_011_3=8
x_001_3=113
x_1.._3=217

x_1.0_3=x_1.._3-x_111_3-x_101_3

x_.11_3=x_111_3+x_011_3

x_.01_3=x_101_3+x_001_3

x_..1_3=x_.11_3+x_.01_3

x0_star_3=x_111_3+x_101_3+x_011_3+x_001_3

NN_3=0.8*N_3

x_.1._0=106

K=x_.1._0-x_.11_2-x_.11_3



########## Gibbs ######################################

## For Manzini

N_2_estimate=array(NA,tot)

p_2=array(NA,c(3,tot))

y111_2=array(NA,tot)

y110_1_2=array(NA,tot)
y110_2_2=array(NA,tot)
y110_3_2=array(NA,tot)

y011_1_2=array(NA,tot)
y011_2_2=array(NA,tot)
y011_3_2=array(NA,tot)

y100_1_2=array(NA,tot)
y100_2_2=array(NA,tot)
y100_3_2=array(NA,tot)

y101_1_2=array(NA,tot)
y101_2_2=array(NA,tot)
y101_3_2=array(NA,tot)

y010_1_2=array(NA,tot)
y010_2_2=array(NA,tot)
y010_3_2=array(NA,tot)

y001_1_2=array(NA,tot)
y001_2_2=array(NA,tot)
y001_3_2=array(NA,tot)

y000_2=array(NA,tot)

z_010_2=array(NA,tot)
z_110_2=array(NA,tot)
z_100_2=array(NA,tot)

## For Mbabane

N_3_estimate=array(NA,tot)

p_3=array(NA,c(3,tot))

p_111_2=array(NA,tot)
p_110_2=array(NA,tot)
p_101_2=array(NA,tot)
p_011_2=array(NA,tot)
p_100_2=array(NA,tot)
p_001_2=array(NA,tot)
p_010_2=array(NA,tot)
p_000_2=array(NA,tot)

p_111_3=array(NA,tot)
p_110_3=array(NA,tot)
p_101_3=array(NA,tot)
p_011_3=array(NA,tot)
p_100_3=array(NA,tot)
p_001_3=array(NA,tot)
p_010_3=array(NA,tot)
p_000_3=array(NA,tot)

y111_3=array(NA,tot)

y110_1_3=array(NA,tot)
y110_2_3=array(NA,tot)
y110_3_3=array(NA,tot)

y011_1_3=array(NA,tot)
y011_2_3=array(NA,tot)
y011_3_3=array(NA,tot)

y100_1_3=array(NA,tot)
y100_2_3=array(NA,tot)
y100_3_3=array(NA,tot)

y101_1_3=array(NA,tot)
y101_2_3=array(NA,tot)
y101_3_3=array(NA,tot)

y010_1_3=array(NA,tot)
y010_2_3=array(NA,tot)
y010_3_3=array(NA,tot)

y001_1_3=array(NA,tot)
y001_2_3=array(NA,tot)
y001_3_3=array(NA,tot)

y000_3=array(NA,tot)

z_010_3=array(NA,tot)
z_110_3=array(NA,tot)
z_100_3=array(NA,tot)

beta_1_2=array(NA,tot)
beta_2_2=array(NA,tot)
beta_3_2=array(NA,tot)
beta_0_2=array(NA,tot)

beta_1_3=array(NA,tot)
beta_2_3=array(NA,tot)
beta_3_3=array(NA,tot)
beta_0_3=array(NA,tot)



# Iteration 1 values for Manzini

beta_1_2[1]=beta01_2
beta_2_2[1]=beta02_2
beta_3_2[1]=beta03_2 

beta_0_2[1]=min((beta_1_2[1]+beta_2_2[1]+beta_3_2[1]),1)

p_2[1,1]=p_1_2
p_2[2,1]=p_2_2
p_2[3,1]=p_3_2

N_2_estimate[1]=N_2

# Iteration 1 values for Mbabane

beta_1_3[1]=beta01_3
beta_2_3[1]=beta02_3
beta_3_3[1]=beta03_3 

beta_0_3[1]=min((beta_1_3[1]+beta_2_3[1]+beta_3_3[1]),1)

p_3[1,1]=p_1_3
p_3[2,1]=p_2_3
p_3[3,1]=p_3_3

N_3_estimate[1]=N_3

z_110_2[1]=round(min(K,x_1.0_2)/3)
z_100_2[1]=round(x_1.0_2/2)


z_110_3[1]=round(min(K,x_1.0_3)/3)
z_100_3[1]=round(x_1.0_3/2)


# for loop starts for Gibbs sampling

for(h in 2:tot){
  
  #### Manzini
  
  p_111_2[h]=((1-beta_0_2[h-1])*p_2[1,h-1]*p_2[2,h-1]*p_2[3,h-1])
  
  p_110_2[h]=((1-beta_0_2[h-1])*p_2[1,h-1]*p_2[2,h-1]*(1-p_2[3,h-1]))+(beta_2_2[h-1]*p_2[1,h-1]*p_2[2,h-1])+(beta_3_2[h-1]*p_2[1,h-1]*p_2[2,h-1])
  Q_110_1_2=((1-beta_0_2[h-1])*p_2[1,h-1]*p_2[2,h-1]*(1-p_2[3,h-1]))/p_110_2[h]
  Q_110_2_2=(beta_2_2[h-1]*p_2[1,h-1]*p_2[2,h-1])/p_110_2[h]
  Q_110_3_2=max(0,1-Q_110_1_2-Q_110_2_2)
  prob_y110_2=c(Q_110_1_2,Q_110_2_2,Q_110_3_2)
  
  
  p_011_2[h]=((1-beta_0_2[h-1])*(1-p_2[1,h-1])*p_2[2,h-1]*p_2[3,h-1])+(beta_1_2[h-1]*(1-p_2[1,h-1])*p_2[3,h-1])+(beta_3_2[h-1]*(1-p_2[1,h-1])*p_2[2,h-1])
  Q_011_1_2=((1-beta_0_2[h-1])*(1-p_2[1,h-1])*p_2[2,h-1]*p_2[3,h-1])/p_011_2[h]
  Q_011_2_2=(beta_1_2[h-1]*(1-p_2[1,h-1])*p_2[3,h-1])/p_011_2[h]
  Q_011_3_2=max(0,1-Q_011_1_2-Q_011_2_2)
  prob_y011_2=c(Q_011_1_2,Q_011_2_2,Q_011_3_2)
  
  
  p_100_2[h]=((1-beta_0_2[h-1])*p_2[1,h-1]*(1-p_2[2,h-1])*(1-p_2[3,h-1]))+(beta_1_2[h-1]*p_2[1,h-1]*(1-p_2[3,h-1]))+(beta_3_2[h-1]*p_2[1,h-1]*(1-p_2[2,h-1]))
  Q_100_1_2=((1-beta_0_2[h-1])*p_2[1,h-1]*(1-p_2[2,h-1])*(1-p_2[3,h-1]))/p_100_2[h]
  Q_100_2_2=(beta_1_2[h-1]*p_2[1,h-1]*(1-p_2[3,h-1]))/p_100_2[h]
  Q_100_3_2=max(0,1-Q_100_1_2-Q_100_2_2)
  prob_y100_2=c(Q_100_1_2,Q_100_2_2,Q_100_3_2)
  
  
  p_001_2[h]=((1-beta_0_2[h-1])*(1-p_2[1,h-1])*(1-p_2[2,h-1])*p_2[3,h-1])+(beta_2_2[h-1]*(1-p_2[1,h-1])*(1-p_2[2,h-1]))+(beta_3_2[h-1]*(1-p_2[1,h-1])*(1-p_2[2,h-1]))
  Q_001_1_2=((1-beta_0_2[h-1])*(1-p_2[1,h-1])*(1-p_2[2,h-1])*p_2[3,h-1])/p_001_2[h]
  Q_001_2_2=(beta_2_2[h-1]*(1-p_2[1,h-1])*(1-p_2[2,h-1]))/p_001_2[h]
  Q_001_3_2=max(0,1-Q_001_1_2-Q_001_2_2)
  prob_y001_2=c(Q_001_1_2,Q_001_2_2,Q_001_3_2)
  
  
  p_101_2[h]=((1-beta_0_2[h-1])*p_2[1,h-1]*(1-p_2[2,h-1])*p_2[3,h-1])+(beta_1_2[h-1]*p_2[1,h-1]*p_2[3,h-1])+(beta_2_2[h-1]*p_2[1,h-1]*(1-p_2[2,h-1]))
  Q_101_1_2=(1-beta_0_2[h-1])*p_2[1,h-1]*(1-p_2[2,h-1])*p_2[3,h-1]/p_101_2[h]
  Q_101_2_2=(beta_1_2[h-1]*p_2[1,h-1]*p_2[3,h-1])/p_101_2[h]
  Q_101_3_2=max(0,1-Q_101_1_2-Q_101_2_2)
  prob_y101_2=c(Q_101_1_2,Q_101_2_2,Q_101_3_2)
  
  
  p_010_2[h]=((1-beta_0_2[h-1])*(1-p_2[1,h-1])*p_2[2,h-1]*(1-p_2[3,h-1]))+(beta_1_2[h-1]*(1-p_2[1,h-1])*(1-p_2[3,h-1]))+(beta_2_2[h-1]*(1-p_2[1,h-1])*p_2[2,h-1])
  Q_010_1_2=((1-beta_0_2[h-1])*(1-p_2[1,h-1])*p_2[2,h-1]*(1-p_2[3,h-1]))/p_010_2[h]
  Q_010_2_2=(beta_1_2[h-1]*(1-p_2[1,h-1])*(1-p_2[3,h-1]))/p_010_2[h]
  Q_010_3_2=max(0,1-Q_010_1_2-Q_010_2_2)
  prob_y010_2=c(Q_010_1_2,Q_010_2_2,Q_010_3_2)
  
  
  p_000_2[h]=(1-beta_0_2[h-1])*(1-p_2[1,h-1])*(1-p_2[2,h-1])*(1-p_2[3,h-1])
  
  
  
  
  
  
  #### Manzini
  
  p_111_3[h]=((1-beta_0_3[h-1])*p_3[1,h-1]*p_3[2,h-1]*p_3[3,h-1])
  
  p_110_3[h]=((1-beta_0_3[h-1])*p_3[1,h-1]*p_3[2,h-1]*(1-p_3[3,h-1]))+(beta_2_3[h-1]*p_3[1,h-1]*p_3[2,h-1])+(beta_3_3[h-1]*p_3[1,h-1]*p_3[2,h-1])
  Q_110_1_3=((1-beta_0_3[h-1])*p_3[1,h-1]*p_3[2,h-1]*(1-p_3[3,h-1]))/p_110_3[h]
  Q_110_2_3=(beta_2_3[h-1]*p_3[1,h-1]*p_3[2,h-1])/p_110_3[h]
  Q_110_3_3=max(0,1-Q_110_1_3-Q_110_2_3)
  prob_y110_3=c(Q_110_1_3,Q_110_2_3,Q_110_3_3)
  
  
  p_011_3[h]=((1-beta_0_3[h-1])*(1-p_3[1,h-1])*p_3[2,h-1]*p_3[3,h-1])+(beta_1_3[h-1]*(1-p_3[1,h-1])*p_3[3,h-1])+(beta_3_3[h-1]*(1-p_3[1,h-1])*p_3[2,h-1])
  Q_011_1_3=((1-beta_0_3[h-1])*(1-p_3[1,h-1])*p_3[2,h-1]*p_3[3,h-1])/p_011_3[h]
  Q_011_2_3=(beta_1_3[h-1]*(1-p_3[1,h-1])*p_3[3,h-1])/p_011_3[h]
  Q_011_3_3=max(0,1-Q_011_1_3-Q_011_2_3)
  prob_y011_3=c(Q_011_1_3,Q_011_2_3,Q_011_3_3)
  
  
  p_100_3[h]=((1-beta_0_3[h-1])*p_3[1,h-1]*(1-p_3[2,h-1])*(1-p_3[3,h-1]))+(beta_1_3[h-1]*p_3[1,h-1]*(1-p_3[3,h-1]))+(beta_3_3[h-1]*p_3[1,h-1]*(1-p_3[2,h-1]))
  Q_100_1_3=((1-beta_0_3[h-1])*p_3[1,h-1]*(1-p_3[2,h-1])*(1-p_3[3,h-1]))/p_100_3[h]
  Q_100_2_3=(beta_1_3[h-1]*p_3[1,h-1]*(1-p_3[3,h-1]))/p_100_3[h]
  Q_100_3_3=max(0,1-Q_100_1_3-Q_100_2_3)
  prob_y100_3=c(Q_100_1_3,Q_100_2_3,Q_100_3_3)
  
  
  p_001_3[h]=((1-beta_0_3[h-1])*(1-p_3[1,h-1])*(1-p_3[2,h-1])*p_3[3,h-1])+(beta_2_3[h-1]*(1-p_3[1,h-1])*(1-p_3[2,h-1]))+(beta_3_3[h-1]*(1-p_3[1,h-1])*(1-p_3[2,h-1]))
  Q_001_1_3=((1-beta_0_3[h-1])*(1-p_3[1,h-1])*(1-p_3[2,h-1])*p_3[3,h-1])/p_001_3[h]
  Q_001_2_3=(beta_2_3[h-1]*(1-p_3[1,h-1])*(1-p_3[2,h-1]))/p_001_3[h]
  Q_001_3_3=max(0,1-Q_001_1_3-Q_001_2_3)
  prob_y001_3=c(Q_001_1_3,Q_001_2_3,Q_001_3_3)
  
  
  p_101_3[h]=((1-beta_0_3[h-1])*p_3[1,h-1]*(1-p_3[2,h-1])*p_3[3,h-1])+(beta_1_3[h-1]*p_3[1,h-1]*p_3[3,h-1])+(beta_2_3[h-1]*p_3[1,h-1]*(1-p_3[2,h-1]))
  Q_101_1_3=(1-beta_0_3[h-1])*p_3[1,h-1]*(1-p_3[2,h-1])*p_3[3,h-1]/p_101_3[h]
  Q_101_2_3=(beta_1_3[h-1]*p_3[1,h-1]*p_3[3,h-1])/p_101_3[h]
  Q_101_3_3=max(0,1-Q_101_1_3-Q_101_2_3)
  prob_y101_3=c(Q_101_1_3,Q_101_2_3,Q_101_3_3)
  
  
  p_010_3[h]=((1-beta_0_3[h-1])*(1-p_3[1,h-1])*p_3[2,h-1]*(1-p_3[3,h-1]))+(beta_1_3[h-1]*(1-p_3[1,h-1])*(1-p_3[3,h-1]))+(beta_2_3[h-1]*(1-p_3[1,h-1])*p_3[2,h-1])
  Q_010_1_3=((1-beta_0_3[h-1])*(1-p_3[1,h-1])*p_3[2,h-1]*(1-p_3[3,h-1]))/p_010_3[h]
  Q_010_2_3=(beta_1_3[h-1]*(1-p_3[1,h-1])*(1-p_3[3,h-1]))/p_010_3[h]
  Q_010_3_3=max(0,1-Q_010_1_3-Q_010_2_3)
  prob_y010_3=c(Q_010_1_3,Q_010_2_3,Q_010_3_3)
  
  
  p_000_3[h]=(1-beta_0_3[h-1])*(1-p_3[1,h-1])*(1-p_3[2,h-1])*(1-p_3[3,h-1])
  
  
  ######  Generation of Z_010 for Manzini & Mbabane 
  
  P010_2_star=p_010_2[h]/(p_010_2[h]+p_000_2[h])
  P010_3_star=p_010_3[h]/(p_010_3[h]+p_000_3[h])
  
  K1=K-z_110_2[h-1]-z_110_3[h-1]
  
  size_param_z_010_2=N_2_estimate[h-1]-x_..1_2-z_110_2[h-1]-z_100_2[h-1]
  
  size_param_z_010_3=N_3_estimate[h-1]-x_..1_3-z_110_3[h-1]-z_100_3[h-1]
  
  lower_z_010_2=max(0,K1-size_param_z_010_3)
  upper_z_010_2=min(size_param_z_010_2,K1)
  
  #lower_z_010_3=0
  #upper_z_010_3
  
  prob2=c()
  
  for(t in lower_z_010_2:upper_z_010_2){
    prob2=c(prob2,(dbinom(t,size_param_z_010_2,P010_2_star)*dbinom((K1-t),size_param_z_010_3,P010_3_star)))   
  }
  
  prob2_star=prob2/sum(prob2)
  
  if (lower_z_010_2==upper_z_010_2) {
    z_010_2[h]=upper_z_010_2 
  } else {
    z_010_2[h]=sample(lower_z_010_2:upper_z_010_2,1,replace =TRUE,prob2_star)
  }
  
  z_010_3[h]=K1-z_010_2[h]
  
  
  ######  Generation of Z_110 for Manzini & Mbabane
  
  P110_2_star=p_110_2[h]/(p_110_2[h]+p_100_2[h])
  P110_3_star=p_110_3[h]/(p_110_3[h]+p_100_3[h])
  
  K2=K-z_010_2[h]-z_010_3[h]
  
  size_param_z_110_2=x_1.0_2
  
  size_param_z_110_3=x_1.0_3
  
  lower_z_110_2=max(0,K2-x_1.0_3)
  upper_z_110_2=min(size_param_z_110_2,K2)
  
  #lower_z_110_3=0
  #upper_z_110_3
  
  prob3=c()
  
  for(v in lower_z_110_2:upper_z_110_2){
    prob3=c(prob3,(dbinom(v,size_param_z_110_2,P110_2_star)*dbinom((K2-v),size_param_z_110_3,P110_3_star)))   
  }
  
  prob3_star=prob3/sum(prob3)
  
  if (lower_z_110_2==upper_z_110_2) {
    z_110_2[h]=upper_z_110_2 
  } else {
    z_110_2[h]=sample(lower_z_110_2:upper_z_110_2,1,replace =TRUE,prob3_star)
  }
  
  z_110_3[h]=K2-z_110_2[h]
  
  z_100_2[h]=x_1.0_2-z_110_2[h]
  
  z_100_3[h]=x_1.0_3-z_110_3[h]
  
  NN_2=x0_star_2 + z_100_2[h] + z_110_2[h] + z_010_2[h]
  
  NN_3=x0_star_3 + z_100_3[h] + z_110_3[h] + z_010_3[h]
  
  ###############################################
  
  y110_2_vec_draw=rmultinom(1,z_110_2[h],prob_y110_2)
  
  y110_1_2[h]=y110_2_vec_draw[1,1]
  y110_2_2[h]=y110_2_vec_draw[2,1]          
  y110_3_2[h]=z_110_2[h]-y110_1_2[h]-y110_2_2[h]
  
  
  
  y011_2_vec_draw=rmultinom(1,x_011_2,prob_y011_2)
  
  y011_1_2[h]=y011_2_vec_draw[1,1]
  y011_2_2[h]=y011_2_vec_draw[2,1]
  y011_3_2[h]=x_011_2-y011_1_2[h]-y011_2_2[h]
  
  
  
  y100_2_vec_draw=rmultinom(1,z_100_2[h],prob_y100_2)
  
  y100_1_2[h]=y100_2_vec_draw[1,1]
  y100_2_2[h]=y100_2_vec_draw[2,1]
  y100_3_2[h]=z_100_2[h]-y100_1_2[h]-y100_2_2[h]
  
  
  y101_2_vec_draw=rmultinom(1,x_101_2,prob_y101_2)
  
  y101_1_2[h]=y101_2_vec_draw[1,1]
  y101_2_2[h]=y101_2_vec_draw[2,1]
  y101_3_2[h]=x_101_2-y101_1_2[h]-y101_2_2[h]
  
  y010_2_vec_draw=rmultinom(1,z_010_2[h],prob_y010_2)
  
  y010_1_2[h]=y010_2_vec_draw[1,1]
  y010_2_2[h]=y010_2_vec_draw[2,1]
  y010_3_2[h]=z_010_2[h]-y010_1_2[h]-y010_2_2[h]
  
  
  y001_2_vec_draw=rmultinom(1,x_001_2,prob_y001_2)
  
  y001_1_2[h]=y001_2_vec_draw[1,1]
  y001_2_2[h]=y001_2_vec_draw[2,1]
  y001_3_2[h]=x_001_2-y001_1_2[h]-y001_2_2[h]
  
  ########################################
  
  y110_3_vec_draw=rmultinom(1,z_110_3[h],prob_y110_3)
  
  y110_1_3[h]=y110_3_vec_draw[1,1]
  y110_2_3[h]=y110_3_vec_draw[2,1]          
  y110_3_3[h]=z_110_3[h]-y110_1_3[h]-y110_2_3[h]
  
  
  
  y011_3_vec_draw=rmultinom(1,x_011_3,prob_y011_3)
  
  y011_1_3[h]=y011_3_vec_draw[1,1]
  y011_2_3[h]=y011_3_vec_draw[2,1]
  y011_3_3[h]=x_011_3-y011_1_3[h]-y011_2_3[h]
  
  
  
  y100_3_vec_draw=rmultinom(1,z_100_3[h],prob_y100_3)
  
  y100_1_3[h]=y100_3_vec_draw[1,1]
  y100_2_3[h]=y100_3_vec_draw[2,1]
  y100_3_3[h]=z_100_3[h]-y100_1_3[h]-y100_2_3[h]
  
  
  y101_3_vec_draw=rmultinom(1,x_101_3,prob_y101_3)
  
  y101_1_3[h]=y101_3_vec_draw[1,1]
  y101_2_3[h]=y101_3_vec_draw[2,1]
  y101_3_3[h]=x_101_3-y101_1_3[h]-y101_2_3[h]
  
  y010_3_vec_draw=rmultinom(1,z_010_3[h],prob_y010_3)
  
  y010_1_3[h]=y010_3_vec_draw[1,1]
  y010_2_3[h]=y010_3_vec_draw[2,1]
  y010_3_3[h]=z_010_3[h]-y010_1_3[h]-y010_2_3[h]
  
  y001_3_vec_draw=rmultinom(1,x_001_3,prob_y001_3)
  
  y001_1_3[h]=y001_3_vec_draw[1,1]
  y001_2_3[h]=y001_3_vec_draw[2,1]
  y001_3_3[h]=x_001_3-y001_1_3[h]-y001_2_3[h]
  
  ##########################################
  
  d1_2=y011_2_2[h]+y100_2_2[h]+y101_2_2[h]+y010_2_2[h]+beta_prior_2[1]
  d2_2=y110_2_2[h]+x_101_2-y101_1_2[h]-y101_2_2[h]+z_010_2[h]-y010_1_2[h]-y010_2_2[h]+y001_2_2[h]+beta_prior_2[2]
  d3_2=z_110_2[h]-y110_1_2[h]-y110_2_2[h]+x_011_2-y011_1_2[h]-y011_2_2[h]+z_100_2[h]-y100_1_2[h]-y100_2_2[h]+x_001_2-y001_1_2[h]-y001_2_2[h]+beta_prior_2[3]
  d4_2=x_111_2+(N_2_estimate[h-1]-NN_2)+(y110_1_2[h]+y011_1_2[h]+y100_1_2[h]+y101_1_2[h]+y010_1_2[h]+y001_1_2[h])+beta_prior_2[4]
  
  beta_vec_draw_2=rdirichlet(1,c(d1_2,d2_2,d3_2,d4_2))
  
  beta_1_2[h]=beta_vec_draw_2[1,1]
  beta_2_2[h]=beta_vec_draw_2[1,2]
  beta_3_2[h]=beta_vec_draw_2[1,3]  
  beta_0_2[h]=min((beta_1_2[h]+beta_2_2[h]+beta_3_2[h]),1)
  
  ##############################################
  
  d1_3=y011_2_3[h]+y100_2_3[h]+y101_2_3[h]+y010_2_3[h]+beta_prior_3[1]
  d2_3=y110_2_3[h]+x_101_3-y101_1_3[h]-y101_2_3[h]+z_010_3[h]-y010_1_3[h]-y010_2_3[h]+y001_2_3[h]+beta_prior_3[2]
  d3_3=z_110_3[h]-y110_1_3[h]-y110_2_3[h]+x_011_3-y011_1_3[h]-y011_2_3[h]+z_100_3[h]-y100_1_3[h]-y100_2_3[h]+x_001_3-y001_1_3[h]-y001_2_3[h]+beta_prior_3[3]
  d4_3=x_111_3+(N_3_estimate[h-1]-NN_3)+(y110_1_3[h]+y011_1_3[h]+y100_1_3[h]+y101_1_3[h]+y010_1_3[h]+y001_1_3[h])+beta_prior_3[4]
  
  beta_vec_draw_3=rdirichlet(1,c(d1_3,d2_3,d3_3,d4_3))
  
  beta_1_3[h]=beta_vec_draw_3[1,1]
  beta_2_3[h]=beta_vec_draw_3[1,2]
  beta_3_3[h]=beta_vec_draw_3[1,3]  
  beta_0_3[h]=min((beta_1_3[h]+beta_2_3[h]+beta_3_3[h]),1)
  
  ########################################################
  
  a_p1_2=x_111_2+z_110_2[h]+z_100_2[h]+x_101_2+a1_2
  b_p1_2=x_011_2+z_010_2[h]+x_001_2+N_2_estimate[h-1]-NN_2+b1_2
  
  a_p2_2=x_111_2+z_110_2[h]+x_011_2-y011_2_2[h]+z_010_2[h]-y010_2_2[h]+a2_2
  b_p2_2=z_100_2[h]-y100_2_2[h]+x_101_2-y101_2_2[h]+x_001_2+N_2_estimate[h-1]-NN_2+b2_2
  
  a_p3_2=x_111_2+y011_1_2[h]+y011_2_2[h]+y101_1_2[h]+y101_2_2[h]+y001_1_2[h]+a3_2
  b_p3_2=y110_1_2[h]+y100_1_2[h]+y100_2_2[h]+y010_1_2[h]+y010_2_2[h]+N_2_estimate[h-1]-NN_2+b3_2
  
  p_2[1,h]=rbeta(1,a_p1_2,b_p1_2)
  p_2[2,h]=rbeta(1,a_p2_2,b_p2_2)
  p_2[3,h]=rbeta(1,a_p3_2,b_p3_2)
  
  ############################################
  
  a_p1_3=x_111_3+z_110_3[h]+z_100_3[h]+x_101_3+a1_3
  b_p1_3=x_011_3+z_010_3[h]+x_001_3+N_3_estimate[h-1]-NN_3+b1_3
  
  a_p2_3=x_111_3+z_110_3[h]+x_011_3-y011_2_3[h]+z_010_3[h]-y010_2_3[h]+a2_3
  b_p2_3=z_100_3[h]-y100_2_3[h]+x_101_3-y101_2_3[h]+x_001_3+N_3_estimate[h-1]-NN_3+b2_3
  
  a_p3_3=x_111_3+y011_1_3[h]+y011_2_3[h]+y101_1_3[h]+y101_2_3[h]+y001_1_3[h]+a3_3
  b_p3_3=y110_1_3[h]+y100_1_3[h]+y100_2_3[h]+y010_1_3[h]+y010_2_3[h]+N_3_estimate[h-1]-NN_3+b3_3
  
  p_3[1,h]=rbeta(1,a_p1_3,b_p1_3)
  p_3[2,h]=rbeta(1,a_p2_3,b_p2_3)
  p_3[3,h]=rbeta(1,a_p3_3,b_p3_3)
  
  #############
  
  w_2=rnbinom(1,NN_2,(1-(1-beta_0_2[h])*(1-p_2[1,h])*(1-p_2[2,h])*(1-p_2[3,h])))
  N_2_estimate[h]=w_2+NN_2
  
  ############
  
  w_3=rnbinom(1,NN_3,(1-(1-beta_0_3[h])*(1-p_3[1,h])*(1-p_3[2,h])*(1-p_3[3,h])))
  N_3_estimate[h]=w_3+NN_3
  
} # end of the 'h' loop





###################################################

t=seq(lw,up,th)

M_2=N_2_estimate[seq(lw,up,th)]


N_2_repli_mean_estimates=round(mean(M_2))
N_2_se=sd(M_2)
N_2_repli_median_estimates=median(M_2)

N_2_Interval_hpd=hdi(M_2, credMass=0.95) #  95% highest posterior credible interval

beta1_2=median(beta_1_2[seq(lw,up, th)])
beta2_2=median(beta_2_2[seq(lw,up, th)])
beta3_2=median(beta_3_2[seq(lw,up, th)])
beta_ind_2=median(1-beta_0_2[seq(lw,up, th)])

p_1_2_estimate=median(p_2[1,seq(lw,up, th)])
p_2_2_estimate=median(p_2[2,seq(lw,up, th)])
p_3_2_estimate=median(p_2[3,seq(lw,up, th)])

p_1.._2=p_111_2[seq(lw,up, th)]+p_110_2[seq(lw,up, th)]+p_101_2[seq(lw,up, th)]+p_100_2[seq(lw,up, th)]
p_.1._2=p_111_2[seq(lw,up, th)]+p_110_2[seq(lw,up, th)]+p_011_2[seq(lw,up, th)]+p_010_2[seq(lw,up, th)]
p_..1_2=p_111_2[seq(lw,up, th)]+p_101_2[seq(lw,up, th)]+p_011_2[seq(lw,up, th)]+p_001_2[seq(lw,up, th)]

est_p_1.._2=median(p_1.._2)
est_p_.1._2=median(p_.1._2)
est_p_..1_2=median(p_..1_2)







M_3=N_3_estimate[seq(lw,up,th)]

N_3_repli_mean_estimates=round(mean(M_3))
N_3_se=sd(M_3)
N_3_repli_median_estimates=median(M_3)

N_3_Interval_hpd=hdi(M_3, credMass=0.95) #  95% highest posterior credible interval

beta1_3=median(beta_1_3[seq(lw,up, th)])
beta2_3=median(beta_2_3[seq(lw,up, th)])
beta3_3=median(beta_3_3[seq(lw,up, th)])
beta_ind_3=median(1-beta_0_3[seq(lw,up, th)])

p_1_3_estimate=median(p_3[1,seq(lw,up, th)])
p_2_3_estimate=median(p_3[2,seq(lw,up, th)])
p_3_3_estimate=median(p_3[3,seq(lw,up, th)])


p_1.._3=p_111_3[seq(lw,up, th)]+p_110_3[seq(lw,up, th)]+p_101_3[seq(lw,up, th)]+p_100_3[seq(lw,up, th)]
p_.1._3=p_111_3[seq(lw,up, th)]+p_110_3[seq(lw,up, th)]+p_011_3[seq(lw,up, th)]+p_010_3[seq(lw,up, th)]
p_..1_3=p_111_3[seq(lw,up, th)]+p_101_3[seq(lw,up, th)]+p_011_3[seq(lw,up, th)]+p_001_3[seq(lw,up, th)]

est_p_1.._3=median(p_1.._3)
est_p_.1._3=median(p_.1._3)
est_p_..1_3=median(p_..1_3)




##################### Trace Plots ########################

par(mgp = c(2, 1.5, 0))

plot(t, N_2_estimate[seq(lw, up, th)], xlab = "Iteration", ylab = "Population Size", type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p_2[1, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[1])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p_2[2, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[2])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p_2[3, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[3])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_1_2[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[1])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_2_2[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[2])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_3_2[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[3])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')

plot(t, N_3_estimate[seq(lw, up, th)], xlab = "Iteration", ylab = "Population Size", type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p_3[1, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[1])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p_3[2, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[2])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p_3[3, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[3])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_1_3[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[1])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_2_3[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[2])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_3_3[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[3])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')






##################### Results ########################################

############### Results of Manzini #####################


#########  KP Sizes ###########

N_2_repli_mean_estimates
N_2_repli_median_estimates
N_2_se

N_2_Interval_hpd[1]  # lower limit of hpd
N_2_Interval_hpd[2]  # upper limit of hpd

########### Dependence parameters ################

beta1_2
beta2_2
beta3_2
beta_ind_2

##############################

p_1_2_estimate
p_2_2_estimate
p_3_2_estimate


################### Capture probability of three sources  ##################################


est_p_1.._2
est_p_.1._2
est_p_..1_2


########### Results of Mbabane ###########################

#########  KP Sizes ###########

N_3_repli_mean_estimates
N_3_repli_median_estimates
N_3_se

N_3_Interval_hpd[1]  # lower limit of hpd
N_3_Interval_hpd[2]  # upper limit of hpd

########### Dependence parameters ################

beta1_3
beta2_3
beta3_3
beta_ind_3


################### Capture probability of three sources  ##################################

p_1_3_estimate
p_2_3_estimate
p_3_3_estimate

est_p_1.._3
est_p_.1._3
est_p_..1_3



############### Density plot KP size ############

#par(mfrow=c(2,2))
plot(density(M_2, bw=15),xlab="Population Size",ylab="Density",col="blue",lwd = 2, main = "")
plot(density(M_3, bw=20),xlab="Population Size",ylab="Density",col="blue",lwd = 2, main = "")




################### Violin plot of capture probabilities of three sources  ##################################


data_Manzini=data.frame(
  Sources=c(rep("uid",length(t)),rep("coupon",length(t)),rep("survey",length(t))),
  Probability=c(p_1.._2,p_.1._2,p_..1_2)
)


par(mfrow=c(1,1))
ggplot(data_Manzini, aes(x=Sources, y=Probability, fill=Sources)) + geom_violin(show.legend = FALSE)+xlab("")+ ylab("Capture probability")+theme(axis.title = element_text(size = 20),axis.text = element_text(size = 20)) 


data_Mbabane=data.frame(
  Sources=c(rep("uid",length(t)),rep("coupon",length(t)),rep("survey",length(t))),
  Probability=c(p_1.._3,p_.1._3,p_..1_3)
)


par(mfrow=c(1,1))
ggplot(data_Mbabane, aes(x=Sources, y=Probability, fill=Sources)) + geom_violin(show.legend = FALSE)+xlab("")+ ylab("Capture probability")+theme(axis.title = element_text(size = 20),axis.text = element_text(size = 20)) 



##############  Plot to test significance of dependence parameters  ################################

w=seq(0.001,0.1,0.005)
#par(mfrow=c(1,3))

Prob1_2=function(v) {
  x1=beta_1_2[t]
  return(length(x1[x1>v])/length(x1))
}

Prob1_2=Vectorize(Prob1_2)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2

plot(w, Prob1_2(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[1]>lambda ~ "|" ~  data)), cex.lab=clab,cex.axis=caxis)


Prob2_2=function(v) {
  x2=beta_2_2[t]
  return(length(x2[x2>v])/length(x2))
}

Prob2_2=Vectorize(Prob2_2)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2

plot(w, Prob2_2(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[2]>lambda ~ "|" ~  data)), cex.lab=clab,cex.axis=caxis)


Prob3_2=function(v) {
  x3=beta_3_2[t]
  return(length(x3[x3>v])/length(x3))
}

Prob3_2=Vectorize(Prob3_2)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2

plot(w, Prob3_2(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[3]>lambda ~ "|" ~  data)) ,cex.lab=clab,cex.axis=caxis)





Prob1_3=function(v) {
  x1=beta_1_3[t]
  return(length(x1[x1>v])/length(x1))
}

Prob1_3=Vectorize(Prob1_3)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2

plot(w, Prob1_3(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[1]>lambda ~ "|" ~  data)), cex.lab=clab,cex.axis=caxis)


Prob2_3=function(v) {
  x2=beta_2_3[t]
  return(length(x2[x2>v])/length(x2))
}

Prob2_3=Vectorize(Prob2_3)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2 

plot(w, Prob2_3(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[2]>lambda ~ "|" ~  data)), cex.lab=clab,cex.axis=caxis)


Prob3_3=function(v) {
  x3=beta_3_3[t]
  return(length(x3[x3>v])/length(x3))
}

Prob3_3=Vectorize(Prob3_3)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2

plot(w, Prob3_3(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[3]>lambda ~ "|" ~  data)) ,cex.lab=clab,cex.axis=caxis)




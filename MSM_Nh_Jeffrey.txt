
##### MSM in Nhalangano with Jeffrey's Prior for N  #####################

rm(list=ls())


################### R packages ##############

library(MCMCpack)
library(HDInterval)


tot=500000  # No of posteriror samples
th=20       # Thining   
lw=tot*0.5  # Burnin
up=tot

##########  Intial value specification ###############

N_estimate=array(NA,tot)

beta01=0.5
beta02=0.3
beta03=0.1 


p_1=0.2
p_2=0.6
p_3=0.3

a1=0.5
b1=0.5

a2=0.5
b2=0.5

a3=0.5
b3=0.5

beta_prior=array(0.5,4) # beta parameters are 0.5 refer Jeffrey's Dirichlet prior


##########################  Data ################################

# TRS (NHALANGANO) #
x_111=2
x_101=41
x_011=4
x_001=23
x_.1.=12
x_1..=106
x_..1=x_111+x_101+x_011+x_001


x_1.0=x_1..-(x_111+x_101)
x_.10=x_.1.-(x_111+x_011)
x0_star=x_111+x_101+x_011+x_001

x0_star=x_111+x_101+x_011+x_001



########## Gibbs ######################################

str=function(x) {
  x*log(x)-x 
}

######################################



beta_1=array(NA,tot)
beta_2= array(NA,tot)
beta_3= array(NA,tot)
beta_0= array(NA,tot)


p=array(NA,c(3,tot))

p_111=array(NA,tot)
p_110=array(NA,tot)
p_101=array(NA,tot)
p_011=array(NA,tot)
p_010=array(NA,tot)
p_001=array(NA,tot)
p_000=array(NA,tot)
p_100=array(NA,tot)

y111=array(NA,tot)

y110_1=array(NA,tot)
y110_2=array(NA,tot)
y110_3=array(NA,tot)

y011_1=array(NA,tot)
y011_2=array(NA,tot)
y011_3=array(NA,tot)

y100_1=array(NA,tot)
y100_2=array(NA,tot)
y100_3=array(NA,tot)

y101_1=array(NA,tot)
y101_2=array(NA,tot)
y101_3=array(NA,tot)

y010_1=array(NA,tot)
y010_2=array(NA,tot)
y010_3=array(NA,tot)

y001_1=array(NA,tot)
y001_2=array(NA,tot)
y001_3=array(NA,tot)

y000=array(NA,tot)

######################################

# initial values

beta_1[1]=beta01
beta_2[1]=beta02
beta_3[1]=beta03  

beta_0[1]=min((beta_1[1]+beta_2[1]+beta_3[1]),1)

p[1,1]=p_1
p[2,1]=p_2
p[3,1]=p_3


p_000[1]=(1-beta_0[1])*(1-p[1,1])*(1-p[2,1])*(1-p[3,1])

lower_trunc_limit_z_110=0 
upper_trunc_limit_z_110=min(x_.10,x_1.0)

z_110=round((lower_trunc_limit_z_110+upper_trunc_limit_z_110)/2)

z_100=x_1.0-z_110

z_010=x_.10-z_110

NN=x0_star + z_100 + z_010 + z_110

w=rnbinom(1,NN,(1-(1-beta_0[1])*(1-p[1,1])*(1-p[2,1])*(1-p[3,1])))
N_estimate[1]=w+NN

# for loop starts for Gibbs sampling

for(h in 2:(tot+1)){
  
  p_111[h-1]=((1-beta_0[h-1])*p[1,h-1]*p[2,h-1]*p[3,h-1])
  
  p_110[h-1]=((1-beta_0[h-1])*p[1,h-1]*p[2,h-1]*(1-p[3,h-1]))+(beta_2[h-1]*p[1,h-1]*p[2,h-1])+(beta_3[h-1]*p[1,h-1]*p[2,h-1])
  Q_110_1=((1-beta_0[h-1])*p[1,h-1]*p[2,h-1]*(1-p[3,h-1]))/p_110[h-1]
  Q_110_2=(beta_2[h-1]*p[1,h-1]*p[2,h-1])/p_110[h-1]
  Q_110_3=1-Q_110_1-Q_110_2
  prob_y110=c(Q_110_1,Q_110_2,Q_110_3)
  
  
  p_011[h-1]=((1-beta_0[h-1])*(1-p[1,h-1])*p[2,h-1]*p[3,h-1])+(beta_1[h-1]*(1-p[1,h-1])*p[3,h-1])+(beta_3[h-1]*(1-p[1,h-1])*p[2,h-1])
  Q_011_1=((1-beta_0[h-1])*(1-p[1,h-1])*p[2,h-1]*p[3,h-1])/p_011[h-1]
  Q_011_2=(beta_1[h-1]*(1-p[1,h-1])*p[3,h-1])/p_011[h-1]
  Q_011_3=1-Q_011_1-Q_011_2
  prob_y011=c(Q_011_1,Q_011_2,Q_011_3)
  
  
  p_100[h-1]=((1-beta_0[h-1])*p[1,h-1]*(1-p[2,h-1])*(1-p[3,h-1]))+(beta_1[h-1]*p[1,h-1]*(1-p[3,h-1]))+(beta_3[h-1]*p[1,h-1]*(1-p[2,h-1]))
  Q_100_1=((1-beta_0[h-1])*p[1,h-1]*(1-p[2,h-1])*(1-p[3,h-1]))/p_100[h-1]
  Q_100_2=(beta_1[h-1]*p[1,h-1]*(1-p[3,h-1]))/p_100[h-1]
  Q_100_3=1-Q_100_1-Q_100_2
  prob_y100=c(Q_100_1,Q_100_2,Q_100_3)
  
  
  p_001[h-1]=((1-beta_0[h-1])*(1-p[1,h-1])*(1-p[2,h-1])*p[3,h-1])+(beta_2[h-1]*(1-p[1,h-1])*(1-p[2,h-1]))+(beta_3[h-1]*(1-p[1,h-1])*(1-p[2,h-1]))
  Q_001_1=((1-beta_0[h-1])*(1-p[1,h-1])*(1-p[2,h-1])*p[3,h-1])/p_001[h-1]
  Q_001_2=(beta_2[h-1]*(1-p[1,h-1])*(1-p[2,h-1]))/p_001[h-1]
  Q_001_3=1-Q_001_1-Q_001_2
  prob_y001=c(Q_001_1,Q_001_2,Q_001_3)
  
  
  p_101[h-1]=((1-beta_0[h-1])*p[1,h-1]*(1-p[2,h-1])*p[3,h-1])+(beta_1[h-1]*p[1,h-1]*p[3,h-1])+(beta_2[h-1]*p[1,h-1]*(1-p[2,h-1]))
  Q_101_1=(1-beta_0[h-1])*p[1,h-1]*(1-p[2,h-1])*p[3,h-1]/p_101[h-1]
  Q_101_2=(beta_1[h-1]*p[1,h-1]*p[3,h-1])/p_101[h-1]
  Q_101_3=1-Q_101_1-Q_101_2
  prob_y101=c(Q_101_1,Q_101_2,Q_101_3)
  
  
  p_010[h-1]=((1-beta_0[h-1])*(1-p[1,h-1])*p[2,h-1]*(1-p[3,h-1]))+(beta_1[h-1]*(1-p[1,h-1])*(1-p[3,h-1]))+(beta_2[h-1]*(1-p[1,h-1])*p[2,h-1])
  Q_010_1=((1-beta_0[h-1])*(1-p[1,h-1])*p[2,h-1]*(1-p[3,h-1]))/p_010[h-1]
  Q_010_2=(beta_1[h-1]*(1-p[1,h-1])*(1-p[3,h-1]))/p_010[h-1]
  Q_010_3=1-Q_010_1-Q_010_2
  prob_y010=c(Q_010_1,Q_010_2,Q_010_3)
  
  
  p_000[h-1]=(1-beta_0[h-1])*(1-p[1,h-1])*(1-p[2,h-1])*(1-p[3,h-1])
  
  ######  Generation of Z_110 
  
  P110_str=p_110[h-1]/(p_110[h-1]+p_100[h-1]+p_010[h-1]+p_000[h-1])
  P100_str=p_100[h-1]/(p_110[h-1]+p_100[h-1]+p_010[h-1]+p_000[h-1])
  P010_str=p_010[h-1]/(p_110[h-1]+p_100[h-1]+p_010[h-1]+p_000[h-1])
  P000_str=p_000[h-1]/(p_110[h-1]+p_100[h-1]+p_010[h-1]+p_000[h-1])
  
  pr_1=P110_str/(P110_str+P100_str)
  
  pr_2=P010_str/(P010_str+P000_str)
  
  lower_z110=max(0,x0_star+x_1.0+x_.10-N_estimate[h-1])
  upper_z110=min(x_1.0,x_.10)
  
  prob1=c()
  
  for(u in lower_z110:upper_z110){
    
    prob1=c(prob1,(dbinom(u,x_1.0,pr_1)*dbinom((x_.10-u),(N_estimate[h-1]-x0_star-x_1.0),pr_2)))   
    
  }
  
  
  prob2=prob1/sum(prob1)
  
  if (lower_z110==upper_z110) {
    z_110=upper_z110 
  } else {
    z_110=sample(lower_z110:upper_z110,1,replace =TRUE,prob2)
  }
  
  z_100=x_1.0-z_110
  
  z_010=x_.10-z_110
  
  NN=x0_star + z_100 + z_010 + z_110
  
  ###############################################
  
  y110_vec_draw=rmultinom(1,z_110,prob_y110)
  
  y110_1[h]=y110_vec_draw[1,1]
  y110_2[h]=y110_vec_draw[2,1]          
  y110_3[h]=z_110-y110_1[h]-y110_2[h]
  
  
  
  y011_vec_draw=rmultinom(1,x_011,prob_y011)
  
  y011_1[h]=y011_vec_draw[1,1]
  y011_2[h]=y011_vec_draw[2,1]
  y011_3[h]=x_011-y011_1[h]-y011_2[h]
  
  
  
  y100_vec_draw=rmultinom(1,z_100,prob_y100)
  
  y100_1[h]=y100_vec_draw[1,1]
  y100_2[h]=y100_vec_draw[2,1]
  y100_3[h]=z_100-y100_1[h]-y100_2[h]
  
  
  y101_vec_draw=rmultinom(1,x_101,prob_y101)
  
  y101_1[h]=y101_vec_draw[1,1]
  y101_2[h]=y101_vec_draw[2,1]
  y101_3[h]=x_101-y101_1[h]-y101_2[h]
  
  y010_vec_draw=rmultinom(1,z_010,prob_y010)
  
  y010_1[h]=y010_vec_draw[1,1]
  y010_2[h]=y010_vec_draw[2,1]
  y010_3[h]=z_010-y010_1[h]-y010_2[h]
  
  
  y001_vec_draw=rmultinom(1,x_001,prob_y001)
  
  y001_1[h]=y001_vec_draw[1,1]
  y001_2[h]=y001_vec_draw[2,1]
  y001_3[h]=x_001-y001_1[h]-y001_2[h]
  
  d1=y011_2[h]+y100_2[h]+y101_2[h]+y010_2[h]+beta_prior[1]
  d2=y110_2[h]+x_101-y101_1[h]-y101_2[h]+z_010-y010_1[h]-y010_2[h]+y001_2[h]+beta_prior[2]
  d3=z_110-y110_1[h]-y110_2[h]+x_011-y011_1[h]-y011_2[h]+z_100-y100_1[h]-y100_2[h]+x_001-y001_1[h]-y001_2[h]+beta_prior[3]
  d4=x_111+(N_estimate[h-1]-NN)+(y110_1[h]+y011_1[h]+y100_1[h]+y101_1[h]+y010_1[h]+y001_1[h])+beta_prior[4]
  
  beta_vec_draw=rdirichlet(1,c(d1,d2,d3,d4))
  
  beta_1[h]=beta_vec_draw[1,1]
  beta_2[h]=beta_vec_draw[1,2]
  beta_3[h]=beta_vec_draw[1,3]  
  beta_0[h]=min((beta_1[h]+beta_2[h]+beta_3[h]),1)
  
  
  
  a_p1=x_111+z_110+z_100+x_101+a1
  b_p1=x_011+z_010+x_001+N_estimate[h-1]-NN+b1
  
  a_p2=x_111+z_110+x_011-y011_2[h]+z_010-y010_2[h]+a2
  b_p2=z_100-y100_2[h]+x_101-y101_2[h]+x_001+N_estimate[h-1]-NN+b2
  
  a_p3=x_111+y011_1[h]+y011_2[h]+y101_1[h]+y101_2[h]+y001_1[h]+a3
  b_p3=y110_1[h]+y100_1[h]+y100_2[h]+y010_1[h]+y010_2[h]+N_estimate[h-1]-NN+b3
  
  p[1,h]=rbeta(1,a_p1,b_p1)
  p[2,h]=rbeta(1,a_p2,b_p2)
  p[3,h]=rbeta(1,a_p3,b_p3)
  
  w=rnbinom(1,NN,(1-(1-beta_0[h])*(1-p[1,h])*(1-p[2,h])*(1-p[3,h])))
  N_estimate[h]=w+NN
  
} # end of the 'h' loop






###################################################3

t=seq(lw,up, th)

par(mgp = c(2, 0.5, 0))
plot(t, N_estimate[seq(lw, up, th)], xlab = "Iteration", ylab = "Population Size", type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p[1, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[1])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p[2, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[2])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, p[3, seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(p[3])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_1[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[1])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_2[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[2])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')
plot(t, beta_3[seq(lw, up, th)], xlab = "Iteration", ylab = expression(bold(beta[3])), type = "l", cex.lab = 2, font.lab = 2, xaxt = 'n', yaxt = 'n')

M=N_estimate[seq(lw,up, th)]

N_mean_estimate=round(mean(M))
N_median_estimate=round(median(M))
N_se=sd(M)




beta1_estimate=median(beta_1[seq(lw,up, th)])
beta2_estimate=median(beta_2[seq(lw,up, th)])
beta3_estimate=median(beta_3[seq(lw,up, th)])
beta_ind_estimate=median(beta_0[seq(lw,up, th)])


p_1_estimate=median(p[1,seq(lw,up, th)])
p_2_estimate=median(p[2,seq(lw,up, th)])
p_3_estimate=median(p[3,seq(lw,up, th)])

Interval_hpd=hdi(M, credMass=0.95) #  95% highest posterior credible interval

p_1..=p_111[seq(lw,up, th)]+p_110[seq(lw,up, th)]+p_101[seq(lw,up, th)]+p_100[seq(lw,up, th)]
p_.1.=p_111[seq(lw,up, th)]+p_110[seq(lw,up, th)]+p_011[seq(lw,up, th)]+p_010[seq(lw,up, th)]
p_..1=p_111[seq(lw,up, th)]+p_101[seq(lw,up, th)]+p_011[seq(lw,up, th)]+p_001[seq(lw,up, th)]

est_p_1..=median(p_1..)
est_p_.1.=median(p_.1.)
est_p_..1=median(p_..1)

##########################################



#######################  Estimates ####################################


#########  KP Sise ###########

N_mean_estimate
N_median_estimate
N_se

Interval_hpd

########### Dependence parameters ################

beta1_estimate
beta2_estimate
beta3_estimate
beta_ind_estimate


############# 

p_1_estimate
p_2_estimate
p_3_estimate



est_p_1..  # List-1 capture probability
est_p_.1.  # List-2 capture probability
est_p_..1  # List-3 capture probability



################### Capture probability of three sources  ##################################



data_Nhalangano=data.frame( uid = p_1.., rainbow = p_.1., survey= p_..1)

 
data_Nhalangano=data.frame(
  Sources=c(rep("uid",length(t)),rep("rainbow",length(t)),rep("survey",length(t))),
  Probability=c(p_1..,p_.1.,p_..1)
)


par(mfrow=c(1,1))
ggplot(data_Nhalangano, aes(x=Sources, y=Probability, fill=Sources)) + geom_violin(show.legend = FALSE)+xlab("")+ ylab("Capture probability")+theme(axis.title = element_text(size = 20),axis.text = element_text(size = 20)) 


############### Density plot KP size ############

par(mfcol=c(1,1))
plot(density(M),xlab="Population Size",ylab="Density",col="blue",lwd = 2, main = "")



##############  Plot to test significance of dependence parameters  ################################

w=seq(0.001,0.1,0.005)
par(mfrow=c(1,1))

Prob1_3=function(v) {
  x1=beta_1[t]
  return(length(x1[x1>v])/length(x1))
}

Prob1_3=Vectorize(Prob1_3)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2

plot(w, Prob1_3(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[1]>lambda ~ "|" ~  data)), cex.lab=clab,cex.axis=caxis)



Prob2_3=function(v) {
  x2=beta_2[t]
  return(length(x2[x2>v])/length(x2))
}

Prob2_3=Vectorize(Prob2_3)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2

plot(w, Prob2_3(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[2]>lambda ~ "|" ~  data)), cex.lab=clab,cex.axis=caxis)


Prob3_3=function(v) {
  x3=beta_3[t]
  return(length(x3[x3>v])/length(x3))
}

Prob3_3=Vectorize(Prob3_3)

par(mar = c(5, 7, 4, 2) + 0.1)
clab = 2
caxis = 2

plot(w, Prob3_3(w), type="h", col="blue",lwd=2, xlim=c(0,0.05), ylim=c(0,1), xlab=expression(~ lambda), ylab=expression(P(beta[3]>lambda ~ "|" ~  data)) ,cex.lab=clab,cex.axis=caxis)







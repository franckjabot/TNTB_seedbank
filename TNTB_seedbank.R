## function to draw from a log-series distribution, 
## that is truncated at Spool species.
lssim<-function(thet,Spool){
  u<-runif(1)
  k<-1
  P<- -(1-thet)/log(thet)
  F<-P
  while ((u>=F)&(k<Spool)){
    P<-P*k*(1-thet)/(k+1)
    k<-k+1
    F<-F+P
  }
  if (k==Spool){
    k=sample(1:Spool,1)
  }
  k
}

draw_ind_logseries<-function(tet,Spool){
  thet=1-exp(-1/tet)
  res=lssim(thet,Spool)
  res
}

## function to draw n seeds from a regional pool
draw_regional_pool<-function(n,theta,Spool){
  res=NULL
  if (n>0){
    res=array(0,Spool)
    for (k in 1:n){
      spec=draw_ind_logseries(theta,Spool)
      res[spec]=res[spec]+1
    }
  }
  res
}

## function to simulate the seeding of the adult community
seed_output2<-function(seed_bank,community,p,fitness){
 res=seed_bank+community*p*fitness
 res
}

## function to simulate adult mortality
adult_mortality<-function(community,d_plant,Spool,J){
 res=community
 rescum=cumsum(community)
 Jdead=rbinom(1,J,d_plant)
 if (Jdead>0){
  pos=sample(1:J,Jdead)
  for (k in 1:Jdead){
   pos_k=min((1:Spool)[rescum>=(pos[k])])
   res[pos_k]=res[pos_k]-1
  }
 }
list("community"=res,"J_dead"=Jdead)
}

## recruitment from the seed bank
recruitment2<-function(seed_bank,J){
 as.numeric(t(rmultinom(1,J,seed_bank)))
}

## function to simulate one step dynamics
step_dynamics2<-function(community,seed_bank,p,d_seed,d_plant,i,theta,Spool,A){
  J=sum(community)
  
  # seeding
  fitness=rlnorm(n=Spool,meanlog=-0.5*(log(A+1)),sdlog=sqrt(log(A+1)))
  new_seed_bank=seed_output2(seed_bank,community,p,fitness)
  
  # immigration from the regional pool
  new_seed_bank=new_seed_bank+draw_regional_pool(i,theta,Spool)

  # seed mortality
  new_seed_bank=new_seed_bank*(1-d_seed)
  
  # adult mortality
  temp=adult_mortality(community,d_plant,Spool,J)
  surviving_community=temp$community
  
  new_community=array(0,Spool)
  if (temp$J_dead>0){
   # adult recruitment
   new_community=recruitment2(new_seed_bank,temp$J_dead)
   new_seed_bank=new_seed_bank-new_community
   for (k in 1:Spool){
    if (new_seed_bank[k]<0){
     new_seed_bank[k]=0
    }
   }
  }
  list("community"=surviving_community+new_community,"seed_bank"=new_seed_bank)
}

## computation of summary statistics
evenness<-function(community,J){
  com=community[community>0]
  res=0
  if (length(com)>1){
    res=-sum((com/J)[(com/J)>0]*log((com/J)[(com/J)>0]))/log(length(com))
  }
res
}

sorensen<-function(table_community1,table_community2,Stot){
  temp=table_community1*table_community2
  C=length(temp[temp>0])
  S1=length(table_community1[table_community1>0])
  S2=length(table_community2[table_community2>0])
  res=2*C/(S1+S2)
res
}

bray_curtis<-function(table_community1,table_community2,Stot){
  num=0
  denomtot=0
  for (i in 1:Stot){
    denom=table_community1[i]+table_community2[i]
    if (denom>0){
      num=num+(2*min(table_community1[i],table_community2[i]))
      denomtot=denomtot+denom
    }
  }
  res=0
  if (denomtot>0){
   res=(num/denomtot)
  }
res
}

make_com<-function(table_com){
  res=array(0,sum(table_com))
  l=length(table_com)
  ind=0
  for (i in 1:l){
   k=table_com[i]
   if (k>0){
    for (j in 1:k){
     res[(ind+j)]=i
    }
    ind=ind+k
   }
  }
res
}

species_individual_relationship<-function(com, table_N){
  l=length(table_N)
  res=array(0,l)
  for (i in 1:l){
   for (j in 1:10){
    samp=sample(com,table_N[i])
    res[i]=res[i]+(length(table(samp)))/3
   }
  }
  as.numeric(summary(lm(log(res)~log(table_N)))$coefficients[,1])
}

## function that simulates a multi-year dynamics 
## and output summary statistics
stochastic_dynamics2<-function(param){
  p=param[1]
  d_seed=param[2]
  d_plant=param[3]
  m=param[4]
  theta=param[5]
  Spool=param[6]
  A=param[7]
  nb_steps=param[8]

  community=com_ini
  seed_bank=seed_bank_ini

  N_bank=array(0,nb_steps)
  S_bank=array(0,nb_steps)
  S_com=array(0,nb_steps)
  H_bank=array(0,nb_steps)
  H_com=array(0,nb_steps)
  Sor=array(0,nb_steps)
  BC=array(0,nb_steps)
  BC_temp=array(0,nb_steps)
  BC_bank_temp=array(0,nb_steps)
  Sor_samp=array(0,nb_steps)
  for (i in 1:nb_steps){
    temp=step_dynamics2(community,seed_bank,p,d_seed,d_plant,m,theta,Spool,A)
    if (i>1){
      BC_bank_temp[i]=bray_curtis(seed_bank/sum(seed_bank),temp$seed_bank/sum(temp$seed_bank),Spool)    
    }
    seed_bank=temp$seed_bank
    N_bank[i]=sum(seed_bank)
    S_bank[i]=length(seed_bank[seed_bank>0])
    H_bank[i]=evenness(seed_bank,N_bank[i])
    BC_temp[i]=bray_curtis(community/J,temp$community/J,Spool)    
    community=temp$community
    S_com[i]=length(community[community>0])
    H_com[i]=evenness(community,J)
    Sor[i]=sorensen(community,seed_bank,Spool)
    seed_bank_samp=as.numeric(t(rmultinom(1,J,seed_bank)))
    Sor_samp[i]=sorensen(community,seed_bank_samp,Spool)
    BC[i]=bray_curtis(community/J,seed_bank/N_bank[i],Spool)
  }
list("community"=community,"seed_bank"=seed_bank,"N_bank"=N_bank,"S_bank"=S_bank,"S_com"=S_com,"H_bank"=H_bank,"H_com"=H_com,"Sor"=Sor,"BC"=BC,"Sor_samp"=Sor_samp,"BC_temp"=BC_temp,"BC_bank_temp"=BC_bank_temp)
}

## function that simulates a multi-year dynamics 
## and output summary statistics
stochastic_dynamics_output2<-function(param){
  p=param[1]
  d_seed=param[2]
  d_plant=param[3]
  m=param[4]
  theta=param[5]
  Spool=param[6]
  A=param[7]
  nb_steps=param[8]

  community=com_ini
  seed_bank=seed_bank_ini

  N_bank=array(0,nb_steps)
  S_bank=array(0,nb_steps)
  S_com=array(0,nb_steps)
  H_bank=array(0,nb_steps)
  H_com=array(0,nb_steps)
  Sor=array(0,nb_steps)
  BC=array(0,nb_steps)
  Sor_samp=array(0,nb_steps)
  Sor_samp2=array(0,nb_steps)
  Sor_samp3=array(0,nb_steps)
  BC_samp=array(0,nb_steps)
  BC_samp2=array(0,nb_steps)
  BC_samp3=array(0,nb_steps)
  S_bank_samp=array(0,nb_steps)
  S_bank_samp2=array(0,nb_steps)
  S_bank_samp3=array(0,nb_steps)
  H_bank_samp=array(0,nb_steps)
  H_bank_samp2=array(0,nb_steps)
  H_bank_samp3=array(0,nb_steps)
  for (i in 1:3000){
    temp=step_dynamics2(community,seed_bank,p,d_seed,d_plant,m,theta,Spool,A)
    seed_bank=temp$seed_bank
    community=temp$community
  }
  for (i in 3001:nb_steps){
   temp=step_dynamics2(community,seed_bank,p,d_seed,d_plant,m,theta,Spool,A)
   seed_bank=temp$seed_bank
   if (i%%10==0){
    N_bank[i]=sum(seed_bank)
    S_bank[i]=length(seed_bank[seed_bank>0])
    H_bank[i]=evenness(seed_bank,N_bank[i])
   }
   community=temp$community
   if (i%%10==0){
    S_com[i]=length(community[community>0])
    H_com[i]=evenness(community,J)
    Sor[i]=sorensen(community,seed_bank,Spool)
    seed_bank_samp=as.numeric(t(rmultinom(1,ceiling(N_bank[i]*0.0001),seed_bank)))
    seed_bank_samp2=as.numeric(t(rmultinom(1,ceiling(N_bank[i]*0.001),seed_bank)))
    seed_bank_samp3=as.numeric(t(rmultinom(1,ceiling(N_bank[i]*0.01),seed_bank)))
    Sor_samp[i]=sorensen(community,seed_bank_samp,Spool)
    Sor_samp2[i]=sorensen(community,seed_bank_samp2,Spool)
    Sor_samp3[i]=sorensen(community,seed_bank_samp3,Spool)
    S_bank_samp[i]=length(seed_bank_samp[seed_bank_samp>0])
    S_bank_samp2[i]=length(seed_bank_samp2[seed_bank_samp2>0])
    S_bank_samp3[i]=length(seed_bank_samp3[seed_bank_samp3>0])
    H_bank_samp[i]=evenness(seed_bank_samp,ceiling(N_bank[i]*0.0001))
    H_bank_samp2[i]=evenness(seed_bank_samp2,ceiling(N_bank[i]*0.001))
    H_bank_samp3[i]=evenness(seed_bank_samp3,ceiling(N_bank[i]*0.01))
    BC_samp[i]=bray_curtis(community/J,seed_bank_samp/ceiling(N_bank[i]*0.0001),Spool)
    BC_samp2[i]=bray_curtis(community/J,seed_bank_samp2/ceiling(N_bank[i]*0.001),Spool)
    BC_samp3[i]=bray_curtis(community/J,seed_bank_samp3/ceiling(N_bank[i]*0.01),Spool)
    BC[i]=bray_curtis(community/J,seed_bank/N_bank[i],Spool)
   }
  }
  samp=3000+10*(1:200)

c(mean(N_bank[samp]),mean(S_bank[samp]),mean(H_bank[samp]),mean(S_com[samp]),mean(H_com[samp]),mean(Sor[samp]),mean(BC[samp]),mean(Sor_samp[samp]),mean(Sor_samp2[samp]),mean(Sor_samp3[samp]),mean(BC_samp[samp]),mean(BC_samp2[samp]),mean(BC_samp3[samp]),mean(S_bank_samp[samp]),mean(S_bank_samp2[samp]),mean(S_bank_samp3[samp]),mean(H_bank_samp[samp]),mean(H_bank_samp2[samp]),mean(H_bank_samp3[samp]),sd(N_bank[samp]),sd(S_bank[samp]),sd(H_bank[samp]),sd(S_com[samp]),sd(H_com[samp]),sd(Sor[samp]),sd(BC[samp]))
}

## function to initialize the community and seed bank
initialization<-function(param){
  com_ini<<-draw_regional_pool(J,param[5],param[6])
  seed_bank_ini<<-array(0,param[6])
}

names_output=c("N_bank","S_seed","H_seed","S_plant","H_plant","Sor","BC","Sor_samp","Sor_samp2","Sor_samp3","BC_samp","BC_samp2","BC_samp3","S_bank_samp","S_bank_samp2","S_bank_samp3","H_bank_samp","H_bank_samp2","H_bank_samp3","sd_N_bank","sd_S_seed","sd_H_seed","sd_S_plant","sd_H_plant","sd_Sor","sd_BC")

J=100
com_ini=NULL
seed_bank_ini=NULL
nb_steps=6000
samp=3000+10*(1:200)

param_test=c(10,0.1,0.1,10,25,50,0,5000)
initialization(param_test)
results_test=stochastic_dynamics_output2(param_test)
names(results_test)=names_output
results_test





rm(list = ls())

library(knitr)
library(abind)

# source("C:\\UCLA\\thesis_ideas\\PhD_thesis\\many-to-many-matching\\large_matching_simulations.R")
source("C:\\UCLA\\thesis_ideas\\PhD_thesis\\many-to-many-matching\\rpm.model.matrix.R")
source("C:\\UCLA\\thesis_ideas\\PhD_thesis\\many-to-many-matching\\Gale_Shapley.R")
source("C:\\UCLA\\thesis_ideas\\PhD_thesis\\many-to-many-matching\\Gale_Shapley_many.R")

# zz=file("C:\\UCLA\\thesis_ideas\\PhD_thesis\\MenzelMatlab-20170816T001358Z-001\\MenzelMatlab\\eta.bin", "rb")
# test_eta = readBin(zz,double(),40000)
# test_eta = matrix(test_eta, nrow=200)
# close(zz)
# 
# zz=file("C:\\UCLA\\thesis_ideas\\PhD_thesis\\MenzelMatlab-20170816T001358Z-001\\MenzelMatlab\\zeta.bin", "rb")
# test_zeta = readBin(zz,double(),40000)
# test_zeta = matrix(test_zeta, nrow=200)
# close(zz)

set.seed(12345)

reme = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

minVal = function(x) {
  temp=x[x>0]
  if (length(temp)==0) {
    z=0
  } else {
    z = min(temp)
  }
  z
}

################## create simulated data ##########################

# specify the formula for utilities
ff = ~ b1cov("f1") + b2cov("f1") + b1absdiff("f1",1) + b2absdiff("f1",1)
# parse the formula
model.terms <- rownames(attr(terms.formula(ff), "factors"))
temp <- strsplit(model.terms, "[(]")
model.terms.names <- unlist(lapply(temp, `[[`, 1))
temp2 <- gsub(")", "", gsub("\"", "", unlist(temp)))
temp2 <- gsub("[ \t\n\r\f\v]", "", temp2)
model.terms.coef.names <- temp2[-which(temp2 %in% model.terms.names)]
model.terms.coef.names <- strsplit(model.terms.coef.names, ",")

# symmetric beta for b1 and b2
symmetric = FALSE
# true parameter for utility generation
alpha = 2
theta_w = c(alpha,0,0)
theta_m = c(alpha,0,0)


n_U_choices = 1 # FCY: number of choice situations (assume same for both sides)
n_V_choices = 6

Gamma = 0.5*(sqrt(1+4*exp(2*alpha))-1) # FCY: solving Gamma using the quadratic formula
prob_gam = 1/(1+Gamma) # FCY: probability of remaining single
# prob_gam = (prob_gam * n_choices)^ n_choices # FCY wrong
# prob_gam = 1/(1+Gamma*n_choices) # FCY wrong

##################################################################
i = 1
# N = 100
N = c(10, 20, 50, 100, 200) #, 1000, 2000, 5000) # in the code below, nw=nm=n
B = 50

variance = matrix(NA, nrow=length(N), ncol=1)
spread = matrix(NA, nrow=length(N), ncol=1)
# corr1 = matrix(NA, nrow=length(N), ncol=1)
# corr2 = matrix(NA, nrow=length(N), ncol=1)
prob_diff = matrix(NA, nrow=length(N), ncol=1)
prob_sin = matrix(NA, nrow=length(N), ncol=1)
freq_sin = matrix(NA, nrow=length(N), ncol=1)
freq_sd = matrix(NA, nrow=length(N), ncol=1)

wc_av1 = matrix(NA, nrow=length(N), ncol=1)
wc_av2 = matrix(NA, nrow=length(N), ncol=1)
I_star = matrix(NA, nrow=length(N), ncol=1)
I_circ = matrix(NA, nrow=length(N), ncol=1)
ch_num = matrix(NA, nrow=length(N), ncol=1)


for (n in N){
  
  print(paste0("N = ", n))
  
  m = matrix(B, ncol=1)
  s = matrix(B, ncol=1)
  v = matrix(B, ncol=1)
  n_m1 = matrix(B, ncol=1)
  n_m2 = matrix(B, ncol=1)
  r11 = matrix(B, ncol=1)
  r12 = matrix(B, ncol=1)
  r22 = matrix(B, ncol=1)
  s11 = matrix(B, ncol=1)
  s12 = matrix(B, ncol=1)
  s22 = matrix(B, ncol=1)
  
  f_sg = matrix(B, ncol=1)
  p_sg = matrix(B, ncol=1)
  df_single = matrix(B, ncol=1)
  n_ch = matrix(B, ncol=1)
  
  
  nw=n
  nm=n
  
  # create X, Z data
  X = matrix(rep(0:1,times=nw/2),ncol=1)
  colnames(X) = "f1"
  Z = matrix(rep(0:1,times=nm/2),ncol=1)
  colnames(Z) = "f1"
  
  I1w = matrix(NA, nrow=dim(X)[1], ncol=B)
  I2w = matrix(NA, nrow=dim(X)[1], ncol=B)
  I1m = matrix(NA, nrow=dim(Z)[1], ncol=B)
  I2m = matrix(NA, nrow=dim(Z)[1], ncol=B)
  
  
  # create utility matrices
  # U_star = matrix(0, nrow=nw, ncol = nm)
  U_star = matrix(0, nrow=dim(X)[1], ncol = dim(Z)[1])
  # V_star = matrix(0, nrow=nm, ncol = nw)
  V_star = matrix(0, nrow=dim(Z)[1], ncol = dim(X)[1])
  modelmat <- rpm.model.matrix(model.terms.names, model.terms.coef.names, X, Z)
  for (ii in 1:dim(modelmat$X)[3]) {
    U_star = U_star + modelmat$X[,,ii] * theta_w[ii]
  }
  for (ii in 1:dim(modelmat$Z)[3]) {
    V_star = V_star + modelmat$Z[,,ii] * theta_m[ii]
  }
  # Menzel's code did t(V_star)
  V_star = t(V_star)
  
  
  for (b in 1:B){
    
    print(paste0("b = ", b))
    
    eta  <- -log(-log(matrix(runif(dim(X)[1] * dim(Z)[1]), dim(X)[1])))
    zeta <- -log(-log(matrix(runif(dim(X)[1] * dim(Z)[1]), dim(X)[1])))
    
    
    
    # new formulation of outside option
    J = round(sqrt(n))
    eta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * dim(X)[1]), J))), 2, max),ncol=1),1,dim(Z)[1])
    zeta0 <- reme(matrix(apply(-log(-log(matrix(runif(J * dim(Z)[1]), J))),2,max),nrow=1), dim(X)[1], 1)
    
    
    eta = eta - eta0
    zeta = zeta - zeta0
    
    U = U_star + eta
    V = V_star + zeta
  
    
    # generate the matching (W-optimal)
    # uses Menzel's GS which allows remaining single
    # mu1 = Gale_Shapley(U,V)
    mu1 = Gale_Shapley_many(U,V,n_U_choices,n_V_choices)
    
    # U_bar = matrix(rowSums(mu1*U),ncol=1) # FCY: each element is the utility from the chosen man
    # V_bar = matrix(colSums(mu1*V),nrow=1) # FCY: each element is the utility from the chosen woman
    # U_bar = matrix(apply(mu1*U, 1, max),ncol=1) # FCY: modifed for many-to-one and many-to-many -- each element is the max utility among the chosen men
    # V_bar = matrix(apply(mu1*V, 2, max),nrow=1) # FCY: modifed for many-to-one and many-to-many -- each element is the max utility among the chosen women
    U_bar = matrix(apply(mu1*U, 1, minVal) ,ncol=1) # FCY: modifed for many-to-one and many-to-many -- each element is the min utility among the chosen men
    V_bar = matrix(apply(mu1*V, 2, minVal),nrow=1) # FCY: modifed for many-to-one and many-to-many -- each element is the min utility among the chosen women
    
    w_av1 = (U - reme(U_bar, 1, dim(Z)[1]))>=0 # FCY: (r,c)=1 means woman r is in man c's opportunity set (man c has utility at least as high as her husband)
    m_av1 = (V - reme(V_bar, dim(X)[1], 1))>=0
    I_w1 = rowSums(m_av1*exp(U_star)) # FCY: I_wl[i] is the sum of woman i's exp(U_star) of her opp. set
    I_m1 = colSums(w_av1*exp(V_star))
    n_m1[b] = mean(colSums(w_av1)) # FCY: the average # of women in a man's opportunity set
    I_m1 = t(I_m1)
    
    
    # choice probabilities
    J_hat = colSums(w_av1) # FCY: J_hat[i] is the # of women in man i's opportunity set
    # p_single = J/(J + I_m1)
    p_single = J/(J + J_hat*exp(alpha)) # FCY: multiply J/J because we are using J_hat? (this is not displayed in Menzel's paper)
    # p_single = (J/(J + J_hat*exp(alpha)))^n_choices # FCY: wrong choosing single option in every choice situation
    # p_single = (J/(J + J_hat*exp(alpha)/n_choices)) ^ n_choices # FCY: wrong: modified to accomodate many-to-one and many-to-many relationships
    # p_single = J/(J+n_m1(b,1)*exp(2));
    
    # f_single = 1 - colSums(mu1) # FCY: f_single[i] is 1 when man i is single
    f_single = 1*(colSums(mu1)==0) # FCY: changed to accomodate many-to-one and many-to-many relationships
    
    f_sg[b] = mean(f_single) # FCY: mean frequency of being single
    # f_sg[b] = f_sg[b] ^(1/n_choices)/n_choices # FCY: added to accomodate many-to-one and many-to-many relationships
    
    p_sg[b] = mean(p_single) # FCY: mean probability of being single?
    df_single[b] = mean(f_single-p_single) # FCY: not included in the table
    
    
    
    
    #mu2 = Gale_Shapley(t(V),t(U))
    mu2 = Gale_Shapley_many(t(V),t(U),n_U_choices,n_V_choices)
    mu2 = t(mu2)

    # U_bar = matrix(rowSums(mu2*U),ncol=1) # FCY: each element is the tility from the chosen man
    # V_bar = matrix(colSums(mu2*V),nrow=1) # FCY: each element is the utility from the chosen woman
    # U_bar = matrix(apply(mu2*U, 1, max),ncol=1) # FCY: modifed for many-to-one and many-to-many -- each element is the max utility among the chosen men
    # V_bar = matrix(apply(mu2*V, 2, max),nrow=1) # FCY: modifed for many-to-one and many-to-many -- each element is the max utility among the chosen women
    U_bar = matrix(apply(mu2*U, 1, minVal),ncol=1) # FCY: modifed for many-to-one and many-to-many -- each element is the min utility among the chosen men
    V_bar = matrix(apply(mu2*V, 2, minVal),nrow=1) # FCY: modifed for many-to-one and many-to-many -- each element is the min utility among the chosen women
    
    w_av2 = (U >= reme(U_bar,1,dim(Z)[1]))
    m_av2 = (V >= reme(V_bar,dim(X)[1],1))
    I_w2 = rowSums(m_av2*exp(U_star))
    I_m2 = colSums(w_av2*exp(V_star))
    n_m2[b] = mean(colSums(w_av2))
    # J1(b,1) = mean(sum(w_av1)<1)
    # J2(b,1) = mean(sum(w_av1)<2)
    # J5(b,1) = mean(sum(w_av1)<5)
    # J10(b,1) = mean(sum(w_av1)<10)
    I_m2 = t(I_m2)
    
    
    n_ch[b] = sum(colSums(w_av2-w_av1)>0) # FCY: check if there are more men in w_av2 than in w_av1...............
    m[b] = mean(I_m2-I_m1)
    v[b] = sd(I_m1)^2
    
    # Imat_w1 = kronecker(matrix(1, 1, n), I_w1)
    # R = cov(array(eta', c(%...)),reshape(eta',n^2,1))
    # r11(b,1) = R(1,1)
    # r12(b,1) = R(2,1)
    # r22(b,1) = R(2,2)
    # 
    # Imat_w2 = kronecker(matrix(1, 1, n), I_w2)
    # R = cov(array(eta', c(%...)),reshape(eta',n^2,1))
    # s11(b,1) = R(1,1)
    # s12(b,1) = R(2,1)
    # s22(b,1) = R(2,2)
    
    I1w[,b] = I_w1
    I2w[,b] = I_w2
    I1m[,b] = I_m1
    I2m[,b] = I_m2
  
  }
  
  variance[i] = mean(v) # FCY: not displayed in Menzel's table
  spread[i] = mean(m) # FCY: not displayed in Menzel's table
  # corr1[i] = mean(r12)./sqrt(mean(r11)*mean(r22))
  # corr2[i] = mean(s12)./sqrt(mean(s11)*mean(s22))
  
  prob_diff[i] = mean(df_single)
  prob_sin[i] = mean(p_sg)
  freq_sin[i] = mean(f_sg)
  freq_sd[i] = sd(f_sg)
  
  # FCY commented out because multdens is undefined
  #     x = 0:0[['05']]:5;
  #     I1m(1) = [];
  #     I2m(1) = [];
  #     f = multdens(I1m,I2m,x');
  #     dens_m(i,:) = f;
  #     I1w(1) = [];
  #     I2w(1) = [];
  #     f = multdens(I1w,I2w,x');
  #     dens_w(i,:) = f;
  
  
  wc_av1[i] = mean(n_m1) # FCY: average # of women in man i's opportunity set
  wc_av2[i] = mean(n_m2)
  
  # lt1(i,1) = mean(J1)
  # lt2(i,1) = mean(J2)
  # lt5(i,1) = mean(J5)
  # lt10(i,1) = mean(J10)
  
  I_star[i] = mean(I2m)
  I_circ[i] = mean(I1m)
  
  ch_num[i] = mean(n_ch)
  
  i = i+1
  
}

w_av_diff = wc_av2 - wc_av1
I_diff = I_star-I_circ
table_1 = data.frame(n=N, diff_size_opp_set_m = w_av_diff, diff_opp_set = ch_num, 
                     ave_incl_val_m_m_prop = I_star/sqrt(N), ave_incl_val_m_w_prop = I_circ/sqrt(N),
                     diff_ave_incl_val_m = I_diff/sqrt(N))

print(kable(table_1, caption="Table 1"))

model_prob = rep(prob_gam, length(N))
# bias = model_prob - freq_sin # FCY: not sure hot to compute the theoretical probability to remain as single in multiple choice situation
bias = prob_sin - freq_sin
table_2 = data.frame(n=N, model_single_prob = model_prob, mean_single_prob_calc_sim = prob_sin,  
                     mean_single_freq_sim = freq_sin, freq_bias = bias)

print(kable(table_2, caption="Table 2"))
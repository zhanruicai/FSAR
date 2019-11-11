
### get power-law network W
getPowerLawW<-function(N, alpha, normalize = T)           
{
  Nfollowers = rpldis(N, 1, alpha)                  ### generate N random numbers following power-law(1, alpha): k1-kN
  A = sapply(Nfollowers, function(n) {              ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  diag(A) = 0
  ind = which(rowSums(A)==0)                        ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1             ### for those node, randomly select 3 followees
  }
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  return(W)
}

### get stochastic block network W
getBlockW<-function(N, Nblock, normalize = T)                                                 
{
  if (N%%Nblock==0){                                                                 ### if N mod Nblock is integer
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)      ### obtain the diagnal block list
    ### generate following relations within the blocks
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),   
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    ### if N mod Nblock is not integer
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1) 
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    ### generate following relations within the blocks
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),     
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.9*N^{-1}),   
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)   ### combine the blocks in matrix
  ### to calculate the index of the off digonal indexes
  offDiag = which(isDiag == 0, arr.ind = T)                                         
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  ### people between blocks have 0.3 prob to follow
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                        
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  ################ transform bA to be a symmetric matrix ##############################################
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  
  ind = which(rowSums(bA)==0)                                 ### in case some row sums are zero
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                      ### for those node, randomly select 3 followees
  }
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                          ### row normalize bA
  #W = as(W, "dgCMatrix")
  return(W)
}

### Data Generating functions
get.func.network = function(n, min.points, max.points, net.type = "SBM", Y_out_point = NULL)
{
  ### Generate network structure, either stochastic block model or power law model
  if(net.type=="SBM") {W = getBlockW(n, 10)} else {W = getPowerLawW(n, alpha = 2.5, normalize = T)}
  
  #### Generate time points according to Yao JASA
  s = seq(from = 0, to = 10, length.out = max.points)
  
  #### Generate functional part
  ep = matrix(rnorm(n*max.points, sd = 0.5), n)
  b1 = matrix((s-5)^2/25, nrow = 1)
  b2 = matrix(s^0.5/4, nrow = 1)
  beta.seq = rbind(b1, b2)

  v1 = matrix(rnorm(n*2), n, 2)
  xi1 = matrix(rnorm(n, 0, 2), n, 1)
  xi2 = matrix(rnorm(n, 0, 1), n, 1)
  phi1.seq = matrix(+cos(2*s*pi/10) / sqrt(5), nrow = 1)
  phi2.seq = matrix(-sin(2*s*pi/10) / sqrt(5), nrow = 1)
  right.mat = v1%*%beta.seq+xi1%*%phi1.seq+xi2%*%phi2.seq+ep
  
  ### Generate rho matrix
  rho.mat = matrix(rep((15-s)/40, each = n), nrow = n, ncol = max.points)
  Y.all = right.mat+rho.mat*(W%*%right.mat)+(rho.mat^2)*(W%*%W%*%right.mat)+(rho.mat^3)*(W%*%W%*%W%*%right.mat)
  
  if (!is.null(Y_out_point))
  {
    Y_true = get.Ytrue(xi1, xi2, s = Y_out_point, W, v1)
  }else{
    Y_true = NULL
  }
  
  #### Sample and store the tragetories
  Y = matrix(NA, n, max.points)
  Tpoints = matrix(NA, n, max.points)
  
  for (i in 1:n) {
    total.points = sample(min.points:max.points, 1, prob = rep(1:(max.points-min.points+1)/(max.points-min.points+1)))
    ind = sort(sample(1:max.points, size = total.points, replace = F))
    Tpoints[i, 1:total.points] = s[ind]
    Y[i, 1:total.points] = Y.all[i,ind]
  }
  list(Y = Y, Tpoints = Tpoints, v1 = v1, W = W, Y.all = Y.all, Y_true = Y_true, time_seq = s)
}

obtain.kernel = function(Y, Tpoints, t)
{
  n = dim(Tpoints)[1]
  max.points = dim(Tpoints)[2]
  Ys = Tpoints[1,] 
  h = 1*density(Ys[Ys<999])$bw
  ### Store kernel
  K = matrix(999, n, max.points)
  for (i in 1:n) {
    ti = Tpoints[i,]
    total.points = sum(ti<999)
    K[i, 1:total.points] = gaussian.kernel.mat(t-ti[1:total.points], h)
  }
  list(K = K, h = h)
}

est_entire_timeline <- function(Y, X, Tpoints, W, num.points = 20)
{
  tlsub = seq(from = 1.5, to = 8.5, length.out = num.points)
  re <- matrix(0, nrow = num.points, ncol = 3)
  for (ii in 1:num.points) {
    #cat(ii, "\r")
    re[ii,] = est_rho_beta(tlsub[ii], Y, X, Tpoints, W)
  }
  list(est_coef = re, timeline = tlsub)
}

est_covariance <- function(est_coef, timeline, Y, X, Tpoints, W)
{
  N = dim(W)[1]
  rhott = est_coef[,1]
  betatt = est_coef[,2:dim(est_coef)[2]]
  num.points = length(timeline)
  Stt = lapply(rhott, function(rho){diag(N)-rho*W})
  CO = matrix(0, num.points, num.points)
  ns = rowSums(!is.na(Tpoints))
  
  for (s in 1:num.points) {
    for(t in 1:num.points){
      ## Calculate the numerator of Vst
      # kernel matrices for timepoint s
      T_diff = Tpoints - timeline[s]
      h = rule.thumb(T_diff[1,])/2.5
      Ker_mat_s = gaussian.kernel(T_diff, h)
      Y_s = rowSums(Y*Ker_mat_s, na.rm = T)
      # kernel matrices for timepoint t
      T_diff = Tpoints - timeline[t]
      h = rule.thumb(T_diff[1,])/2.5
      Ker_mat_t = gaussian.kernel(T_diff, h)
      Y_t = rowSums(Y*Ker_mat_t, na.rm = T)
      Y2_st = rowSums(Y^2*Ker_mat_s*Ker_mat_t, na.rm = T)
      Vst.n = Y_s%*%t(Y_t) - diag(Y2_st)
      ## Calculate the denominator of Vst
      K_s = rowSums(Ker_mat_s, na.rm = T)
      K_t = rowSums(Ker_mat_t, na.rm = T)
      K_st = rowSums(Ker_mat_s*Ker_mat_t, na.rm = T)
      Vst.d = K_s%*%t(K_t) - diag(K_st)
      ## Obtain Vst, which is the kernel approximation to Y(s)%*%t(Y(t))
      Vst = Vst.n/Vst.d
      suu=0
      for(i in 1:N)
      {
        suu = suu+t(Stt[[s]][i,])%*%Vst%*%(Stt[[t]][i,])
      }
      CO[s,t] = (suu - t(X%*%betatt[s,])%*%(X%*%betatt[t,]))/N
    }
  }
  list(cov = CO)
}

est_rho_beta <- function(t_est, Y, X, Tpoints, W, Ytrue = NULL, verbose = F)
{
  # dimensions
  N = nrow(Y)
  
  # kernel matrices
  T_diff = Tpoints - t_est
  h = rule.thumb(T_diff[1,])/3
  Ker_mat = gaussian.kernel(T_diff, h)
  ns = rowSums(!is.na(Tpoints))
  fs = rowSums(Ker_mat, na.rm = T)/ns # density estimation
  #fs = sum(Ker_mat, na.rm = T)/sum(ns)
  
  Y_Ker_mat = Y*Ker_mat/ns/fs
  Y2_Ker_mat = Y^2*Ker_mat/ns/fs
  Y_Ker_vec = rowSums(Y_Ker_mat, na.rm = T) # Y_Ker_vec is kernel approximation to Y(t)
  Y2_Ker_vec = rowSums(Y2_Ker_mat, na.rm = T) # Y2_Ker_vec is kernel approximation to Y^2
  
  zeta = apply(Y*Ker_mat, 1, function(s){return(mean(s, na.rm = T))})
  zetat = t((zeta/fs)%*%t(matrix(1, nrow = N, ncol = 1)))

  # nu_mat is kernel approximation to Y%*%t(Y)
  nu_mat = Y_Ker_vec%*%t(Y_Ker_vec)
  #diag(nu_mat) = Y2_Ker_vec
  
  # some basic matrices
  I_N = Diagonal(n = N, x = 1)
  tWW = t(W)%*%W
  d_tWW = diag(tWW)
  tW_W = t(W) + W
  
  rho = 0.5
  del = 1
  iter = 0
  while (max(abs(del))>10^{-4} &  iter < 20)
  {
    #### Update beta
    St = diag(N)-rho*W
    Mt = t(St)%*%St
    Dt = diag(diag(Mt)^{-1})
    Ut = Dt%*%Mt ### Ui is the i-th column of Ut
    A = Dt%*%t(St)%*%X
    beta = solve(t(A)%*%A)%*%(t(A)%*%rowSums(Ut*zetat))   
    
    if (verbose)
      cat(del, "\n")
    
    #### Update rho
    # S matrix and related derivatives
    S = I_N - rho*W
    tSS = t(S)%*%S
    tSS_rho = -tW_W + 2*rho*tWW # 1st order derivative of tSS
    tSS_rho2 = 2*tWW # 2nd order derivative of tSS
    
    # D matrix and related derivatives
    D = 1/diag(tSS)
    D1 = -2*rho*D^2*diag(tWW)
    D2 = -2*D^2*d_tWW + 8*rho^2*D^3*d_tWW^2
    
    # U_mat and u_mat u1_mat
    U_mat = D*tSS
    U1_mat = D1*tSS - D*tW_W + 2*rho*D*tWW
    u_mat = t(U1_mat)%*%U_mat
    U2_mat = D2*tSS + 2*D1*tSS_rho+D*tSS_rho2
    u1_mat = t(U2_mat)%*%U_mat + t(U1_mat)%*%U1_mat
    
    # xi vectors and v_vec
    Xbeta = X%*%beta
    xi = D*t(S)%*%Xbeta
    xi_rho = (D1*t(S) - D*t(W))%*%Xbeta
    xi_rho2 = (D2*t(S) - 2*D1*t(W))%*%Xbeta
    xi_beta = D*t(S)%*%X
    v_vec_rho = t(U1_mat)%*%xi + t(U_mat)%*%xi_rho
    v1_vec_rho = t(U2_mat)%*%xi + 2*t(U1_mat)%*%xi_rho + t(U_mat)%*%xi_rho2
    
    # first order derivatives
    # nu_mat = Ytrue%*%t(Ytrue)
    # Y_Ker_vec = Ytrue
    G1_rho = 2/N*sum(nu_mat*u_mat)
    G2_rho = -2/N*sum(Y_Ker_vec*v_vec_rho)
    G3_rho = 2/N*sum(xi_rho*xi)
    G_rho = G1_rho + G2_rho + G3_rho
    
    # # use F matrix
    # F_vec = D*t(S)%*%(S%*%Ytrue - Xbeta)
    # F_rho = (D1*t(S) - D*t(W))%*%(S%*%Ytrue - Xbeta)-D*t(S)%*%W%*%Ytrue
    # mean(F_vec*F_rho)*2
    
    # second order derivative
    
    H1_rho = 2/N*sum(nu_mat*u1_mat)
    H2_rho = -2/N*sum(Y_Ker_vec*v1_vec_rho)
    H3_rho = 2/N*sum(xi_rho*xi_rho + xi_rho2*xi)
    H_rho = H1_rho + H2_rho + H3_rho
    del = H_rho^{-1}*G_rho
    rho = rho - del
    iter = iter + 1
    if(rho<0 | rho>1) {rho = runif(1)}
    
  }
  
  return(c(rho, beta))
  
}

get.Ytrue <-function(xi1, xi2, s, W, v1)
{
  max.points = length(s)
  n = nrow(W)
  beta.seq = matrix(s/5, nrow = 1)
  phi1.seq = matrix(-cos(pi*s/10)/sqrt(5), nrow = 1)
  phi2.seq = matrix(sin(pi*s/10)/sqrt(5), nrow = 1)
  right.mat = v1%*%beta.seq+xi1%*%phi1.seq+xi2%*%phi2.seq
  ### Generate rho matrix
  rho.mat = matrix(rep((15-s)/40, each = n), nrow = n, ncol = max.points)
  Ytrue = right.mat+rho.mat*(W%*%right.mat)+(rho.mat^2)*(W%*%W%*%right.mat)+(rho.mat^3)*(W%*%W%*%W%*%right.mat)
  return(Ytrue)
}

### Get true covariance C(s,t)
getCovCST = function(s)
{
  phi1.seq = matrix(+cos(2*s*pi/10) / sqrt(5), ncol = 1)
  phi2.seq = matrix(-sin(2*s*pi/10) / sqrt(5), ncol = 1)
  TCST = 4*phi1.seq%*%t(phi1.seq)+phi2.seq%*%t(phi2.seq)
  TCST
}

### Define kernel functions
gaussian.kernel.mat = function(u, h)
{
  ### Input t is a vector
  u = u/h
  return(exp(-0.5*u^2)/sqrt(2*pi))
}


rule.thumb = function(x)
{
  1.06*sd(x, na.rm = T)*length(x)^{-1/5}
}

gaussian.kernel = function(u, h)
{
  ### Input t is a vector
  u = u/h
  return(exp(-0.5*u^2)/sqrt(2*pi))
}

### format function for keeping 2 digits after
specify_decimal = function(x, k) gsub('\\s+','',format(round(x, k), nsmall=k))    

rowSd.K<-function(x, K = 4){
  specify_decimal(apply(x, 1, sd), K)
}

rowMeans.K<-function(x, K = 4){
  specify_decimal(rowMeans(x, na.rm = T), K)
}
rowMedian.K<-function(x, K = 4){
  specify_decimal(apply(x, 1, median, na.rm = T), K)
}






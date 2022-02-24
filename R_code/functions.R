triu = function(A, diag = FALSE){
  n = nrow(A)
  Iu = matrix(1, nrow = n, ncol = n)
  if(diag){
    for(i in 2:n){
      Iu[i, 1:(i-1)] = 0
    }
  }else{
    for(i in 1:n){
      Iu[i, 1:i] = 0
    }
  }
  A = Iu*A
  return(A)
}

choosek = function(A, D, Kmax, d, At){
  # A: smoothed adjacency matrix from iteration of PisCES with dimension N*N
  # D: degree matrix with dimension N*N
  # Kmax: maximum number of modules
  # At: original adjacency matrix with dimension N*N
  
  # output: K: number of modules
  
  S = dim(A)
  n = S[1]
  if(identical(D, diag(n))){
    ad = diag(rowSums(A) + 1e-10)
    erw = eigen(D - (sqrt(solve(ad))%*%A%*%sqrt(solve(ad))), only.values = T)$values
  }else{
    erw = eigen(diag(n) - A, only.values = T)$values
  }
  erw = sort(erw)
  gaprw = erw[2:(Kmax+1)] - erw[1:Kmax]
  pin = sum( sqrt(D)%*%At%*%sqrt(D) )/choose(n, 2)/2
  
  threshold = 3.5/(pin^0.58)/(n^1.15)
  
  i = which(gaprw > threshold);
  
  K = max(i)
  if(length(i) == 0){
    K = 1
  }
  return(K)
}


PisCES = function(A, degree = TRUE, alpha = NULL, Kmax = NULL, T_n = 50){
  # A: adjacency matrices with dimention N*N*M, N is the number of nodes and 
  # M is the number of networks
  # degree: "TRUE" for degree correction
  # alpha: tuning parameter, default as 0.1
  # Kmax: maximum number of communities, default as N/10
  # T_n: number of iteration of pisces, default is 50
  
  S = dim(A)
  M = S[3]; N = S[2]
  if(is.null(alpha)){
    alpha = 0.1*matrix(data = 1, nrow = M, ncol = 2)
  }
  if(is.null(Kmax)){
    Kmax = floor(N/10)
  }
  
  Z = matrix(0, nrow = M, ncol = N)
  k = matrix(0, nrow = M, ncol = 1) + Kmax;
  obj = matrix(0, nrow = T_n, 1)
  V = array(0, c(N, Kmax, M));
  R = array(0, c(N, Kmax, M));
  D = A;
  Aori = A;
  
  # Degree correction
  if(degree){
    for(t in 1:M){
      a = A[, , t]
      d = diag(rowSums(abs(a)) + 1e-10)
      A[, , t] = sqrt(solve(d))%*%a%*%sqrt(solve(d))
      D[, , t] = d
    }
  }else{
    for(t in 1:M){
      D[, , t] = diag(N)
    }
  }
  
  # initial V and R

  for(t in 1:M){
    a = A[, , t]
    k[t] = choosek(A = a, D = D[, , t], Kmax = Kmax, d = 0.01, At = a)
    V[,1:k[t], t] = eigen(a)$vectors[, 1:k[t]]
  }
  for(trial in 1:T_n){
    R_old = R;
    V_old = V;
    for(t in 1:M){
      if(t == 1){
        X = V_old[, 1:k[t+1], t+1];
        At = A[, , t]
        S = A[, , t] + alpha[t, 2]*X%*%t(X)
        k[t] = choosek(S, D[, , t], Kmax, 0.01, At)
        V[ , 1:k[t], t] = eigen(At)$vectors[, 1:k[t]]
        obj[trial] = obj[trial] + 
          sum(abs(eigen(t(V[, 1:k[t], t]) %*% V_old[, 1:k[t], t])$values))
      }else{
        if(t == M){
          X = V_old[, 1:k[t-1], t-1]
          At = A[, , t]
          S = At + alpha[t, 1]*X%*%t(X)
          Ss = At%*%t(At) + 2*alpha[t, 1]*X%*%t(X)
          k[t] = choosek(S, D[, , t], Kmax, 0.01, At)
          eigen.temp = eigen(S);
          V[, 1:k[t], t] = eigen.temp$vectors[, 1:k[t]]
          e = eigen.temp$values[1:k[t]]
          obj[trial] = obj[trial] + 
            sum(abs(eigen(t(V[, 1:k[t], t]) %*% V_old[, 1:k[t], t])$values))
        }else{
          X1 = V_old[, 1:k[t-1], t-1]
          X2 = V_old[, 1:k[t+1], t+1]
          At = A[, , t]
          S = At + alpha[t, 1]*X1%*%t(X1) + alpha[t, 2]*X2%*%t(X2)
          Ss = At%*%t(At) + alpha[t, 1]*X1%*%t(X1) + alpha[t, 2]*X2%*%t(X2)
          k[t] = choosek(S, D[, , t], Kmax, 0.01, At)
          eigen.temp = eigen(S);
          V[, 1:k[t], t] = eigen.temp$vectors[, 1:k[t]]
          e = eigen.temp$values[1:k[t]]
          obj[trial] = obj[trial] + 
            sum(abs(eigen(t(V[, 1:k[t], t]) %*% V_old[, 1:k[t], t])$values))
        }
      }
    }
    obj_temp = obj[trial]
    cat('trial ', trial, ': have ', obj_temp, '\n')
    if(trial > 1){
      if(abs(obj[trial] - obj[trial - 1]) < 0.001){
        break;
      }
    }
  }
  if(obj[trial] - obj[trial-1] >= 0.001){
    stop('Method does not converage for alpha, please try smaller alpha')
  }
  Z[1, ] = kmeans(V[, 1:k[1], 1], k[1])$cluster
  for(i in 2:(M-1)){
    Z[i, ] = kmeans(V[, 1:k[i], i], k[i])$cluster
  }
  Z[M, ] = kmeans(V[, 1:k[M], M], k[M])$cluster
  return(Z)
}


eigen_complet = function(A, cv.idx, epsilon, k){
  M = matrix(0, nrow = nrow(A), ncol = ncol(A))
  while(TRUE){
    A2 = A*(1-cv.idx) + M*cv.idx
    temp = svd(A2)
    U = temp$u; V = temp$v; s = temp$d;
    s[(k+1):length(s)] = 0
    S = diag(s);
    M2 = U%*%S%*%t(V)
    M2[M2>1] = 1; 
    M2[M2<0] = 0;
    if(sqrt(sum((M-M2)^2)) < epsilon){
      break;
    }
    M = M2
  }
  return(A2)
}

wLoss = function(Adjtest, Adjtrain, zhat, type, cv.idx){
  # Adjtest: test matrix with dimension N*N; training edges with value 0;
  # Adjtrain: training matrix with dimension N*N; test edges with value 0;
  # zhat: estimated community assignment with dimension 1*N
  # type: modulatiry for 1 and log-likelihood for 2;
  # cv.idx: N*N matrix indicates the index of test edges: 1 for test and 0 for training
  
  # loss: modulatiry or log-likelihood on test data
  
  k = max(zhat)
  loss = 0
  M = mean(Adjtrain)
  minvalue = 0
  maxvalue = max(Adjtrain)
  
  n = length(zhat)
  W = sum(Adjtest)
  if(type){
    Abin = Adjtrain
    Kts = colSums(Abin);
    W = sum(Kts)
    
    row.col = which(cv.idx==0, arr.ind = TRUE)
    Num = hist(row.col[, 2], n, plot = FALSE)$counts
    NE = sum(Num)
    
    for(k1 in 1:n){
      for(k2 in 1:n){
        if((cv.idx[k1, k2] > 0) && (zhat[k1] == zhat[k2]) && (k1 != k2) ){
          loss = loss + (Adjtest[k1, k2] - Kts[k1]/Num[k1] * Kts[k2]/Num[k2])/W * NE
        }
      }
    }
  }else{
    Abin = Adjtrain
    Kts = colSums(Abin)
    W = sum(Kts)
    hat0 = matrix(0, nrow = k, ncol = k)
    theta = matrix(0, nrow = n, ncol = 1)
    for(i in 1:k){
      for(j in 1:k){
        Ak = Abin[zhat == i, zhat == j]
        Aj = Adjtrain[zhat == i, zhat == j]
        ajx = cv.idx[zhat == i, zhat == j]
        hat0[i, j] = sum(Ak)
      }
    }
    for(k1 in 1:n){
      kk = zhat[k1]
      theta[k1] = Kts[k1]/sum(hat0[kk, ])
    }
    for(k1 in 1:n){
      for(k2 in 1:n){
        if(cv.idx[k1, k2] > 0){
          prob = theta[k1]*theta[k2]*hat0[zhat[k1], zhat[k2]]/4*5
          
          if(is.na(prob)|prob == 0){
            prob = 1e-5
          }
          if(prob >= 1){
            prob = 1 - 1e-5
          }
          loss = loss - log(prob)*(Adjtest[k1, k2]>0.7^6) - 
            log(1-prob)*(Adjtest[k1, k2] <= 0.7^6)
        }
      }
    }
    loss = loss/W
  }
  return(loss)
}


PisCES.cv = function(A, alphalist = NULL, degree = TRUE, Kmax = NULL, T_n = 50, Kfold = 5){
  S = dim(A)
  M = S[3]; N = S[2];
  
  if(is.null(alphalist)){
    alphalist = c(0.05, 1)
  }
  
  if(is.null(Kmax)){
    Kmax = floor(N/10)
  }
  
  n.idx = sample(N^2)
  Atrain = array(0, c(N, N, M, Kfold))
  Atrain2 = Atrain
  
  for(t in 1:M){
    cv.idx = array(0, c(N, N, Kfold))
    for(k in 1:Kfold){
      at = A[, , t]
      cv.idx.temp = matrix(0, ncol = N, nrow = N)
      test = n.idx[((k-1)*floor(N^2/Kfold) + 1):floor(N^2/Kfold)*k]
      at[test] = 0;
      cv.idx.temp = triu(cv.idx.temp) + t(triu(cv.idx.temp))
      cv.idx[, , k] = cv.idx.temp;
      at = triu(at) + t(triu(at))
      Atrain[, , t, k] = at;
      Atrain2[, , t, k] = eigen_complet(A = at, cv.idx = cv.idx.temp, epsilon = 10, k = 10)
    }
  }
  l = length(alphalist)
  modu = matrix(0, nrow = 1, ncol = l)
  like = matrix(0, nrow = 1, ncol = l)
  
  for(a in 1:l){
    alpha = matrix(1, nrow = M, ncol = 2) * alphalist[a]
    Z = array(0, c(M, N, Kfold))
    for(k in 1:Kfold){
      PisCES(Atrain2[, , , k], degree, alpha, Kmax, T_n)
      for(t in 1:M){
        modu[a] = modu[a] + wLoss(A[, , t], Atrain[, , t, k], Z[t, , k], 1, cv.idx[, , k])
        like[a] = like[a] + wLoss(A[, , t], Atrain[, , t, k], Z[t, , k], 2, cv.idx[, , k])
      }
    }
    cat('modularity for alpha = ', alphalist[a], ' is ', modu[a], '\n')
    cat('loglikelihood for alpha = ', alphalist[a], ' is ', like[a], '\n')
  }
  return(list(modularity = modu, log.likelihood = like))
}


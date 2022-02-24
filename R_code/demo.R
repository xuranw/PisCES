N <- 200; M <- 6; K <- 4;
P <- diag(K)*0.3 + matrix(1, K, K)*0.1
A <- array(0, c(N, N, M));  
z = rep(1:K, floor(N/K))
for(t in 1:M){
  a = matrix(runif(N^2), N, N) < P[z, z] + 0
  A[, , t] = triu(a) + t(triu(a))
}

Z = PisCES(A, alpha = matrix(0.1, nrow = M, ncol = 2)) 

alphalist = c(0.05, 0.1, 0.15)
cv.temp = PisCES.cv(A = A, alphalist = alphalist)

library(earth) 
num_split <- 50
n <-400
p <- 500
p0 <- 10#25
q <- 0.25
delta <- 10
amountTest=0.5
amountTrain=0.5
signal_index <- sample(c(1:p), size = p0, replace = F)
# simulate data
n1 <- floor(n/2); n2 <- n - n1
X1 <- matrix(rnorm(n1*p, mean= 1), n1, p)
X2 <- matrix(rnorm(n2*p, mean=-1), n2, p)
X  <- rbind(X1, X2)
beta_star <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*1000
y <- (X^2 %*% beta_star+ rnorm(n))
n <- nrow(X); p <- ncol(X)
data <- data.frame(cbind(y, X))
colnames(data) <- c("y", paste0("X", 1:p))
size_half <- floor((amountTest/2) * n)
train_index <- sample.int(n, size = floor(amountTrain * n), replace = FALSE)
remaining_index <- setdiff(seq_len(n), train_index)
sample_index1 <- sample(remaining_index, size = size_half, replace = FALSE)
sample_index2 <- setdiff(remaining_index, sample_index1)
dataTrain <- data[train_index, , drop = FALSE]
#mars_poly= earth(
#  y ~ .,
#  data    = dataTrain,
#  nk      = 200,      # MUCH larger cap on #basis functions
#  fast.k  = 0,        # turn off fast MARS (more exhaustive forward step)
# pmethod = "seqrep", 
#  nprune  = 150,       # target a richer final model
#  penalty = -1,       # negative/low penalty -> prefers *more* terms
# thresh  = 0.00001,     # small threshold -> easier to add new terms
#  trace   = 0,
#  minspan=-1
#)
mars_poly= earth(
  y ~ .,
  data    = dataTrain
)
lm <- mars_poly
pred1 <- predict(lm, newdata = as.data.frame(X[sample_index1, , drop = FALSE]))
pred2 <- predict(lm, newdata = as.data.frame(X[sample_index2, , drop = FALSE]))
y1 <- y[sample_index1]; y2 <- y[sample_index2]
R2orig1 <- 1 - sum((y1 - pred1)^2) / sum((y1 - mean(y1))^2)
R2orig2 <- 1 - sum((y2 - pred2)^2) / sum((y2 - mean(y2))^2)
R2orig1
R2orig2



Rnew1 <- sapply(seq_len(p), function(j)
  permR2Mars(as.data.frame(X[sample_index1, , drop = FALSE]), Y = y1, j = j, model = lm))
Rnew2 <- sapply(seq_len(p), function(j)
  permR2Mars(as.data.frame(X[sample_index2, , drop = FALSE]), Y = y2, j = j, model = lm))
beta1 <- R2orig1 - Rnew1
beta2 <- R2orig2 - Rnew2
mirror <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
hist(mirror)
selected_index <- SelectFeatures(mirror, abs(mirror), q)
res <- CalculateFDP_Power(selected_index, signal_index)
fdp_val <- res$fdp
pow_val <- res$power

fdp_val
pow_val
sum(signal_index%in%evimp(mars_poly, trim=T)[,1])/p0

#summary(mars_poly, digits = 2, style = "pmax")




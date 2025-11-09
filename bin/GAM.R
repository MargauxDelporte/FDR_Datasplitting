library(mgcv)
### algorithmic settings
num_split <-5# 50
n <-400
p <- 500
p0 <- 25
q <- 0.1
delta=10
#set.seed(124)(123) i=5
set.seed(456)
signal_index <- sample(c(1:p), size = p0, replace = F)
n1 <- floor(n/2); n2 <- n - n1
X1 <- matrix(rnorm(n1*p, mean= 1), n1, p)
X2 <- matrix(rnorm(n2*p, mean=-1), n2, p)
X  <- rbind(X1, X2)
beta_star <- numeric(p)
beta_star[signal_index] <- rnorm(p0, 0, delta*sqrt(log(p)/n))*1000
y <- (X^2 %*% beta_star + rnorm(n))
# Optional: nonlinear screening to keep model size sane with n=30
# (does not change data; just selects variables)
library(energy)
dc <- sapply(as.data.frame(X), function(z) dcor(z, y))   # distance correlation
keep <- order(dc, decreasing = TRUE)[1:min(30, ncol(X))] # keep top ~20–40
Xd <- as.data.frame(X[, keep, drop=FALSE])
dat <- data.frame(y=y, Xd)

# Build additive model; 'bs="ts"' gives shrinkage so terms can be zeroed out
form <- as.formula(
  paste("y ~", paste(sprintf("s(%s, k=3, bs='ts')", names(Xd)), collapse=" + "))
)
fit <- gam(form, data=dat, method="REML", select=TRUE)  # select=TRUE adds extra sparsity
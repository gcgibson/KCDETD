
#a
data <- read.table("data.dat",header = T)

X <- cbind(data$x1,data$x2,data$x3,data$x4)
y <- data$y


beta_hat <- solve(t(X)%*%X)%*%t(X)%*%y
y_bar <- mean(y)
n <- nrow(X)
k <- 4

SSR <- t(beta_hat)%*%t(X)%*%y - n*y_bar^2
SSE <- t(y)%*%y- t(beta_hat)%*%t(X)%*%y


F_stat <- (SSR/k)/(SSE/(n-k-1))

pf(F_stat,df1=k,df2=n-k-1,lower.tail = F)


#b


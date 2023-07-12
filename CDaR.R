rm(list=ls())

library(parma)
library(xts)
library(CVXR)

data("etfdata")
Data <- etfdata/lag(etfdata) - 1
Data <- na.omit(Data)

#parma code

spec <- parmaspec(scenario = Data, forecast = colMeans(Data), risk = "CDaR", target = mean(colMeans(Data)), 
                  targetType = "equality", riskType =  "minrisk", LB = rep(0, 15), UB = rep(1, 15), budget = 1)
solution <- parmasolve(spec, type = "LP")
weight_parma <-  solution@solution$weights

#CVXR code

a <- 0.05
A <- diag(2271)
for (i in 1:2270){
  A[i+1, i] = -1
}

weight <- Variable(15)
z <- Variable(2271)
v <- Variable(1)
u <- Variable(2271)

objective <- Minimize(v+1/(2271*a)*sum(z))
constraints <- list(z-u+v>=0, as.matrix(Data)%*%weight+A%*%u>=0, z>=0, u>=0, weight >= 0, sum(weight) == 1, t(colMeans(Data))%*%weight == mean(colMeans(Data)))
problem <- Problem(objective, constraints = constraints)

weight_cvxr <- solve(problem, solver = "GLPK_MI")
weight_cvxr <- weight_cvxr$getValue(weight)

#comparison

risk_parma <-  riskfun(weights = weight_parma, Data = Data, risk = c("CDaR"))
risk_cvxr <-  riskfun(weights = weight_cvxr, Data = Data, risk = c("CDaR"))

target_parma <-  t(colMeans(Data))%*%weight_parma
target_cvxr <-  t(colMeans(Data))%*%weight_cvxr

comparison <-  data.frame(matrix(c(weight_parma, risk_parma, target_parma, weight_cvxr, risk_cvxr, target_cvxr), nrow = 2, byrow = TRUE))

colnames(comparison) <-  c(colnames(Data), c("Risk", "Target Value"))
rownames(comparison) <-  c("parma","CVXR")

print(comparison)

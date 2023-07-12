rm(list=ls())

library(parma)
library(xts)
library(CVXR)

data("etfdata")
Data <- etfdata/lag(etfdata) - 1
Data <- na.omit(Data)

#parma code

spec <- parmaspec(scenario = Data, forecast = colMeans(Data), risk = "LPM", target = mean(colMeans(Data)), 
                  targetType = "equality", riskType =  "minrisk", options = list(moment = 1, threshold = 999), LB = rep(0, 15), UB = rep(1, 15), budget = 1)
solution <- parmasolve(spec)
weight_parma <-  solution@solution$weights

#CVXR code

ones <- matrix(1, 2271, 1)

a <- 1
tau <- mean(colMeans(Data))

v <- Variable(2271)
weight <- Variable(15)

objective <- Minimize(sum(v^a))
constraints <- list(weight >= 0, sum(weight) == 1, t(colMeans(Data))%*%weight == mean(colMeans(Data)), v >= 0, v >= tau*ones - as.matrix(Data)%*%weight)
problem <- Problem(objective, constraints = constraints)

weight_cvxr <- solve(problem)
weight_cvxr <- weight_cvxr$getValue(weight)
print(weight_cvxr)
#comparison

risk_parma <-  riskfun(weights = weight_parma, Data = Data, risk = c("LPM"))
risk_cvxr <-  riskfun(weights = weight_cvxr, Data = Data, risk = c("LPM"))

target_parma <-  t(colMeans(Data))%*%weight_parma
target_cvxr <-  t(colMeans(Data))%*%weight_cvxr

comparison <-  data.frame(matrix(c(weight_parma, risk_parma, target_parma, weight_cvxr, risk_cvxr, target_cvxr), nrow = 2, byrow = TRUE))

colnames(comparison) <-  c(colnames(Data), c("Risk", "Target Value"))
rownames(comparison) <-  c("parma","CVXR")

print(comparison)

rm(list=ls())

library(parma)
library(xts)
library(CVXR)

data("etfdata")
Data <- etfdata/lag(etfdata) - 1
Data <- na.omit(Data)

#parma code

spec <- parmaspec(scenario = Data, forecast = colMeans(Data), risk = "MAD", target = mean(colMeans(Data)), targetType = "equality", riskType = "minrisk", LB = rep(0, 15), UB = rep(1, 15), budget = 1)
solution <- parmasolve(spec, type = "LP")
weight_parma <-  solution@solution$weights

#CVXR code

weight <- Variable(15)
y <- Variable(2271)

Data_normal <-  as.matrix(sweep(Data, 2, t(colMeans(Data))))

objective <- Minimize(CVXR::norm1(Data_normal %*% weight))
constraints <- list(weight >= 0, sum(weight) == 1, t(colMeans(Data)) %*% weight == mean(colMeans(Data)))
problem <- Problem(objective, constraints = constraints)
result <- psolve(problem, solver = "GLPK")

weight_cvxr <-  result$getValue(weight)

#Comparison

risk_parma <-  riskfun(weights = weight_parma, Data = Data, risk = c("mad"))
risk_cvxr <-  riskfun(weights = weight_cvxr, Data = Data, risk = c("mad"))

target_parma <-  t(colMeans(Data))%*%weight_parma
target_cvxr <-  t(colMeans(Data))%*%weight_cvxr

comparison <-  data.frame(matrix(c(weight_parma, risk_parma, target_parma, weight_cvxr, risk_cvxr, target_cvxr), nrow = 2, byrow = TRUE))

colnames(comparison) <-  c(colnames(Data), c("Risk", "Target Value"))
rownames(comparison) <-  c("parma","CVXR")

print(comparison)

### source functions
source("twoStage.R") ## data generation used in Section 4.3
source("twoStage_nontree.R") ## data generation used in Section 4.3
source("oneStage.R") ## data generation used in Section 4.2
source("oneStage_nontree.R") ## data generation used in Section 4.2
source("twoStage2.R") ## data generation used in Section 4.1
source("twoStage2_nontree.R") ## data generation used in Section 4.1
source("general_functions.R") ## general functions used in simulations
source("restl_functions.R") ## functions used to implement ReST-L
source("trl_functions.R") ## functions used to implement T-RL

## This is an example of function calls for the simulations with the following specifications:
## Sample size, N = 350
## Test set size, N2 = 1000
## Monte carlo iterations, B = 500
## Correlation coefficient, rho = 0.2, reflecting how the H covariates were generated
## Number of variables in H, No.Var.H = 20. Options include (No.Var.H = 20, 50, 100)
## Number of variables in Hsub, No.Var.Hsub = 7. Options include (No.Var.Hsub = 7, 10, 20)
## Method used, Method = "REST-L". Options include ("REST-L", "T-RL", "Q-Linear", "Q-Linear-REST", "Q-NP", and "Q-NP-REST")
## Propensity model specification, prop.model = "correct". Options include ("correct", "incorrect")

set.seed(110)
results1 <- twoStage(N = 350, N2 = 1000, iter = 500, rho = 0.2, No.Var.H=20, No.Var.Hsub=7, Method="REST-L",prop.model = "correct")

set.seed(110)
results1b <- twoStage(N = 350, N2 = 1000, iter = 500, rho = 0.2, No.Var.H=20, No.Var.Hsub=7, Method="T-RL",prop.model = "correct")

set.seed(110)
results1c <- twoStage(N = 350, N2 = 1000, iter = 500, rho = 0.2, No.Var.H=20, No.Var.Hsub=7, Method="Q-Linear",prop.model = "correct")

set.seed(110)
results1d <- twoStage(N = 350, N2 = 1000, iter = 500, rho = 0.2, No.Var.H=20, No.Var.Hsub=7, Method="Q-Linear-REST",prop.model = "correct")

set.seed(110)
results1e <- twoStage(N = 350, N2 = 1000, iter = 500, rho = 0.2, No.Var.H=20, No.Var.Hsub=7, Method="Q-NP",prop.model = "correct")

set.seed(110)
results1f <- twoStage(N = 350, N2 = 1000, iter = 500, rho = 0.2, No.Var.H=20, No.Var.Hsub=7, Method="Q-NP-REST",prop.model = "correct")

results <- rbind(results1, results1b, results1c, results1d, results1e, results1f)
write.csv(results, file = "results_File101.csv")

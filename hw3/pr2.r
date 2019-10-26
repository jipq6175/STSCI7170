
library('lme4');

n = 5;
d = 2;
t = 4;

Y = c(1.08, 1.99, 1.46, 1.21, 1.48, 2.50, 2.62, 1.95,
      1.19, 2.10, 1.21, 0.96, 0.62, 0.88, 0.68, 0.48,
      1.22, 1.91, 1.36, 0.90, 0.65, 1.52, 1.32, 0.95,
      0.60, 1.10, 1.03, 0.61, 0.32, 2.12, 1.48, 1.09,
      0.55, 1.00, 0.82, 0.52, 1.48, 0.90, 0.75, 0.44)



Jn = rep(1, n) %*% t(rep(1, n)) / n
Cn = diag(n) - Jn
Jd = rep(1, d) %*% t(rep(1, d)) / d
Cd = diag(d) - Jd
Jt = rep(1, t) %*% t(rep(1, t)) / t
Ct = diag(t) - Jt

A1 = Jn %x% Jd %x% Jt
A2 = Cn %x% Jd %x% Jt
A3 = Jn %x% Cd %x% Jt
A4 = Cn %x% Cd %x% Jt
A5 = Jn %x% Jd %x% Ct
A6 = Cn %x% Jd %x% Ct
A7 = Jn %x% Cd %x% Ct
A8 = Cn %x% Cd %x% Ct



serum =  data.frame(matrix(ncol = 4, nrow = n*d*t))
colnames(serum) = c("response", "subject", "drug", "time")
serum$response = Y
serum$subject = 1:n %x% rep(1, d) %x% rep(1,t)
serum$drug = rep(1, n) %x% 1:d %x% rep(1,t)
serum$time = rep(1, n) %x% rep(1, d) %x% 1:t

ml = lmer(response ~ 1 + (1|subject) + (drug) + (time)
                       + (1|subject:drug) 
                       + (1|subject:time)
                       + (drug:time)
                       , data=serum, REML=FALSE)


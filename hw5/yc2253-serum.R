source('yc2253-em.r');

# For blood serem data

n = 5; # subject
d = 2; # treatment
t = 4; # time

Y = rep(0.0, n*d*t);

# load the serum.txt and organize it into Y vector
lines = readLines('serum.txt');
for(i in 1:n){
  s = as.numeric(unlist(strsplit(lines[i+4], ' ')));
  Y[(8*(i-1)+1):(8*i)] = s[2:length(s)];
}


# Construct the "vectors" for the regression matrix X = [x1 x2 x3 x4]
# And 
Jn = rep(1, n) %*% t(rep(1, n)) / n; 
Cn = diag(n) - Jn;
Jd = rep(1, d) %*% t(rep(1, d)) / d;
Cd = diag(d) - Jd;
Jt = rep(1, t) %*% t(rep(1, t)) / t;
Ct = diag(t) - Jt;

A1 = Jn %x% Jd %x% Jt;
A2 = Cn %x% Jd %x% Jt;
A3 = Jn %x% Cd %x% Jt;
A4 = Cn %x% Cd %x% Jt;
A5 = Jn %x% Jd %x% Ct;
A6 = Cn %x% Jd %x% Ct;
A7 = Jn %x% Cd %x% Ct;
A8 = Cn %x% Cd %x% Ct;


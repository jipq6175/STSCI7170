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


# Construct the "vectors" for the regression matrix X = [x1 x2 x3 x4] of fixed effects
# And Z = [z1 z2 z3 z4] of the random effects
x1 = rep(1, n) %x% rep(1, d) %x% rep(1, t); 
x2 = rep(1, n) %x% diag(d) %x% rep(1, t);
x3 = rep(1, n) %x% rep(1, d) %x% diag(t);
x4 = rep(1, n) %x% diag(d) %x% diag(t);

X = cbind(x1, x2, x3, x4);

z1 = diag(n) %x% rep(1, d) %x% rep(1, t);
z2 = diag(n) %x% diag(d) %x% rep(1, t);
z3 = diag(n) %x% rep(1, d) %x% diag(t);
z4 = diag(n) %x% diag(d) %x% diag(t); 

zlist = list(z1, z2, z3, z4);


result = vcem(x4, Y, zlist, runif(8, 1, 5), runif(4))
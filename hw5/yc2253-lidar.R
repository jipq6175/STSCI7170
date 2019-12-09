
source('yc2253-em.r');
# For the lidar data 
x = rep(0, 221); 
y = rep(0, 221); 
lines = readLines('lidar-1.txt'); 
for (i in 1:221){
  s = as.numeric(unlist(strsplit(lines[i+3], ' ')));
  x[i] = s[3];
  y[i] = s[4]; 
}


knots = sort(sample(x, 9)); 


# construct X
X = cbind(rep(1, length(x)), x, x^2); 

# construct z

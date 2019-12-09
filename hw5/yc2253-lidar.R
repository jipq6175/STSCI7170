
source('yc2253-em.r');
# For the lidar data 
x = rep(0, 221); 
y = rep(0, 221); 
lines = readLines('lidar-1.txt'); 
for (i in 1:221){
  s = as.numeric(unlist(strsplit(lines[i+3], ' ')));
  x[i] = s[3];
  if (is.na(s[4])){
    y[i] = s[5];
  }else{
    y[i] = s[4];
  }
}


knots = sort(sample(x, 9)); 


# construct X
X = cbind(rep(1, length(x)), x, x^2); 
XG = X; 
# construct z
zlist = list();
for (i in 1:9){
  tmp = x - knots[i]; 
  tmp[which(tmp <= 0)] = 0.0; 
  zlist[i] = list(tmp^2); 
  XG = cbind(XG, tmp^2); 
}
zlist[10] = list(diag(length(y))); 
yfit = XG %*% ginv(t(XG) %*% XG) %*% t(XG) %*% y; 
result = vcem(X, y, zlist, c(1, 1e-5, 1e-10), rep(1, 10), maxiter=1000);



plot(x, y);
plot(x, yfit, col='blue')
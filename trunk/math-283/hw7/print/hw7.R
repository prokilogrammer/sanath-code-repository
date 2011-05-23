# Read data from file
data = read.delim('hw7data.dat.txt', header=FALSE)
nr = nrow(data)
nc = ncol(data)
n_sample = 9
before = data[, 2:10];
after = data[, 11:19];

# Problem-a
### 1.a) Estimating value of SIGMA for Two Sample Z-test
###       - Find sample variance for each spot and average it to get SIGMA 
### 
z_xbar=1
z_ybar=1
for( r in c(1:nr))
{
  z_xbar[r] = mean( t(before[r,]) );
  z_ybar[r] = mean( t(after[r,]) ); 

}
sample_var_x = sum(  (before - z_xbar)^2 ) / (n_sample-1);
sample_var_y = sum(  (after - z_ybar)^2 ) / (n_sample-1);
z_var_x = sum(sample_var_x) / nr;
z_var_y = sum(sample_var_y) / nr;

print("Values of Standard deviation for patients BEFORE treatment is:")
print(sqrt(z_var_x));
print("Values of Standard deviation for patients AFTER treatment is:")
print(sqrt(z_var_y));

### 4.a) Estimating value of SIGMA for Paired Sample Z-test
###       - D = X - Y. Run test on D
###       - Find sample variance for each spot and average it to get SIGMA 
### 
D = before - after;
pz_xbar = c(1,nr);
for( r in c(1:nr))
{
  pz_xbar[r] = mean( t(D[r, 1:9]) );
}
sample_var_x = sum(  ( D[,1:9] - pz_xbar)^2 ) / (9-1);
pz_var_x = sum(sample_var_x) / nr;

print("Values of Standard deviation for Paired Data is:")
print(sqrt(pz_var_x));

# Problem-b
### 1.b) Two-sample Z-test
temp_denom =  sqrt( (z_var_x + z_var_y) / 9)
zscore = (z_xbar - z_ybar) / temp_denom;
z_pvalue = 2* pmin( pnorm(zscore), 1 - pnorm(zscore) );

### 2.b) Two sample t-test
t_pvalue=1;
t_stat = 1;
for( r in c(1:nr) )
{
  t_test = t.test( as.double(before[r,]), as.double(after[r,]), var.equal = TRUE);
  t_pvalue[r] = t_test$p.value;
  t_stat[r] = t_test$statistic;
}
### 3.b) Mann-Whitney Test
mw_pvalue=1;
mw_stat = 1;
for( r in c(1:nr) )
{
  mw_test = wilcox.test( as.double(before[r,]), as.double(after[r,]));
  mw_pvalue[r] = mw_test$p.value;
  mw_stat[r] = mw_test$statistic;
}

### 4.b) Paired Z Test
pz_score = pz_xbar / sqrt( pz_var_x / 9 )
pz_pvalue = 2* pmin( pnorm(pz_score), 1 - pnorm(pz_score) );

### 5.b) Paired t-test
pt_pvalue=1;
pt_stat = 1;
for( r in c(1:nr) )
{
  pt_test = t.test( as.double(before[r,]), as.double(after[r,]), paired=TRUE);
  pt_pvalue[r] = pt_test$p.value;
  pt_stat[r] = pt_test$statistic;
}

### 6.b) Wilcoxon Test
w_pvalue=1;
w_stat = 1;
for( r in c(1:nr) )
{
  #w_test = wilcox.test( as.double(after[r,]) - as.double(before[r,]), 0);
  w_test = wilcox.test(as.double(after[r,]), as.double(before[r,]), paired=TRUE);
  w_pvalue[r] = w_test$p.value;
  w_stat[r] = w_test$statistic;
}

## Combining pvalue into single arrays for convenience
z = cbind( id=data[,1], p=z_pvalue, stat=zscore);
t = cbind( id=data[,1], p=t_pvalue, stat=t_stat);
mw = cbind( id=data[,1], p=mw_pvalue, stat=mw_stat);
pz = cbind( id=data[,1], p=pz_pvalue, stat=pz_score);
pt = cbind( id=data[,1], p=pt_pvalue, stat=pt_stat);
w = cbind( id=data[,1], p=w_pvalue, stat=w_stat);


######### Problem-b : OUTPUT #######
printOutput <- function(mat)
{
  pos1 = which( data[,1]==1020, arr.ind=TRUE)
  pos2 = which( data[,1]==1021, arr.ind=TRUE)
  pos3 = which( data[,1]==1022, arr.ind=TRUE)
  pos4 = which( data[,1]==1023, arr.ind=TRUE)
  pos5 = which( data[,1]==1024, arr.ind=TRUE)
  print( c(1020, 1021, 1022, 1023, 1024) )
  print( c( mat[pos1,2], mat[pos2,2], mat[pos3,2], mat[pos4,2], mat[pos5,2]) );
  print( c( mat[pos1,3], mat[pos2,3], mat[pos3,3], mat[pos4,3], mat[pos5,3]) );
}

print( "-- Two Sample Z Test - P Value & Statistic--")
printOutput(z)

print( "-- Two Sample T Test - P Value & Statistic--")
printOutput(t)

print( "-- Mann-Whitney Test - P Value & Statistic--")
printOutput(mw)

print( "-- Paired Z Test - P Value & Statistic --")
printOutput(pz)

print( "-- Paired T Test - P Value & Statistic --")
printOutput(pt)

print( "-- Wilcoxon Signed Test - P Value & Statistic --")
printOutput(w)


#### Problem-c
print(" Top 12 significant P-values ")
print_top_dozen <- function(mat)
{
 perm = order( mat[,2] )
 print( (mat[perm, ])[1:12,1:2] )
}

print( "-- Two Sample Z Test --")
print_top_dozen(z)

print( "-- Two Sample T Test --")
print_top_dozen(t)

print( "-- Mann-Whitney Test --")
print_top_dozen(mw)

print( "-- Paired Z Test --")
print_top_dozen(pz)

print( "-- Paired T Test --")
print_top_dozen(pt)

print( "-- Wilcoxon Signed Test --")
print_top_dozen(w)


#### Problem-d
##### 1.d.i 

##### 1.d.ii
count_significance <- function(mat, alpha)
{
 print( alpha );
 print( length( which( mat[,2] <= alpha ) ) ); 
}
print(" # of genes at different significance levels");
print( "-- Two Sample Z Test --")
count_significance(z,0.01)
count_significance(z,0.05)
count_significance(z,0.10)

print( "-- Two Sample T Test --")
count_significance(t,0.01)
count_significance(t,0.05)
count_significance(t,0.10)


print( "-- Mann-Whitney Test --")
count_significance(mw,0.01)
count_significance(mw,0.05)
count_significance(mw,0.10)


print( "-- Paired Z Test --")
count_significance(pz,0.01)
count_significance(pz,0.05)
count_significance(pz,0.10)


print( "-- Paired T Test --")
count_significance(pt,0.01)
count_significance(pt,0.05)
count_significance(pt,0.10)


print( "-- Wilcoxon Signed Test --")
count_significance(w,0.01)
count_significance(w,0.05)
count_significance(w,0.10)

##### 1.d.iii
print(" # of genes at different significance levels with Bonferroni Correction");
print( "-- Two Sample Z Test --")
count_significance(z,0.01/nr)
count_significance(z,0.05/nr)
count_significance(z,0.10/nr)

print( "-- Two Sample T Test --")
count_significance(t,0.01/nr)
count_significance(t,0.05/nr)
count_significance(t,0.10/nr)


print( "-- Mann-Whitney Test --")
count_significance(mw,0.01/nr)
count_significance(mw,0.05/nr)
count_significance(mw,0.10/nr)


print( "-- Paired Z Test --")
count_significance(pz,0.01/nr)
count_significance(pz,0.05/nr)
count_significance(pz,0.10/nr)


print( "-- Paired T Test --")
count_significance(pt,0.01/nr)
count_significance(pt,0.05/nr)
count_significance(pt,0.10/nr)


print( "-- Wilcoxon Signed Test --")
count_significance(w,0.01/nr)
count_significance(w,0.05/nr)
count_significance(w,0.10/nr)


##### 1.d.iv
print(" # of genes at different significance levels with Sidak Correction");
print( "-- Two Sample Z Test --")
count_significance(z,1- (1-0.01)^(1/nr) )
count_significance(z,1- (1-0.05)^(1/nr) )
count_significance(z,1- (1-0.10)^(1/nr) )

print( "-- Two Sample T Test --")
count_significance(t,1- (1-0.01)^(1/nr) )
count_significance(t,1- (1-0.05)^(1/nr) )
count_significance(t,1- (1-0.10)^(1/nr) )


print( "-- Mann-Whitney Test --")
count_significance(mw,1- (1-0.01)^(1/nr) )
count_significance(mw,1- (1-0.05)^(1/nr) )
count_significance(mw,1- (1-0.10)^(1/nr) )


print( "-- Paired Z Test --")
count_significance(pz,1- (1-0.01)^(1/nr) )
count_significance(pz,1- (1-0.05)^(1/nr) )
count_significance(pz,1- (1-0.10)^(1/nr) )


print( "-- Paired T Test --")
count_significance(pt,1- (1-0.01)^(1/nr))
count_significance(pt,1- (1-0.05)^(1/nr))
count_significance(pt,1- (1-0.10)^(1/nr))


print( "-- Wilcoxon Signed Test --")
count_significance(w,1- (1-0.01)^(1/nr) )
count_significance(w,1- (1-0.05)^(1/nr) )
count_significance(w,1- (1-0.10)^(1/nr) )


#### Problem - e
range = which( pt[,2] <= 0.1 )
true_positives = pt[range, ]
fdr <- function(mat)
{
  rang = which( mat[,2] <= 0.1 ) 
  mismatch = setdiff( mat[rang,1],true_positives[,1] )
  print( length(mismatch)/ length(true_positives[,1]) )
}
print( "-- Two Sample Z Test: FDR--")
fdr(z)
print( "-- Two Sample T Test: FDR --")
fdr(t)
print( "-- Mann-Whitney Test: FDR --")
fdr(mw)
print( "-- Paired Z Test: FDR --")
fdr(pz)
print( "-- Paired T Test: FDR --")
fdr(pt)
print( "-- Wilcoxon Signed Test: FDR --")
fdr(w)


#### Problem-f
pdf("output-f.pdf")
hist(z[,3] ,col="blue1", main="Histogram of test statistic for two sample z-test", xlab='Test Statistic', ylab='number of occurrence', labels=TRUE)

hist(t[,3] ,col="blue1", main="Histogram of test statistic for two sample t-test", xlab='Test Statistic', ylab='number of occurrence', labels=TRUE)

hist(mw[,3] ,col="blue1", main="Histogram of test statistic for Mann Whitney U test", xlab='Test Statistic', ylab='number of occurrence', labels=TRUE)

hist(pz[,3] ,col="blue1", main="Histogram of test statistic for paired z-test", xlab='Test Statistic', ylab='number of occurrence', labels=TRUE)

hist(t[,3] ,col="blue1", main="Histogram of test statistic for paired t-test", xlab='Test Statistic', ylab='number of occurrence', labels=TRUE)

hist(w[,3] ,col="blue1", main="Histogram of test statistic for Wilcoxon test", xlab='Test Statistic', ylab='number of occurrence', labels=TRUE)
dev.off()

#### Problem-g
pdf("output-g.pdf")
hist(z[,2] ,col="blue1", main="Histogram of pvalue for two sample z-test", xlab='pvalue', ylab='number of occurrence', labels=TRUE)

hist(t[,2] ,col="blue1", main="Histogram of pvalue for two sample t-test", xlab='pvalue', ylab='number of occurrence', labels=TRUE)

hist(mw[,2] ,col="blue1", main="Histogram of pvalue for Mann Whitney U test", xlab='pvalue', ylab='number of occurrence', labels=TRUE)

hist(pz[,2] ,col="blue1", main="Histogram of pvalue for paired z-test", xlab='pvalue', ylab='number of occurrence', labels=TRUE)

hist(t[,2] ,col="blue1", main="Histogram of pvalue for paired t-test", xlab='pvalue', ylab='number of occurrence', labels=TRUE)

hist(w[,2] ,col="blue1", main="Histogram of pvalue for Wilcoxon test", xlab='pvalue', ylab='number of occurrence', labels=TRUE)
dev.off()

#### Problem-h
pdf("output-h.pdf")
plot(z[,2], pt[,2], cex=0.1, main="Z vs Paired T test", xlab="Z pvalue", ylab="Paired T pvalue")

plot(t[,2], pt[,2], cex=0.1, main="T vs Paired T test", xlab="T pvalue", ylab="Paired T pvalue")

plot(mw[,2], pt[,2], cex=0.1, main="Mann Whitney vs Paired T test", xlab="Mann Whitney pvalue", ylab="Paired T pvalue")

plot(pz[,2], pt[,2], cex=0.1, main="Paired Z Test vs Paired T test", xlab="Paired Z pvalue", ylab="Paired T pvalue")

plot(pt[,2], pt[,2], cex=0.1, main="Paired T test vs Paired T test", xlab="Paired T test pvalue", ylab="Paired T pvalue")

plot(w[,2], pt[,2], cex=0.1, main="Wilcoxon test vs Paired T test", xlab="Wilcoxon test pvalue", ylab="Paired T pvalue")
dev.off()



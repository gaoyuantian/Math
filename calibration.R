
install.packages("stabledist") 
library("stabledist")
# Empirical function
emp.fun = function(x){
	mu = mean(x)
	x = sort(x)
	n = length(x)
      
      q010 = quantile(x,0.1)
	q020 = quantile(x,0.2)
	q030 = quantile(x,0.3)
	q040 = quantile(x,0.4)
	q050 = quantile(x,0.5)
      q060 = quantile(x,0.6)
	q070 = quantile(x,0.7)
	q080 = quantile(x,0.8)
	q090 = quantile(x,0.9)


		
	g.x = c(q010,q020,q030,q040,q050,q060,q070,q080,q090)
	return(g.x)
}


# theoretical function
theo.fun = function(theta){
	alpha = theta[1]
	beta = theta[2]
	gamma = theta[3]
	delta = theta[4]
	g.theta = matrix(NA,h,9)
	for (i in 1:h){
		set.seed(1024 + i)
		x.star = rstable(n, alpha = alpha, beta = beta, gamma = gamma, delta = delta)
		g.theta[i,] = emp.fun(x.star)
	}
	g.return = apply(g.theta,2,mean)
	return(g.return)
}

# Objective function
obj.fun = function(theta){
	theo = theo.fun(theta)
	dif = theo - emp.estim
	obj = (t(dif))%*%Omega%*%dif
	return(obj)
}


# Globle variable


n <<- 500
h <<- 1
Omega <<- diag(c(1,1,1,1,1,1,1,1,1))
B = 100

theta0 = c(1.75,-0.2, 1, 0)


# create progress bar
pb <- txtProgressBar(min = 0, max = B, style = 3)
theta.hat1=matrix(NA,1000,4)
for (i in 1:B){
	set.seed(i)
	x = rstable(n, alpha = theta0[1], beta = theta0[2], gamma = theta0[3], delta = theta0[4])
	emp.estim <<- emp.fun(x)
        theta.start = theta0
	Omega <<- diag(c(1,1,1,1,1,1,1,1,1))
	theta.hat1[i,] = optim(theta.start,obj.fun, method = "L-BFGS-B", lower = c(1.00001,-0.9999,0.00001,-100), upper = c(1.9999,0.9999,100,100))$par
	
	
	par(mfrow = c(2,2))
	boxplot(theta.hat1[,1],col = "lightgrey", main = expression(alpha), ylim = c(1.25,1.99))
      abline(h = theta0[1])
	boxplot(theta.hat1[,2], col = "lightgrey", main = expression(beta), ylim = c(-0.5,0))
      abline(h = theta0[2])
      boxplot(theta.hat1[,3], col = "lightgrey", main = expression(gamma), ylim = c(0.5,1.5))
      abline(h = theta0[3])
      boxplot(theta.hat1[,4], col = "lightgrey", main = expression(delta), ylim = c(-0.2,0.2))
      abline(h = theta0[4])
	setTxtProgressBar(pb, i)
    
}

close(pb)


### omega change to inverse of covariance matrix of quanitle estimator

boot.Var = function(theta, B = 100){
	alpha = theta[1]
	beta = theta[2]
      gamma= theta[3]
      delta= theta[4]
	emp.boot = matrix(NA,B,9)
	for (i in 1:B){
		set.seed(i)
		x.star = rstable(n, alpha = alpha, beta = beta, gamma = gamma, delta = delta)
		emp.boot[i,]= emp.fun(x.star)
	}
	return(cov(emp.boot))
}

V = boot.Var(theta.hat1)
omega=solve(V)




# create progress bar

pb <- txtProgressBar(min = 0, max = B, style = 3)
theta.hat2=matrix(NA,B=1000,4)
for (i in 1:B){
	set.seed(i)
	x = rstable(n, alpha = theta0[1], beta = theta0[2], gamma = theta0[3], delta = theta0[4])
	emp.estim <<- emp.fun(x)
        theta.start = theta0
        v=boot.var(theta.hat1[i,])
	Omega <<-solve(v)
	theta.hat2[i,] = optim(theta.start,obj.fun, method = "L-BFGS-B", lower = c(1.00001,-0.9999,0.00001,-100), upper = c(1.9999,0.9999,100,100))$par
	
	
	par(mfrow = c(2,2))
	boxplot(theta.hat2[,1],col = "lightgrey", main = expression(alpha), ylim = c(1.35,1.7))
      abline(h = theta0[1])
	boxplot(theta.hat2[,2], col = "lightgrey", main = expression(beta), ylim = c(-0.7,0))
      abline(h = theta0[2])
      boxplot(theta.hat2[,3], col = "lightgrey", main = expression(gamma), ylim = c(2.1,2.8))
      abline(h = theta0[3])
      boxplot(theta.hat2[,4], col = "lightgrey", main = expression(delta), ylim = c(-0.1,2))
      abline(h = theta0[4])
	setTxtProgressBar(pb, i)
}

close(pb)






#STA457 Assn 1 Part A
#Ahmed N Amanullah; 1001325773

#Step 1: Get the variance of the predictor
# d: vector of dj coefficients (j = 0, ..., m-2)
# X: log returns
#Will use quadratic form... 
varF <- function(d, X){
	#maximum lag
	M <- length(d)-1
	#get auto-covariance
	acfs <- acf(X, plot = F, type="covariance", lag.max=M)$acf
	#get toeplitz matrix of acfs
	#name used to refer to:
	#[[\gamma_0, ..., \gamma_{m-2}],
	# [\gamma_1, \gamma_0, ..., \gamma_{m-3}], ...
	# [\gamma_{m-2}, ..., \gamma_0]]
	Gamma <- toeplitz(as.vector(acfs))
	#quadratic form
	varF <- d%*%Gamma%*%as.vector(d)
	varF
}

#Step 2: Get the expectation of the predictor
# d: vector of dj coefficients (j = 0, ..., m-2)
# X: log returns
muF <- function(d, X){
	muF <- mean(X)*sum(d)
	muF
}

#Step 3: Computing ACF(1) of forecaster
#will use a quadratic form
# d: vector of dj coefficients (j = 0, ..., m-2)
# X: log returns
rhoF <- function(d, X){
	#max lag
	M <- length(d)-1
	#get acfs
	acfs <- acf(X, plot = F, type = "covariance", lag.max = M+2)$acf
	#get the covariance
	#print(acfs)
	#print(M)
	#print(matrix(acfs[abs(outer(0:M,1:(M+1), "-"))+1], 
	#M+1, M+1))
	temp <- d%*%matrix(acfs[as.vector(
	abs(outer(0:M,1:(M+1), "-"))+1)], M+1, M+1)%*%as.vector(d)
	rhoF <- temp/varF(d, X)
	rhoF
}

#Step 4: expected rule returns and holding period
#first correlation between X_t and F_{t-1}
corXF <- function(d, X){
	Mp <- length(d)
	acfs <- acf(X, plot = F, type = "covariance", lag.max = Mp)$acf
	corXF <- sum(d*acfs[-1])/sqrt(acfs[1]*varF(d, X))
	corXF
}

#expected holding period
#rho is lag one autocorrelation
#provided by step 3
Hold <- function(rho){
	Hold <- pi/acos(rho)
}

#next is the double MA co-efficients
d <- function(m, r){
	d <- c((m-r)*((0:(r-1))+1), r*(m-(r: (m-1))-1))
	d
}

# retX: log asset return
# m: long-term MA
# r: short-term MA

ruleReturn <- function(retX, m, r){
	#returns variance
	vX <- sd(retX)
	#returns mean
	mX <- mean(retX)
	#predictor mean
	mF <- muF(d(m, r), retX)
	#predictor sd
	vF <- sqrt(varF(d(m, r), retX))
	#lag 1 correlation (rho) between X_t and F_{t-1}
	rXF <- corXF(d(m, r), retX)
	#lag 1 correlation of predictor
	rF <- rhoF(d(m, r), retX)
	#expected return
	ER <- sqrt(2/pi)*vX*rXF*exp(-mF*mF/(vF*vF)) + 
	mX*(1-2*pnorm(-mF/vF))
	#expected hold duration
	H <- Hold(rF)
	#return results
	list("ER" = ER, "H" = H, "rhoF" = rF, "VF"= vF, "muF"=mF, 
	"corXF"=rXF)
}

#This function is to download data from folder...
#idea: return three arrays - one containing the stock tickers...
#	another the adjusted close (using hint from prev assn) 
#	and finally one containing the log returns

# dataDir: directory containing csv files containing data
# currDir: current working directory (to change back to after 
#	loading data)
getAdC <- function(dataDir, currDir){
	#switch to directory containing price data
	setwd(dataDir)
	#get list of files in directory
	priceFiles <- list.files(dataDir, pattern="*.csv", 
	full.names = F, recursive = F)
	#to hold tickers
	tickers <- c()
	#to hold prices
	prices <- c()
	#to hold log returns
	logrets <- c()
	#iterate through list above
	for (file in seq(1, length(priceFiles))){
		tickers <- c(tickers, substr(priceFiles[file],
		1,nchar(priceFiles[file])-4))
		#use rbind for the adjusted close data...
		#first read the file
		stockdat <- read.csv(file = priceFiles[file],
		 header = T, sep=",")
		#then get adjusted close prices from file
		price <- stockdat[1:dim(stockdat)[1], 6] 
		if (length(prices) < 1){
			prices <- price
		}else{
			prices <- rbind(prices, price)
		}
		#also the log returns
		if (length(logrets) < 1){
			logrets <- log(price[2:length(price)]/
			price[1:(length(price)-1)])
		}else{
			logrets <- rbind(logrets, log(price[2:length(price)]/
			price[1:(length(price)-1)]))
		}
	}
	#return to original directory
	setwd(currDir)
	
	#return tickers, prices and log returns
	list("tickers"=tickers, "prices"=prices, "logrets" = logrets)
}

#Function to choose the optimal daily and monthly MA trading rules
#	(i.e. that maximize expected rule returns)
#will test function against the one in prev assn after 
#	implementation...

#first optimal monthly
#	retX: vector of log returns
monthlyoptimal_dma <- function(retX){
	#to hold optimal m and r
	optimal_m <- 2
	optimal_r <- 1
	#get ruleReturn for this setting
	currER <- ruleReturn(retX, optimal_m, optimal_r)$ER
	#iterate up to 12 max as we are doing monthly
	#loop through r
	for (i in seq(1, 11)){
		for (j in seq(i+1, 12)){
			ERij <- ruleReturn(retX, j, i)$ER
			if (ERij > currER){
				optimal_m <- j
				optimal_r <- i
				currER <- ERij
			} 
		}
	}
	#return optimal double MA trading rules
	list("monthlyoptimal_m"=optimal_m, "monthlyoptimal_r"=optimal_r)
}

#next optimal daily - will try only m and r values as in prev
#	assn due to efficiency concerns

#	retX: vector of log returns
dailyoptimal_dma <- function(retX){
	#to hold optimal m and r
	optimal_m <- 2
	optimal_r <- 1
	#get ruleReturn for this setting
	currER <- ruleReturn(retX, optimal_m, optimal_r)$ER
	#iterate up to 12 max as we are doing monthly
	#loop through r
	for (i in c(1,5,10,20,60,249)){
		for (j in c(5, 10, 20, 60, 120, 250)){
			if (j > i){
				ERij <- ruleReturn(retX, j, i)$ER
				if (ERij > currER){
					optimal_m <- j
					optimal_r <- i
					currER <- ERij
				} 	
			}
		}
	}
	#return optimal double MA trading rules
	list("dailyoptimal_m"=optimal_m, "dailyoptimal_r"=optimal_r)
}

#functions that will implement the above for all constituents...

#	retsX: the log returns of the constituents 

monthlyoptimals_dma <- function(retsX){
	#to hold the optimal trading rules
	m <- c()
	r <- c()
	#amount of tickers to go through
	numstocks <- dim(retsX)[1]
	#number of periods/ months
	months <- dim(retsX)[2]
	for (i in seq(1, numstocks)){
		optimals <- monthlyoptimal_dma(retsX[i, 1:months])
		#add optimals to list
		m <- c(m, optimals$monthlyoptimal_m)
		r <- c(r, optimals$monthlyoptimal_r)
	} 
	
	#return optimal rules
	list("monthlyoptimals_m"=m, "monthlyoptimals_r"=r)
}

#the version for the getting the optimal rules for each window
windowsmonthlyoptimals_dma <- function(retX){
	#number of stock
	num_stocks <- dim(retX)[1]
	#the number of periods
	months <- dim(retX)[2]
	#hold optimal m's and r's
	m <- matrix(,nrow = num_stocks, ncol = 1)
	r <- matrix(,nrow = num_stocks, ncol = 1)
	for (i in seq(1, months - 60, 12)){
		#get optimal dma within 60 month window for all stocks
		optimal_dma <- monthlyoptimals_dma(
		retX[1:num_stocks, i:(i+60-1)])
		if(i == 1){
			m[1:num_stocks, 1] <- optimal_dma$monthlyoptimals_m
			r[1:num_stocks, 1] <- optimal_dma$monthlyoptimals_r	
		}else{
			m <- cbind(m, optimal_dma$monthlyoptimals_m)
			r <- cbind(r, optimal_dma$monthlyoptimals_r)
		}
	}
	
	list("m"=m, "r"=r)
}

#functions for computing in-sample trading statistics 
#	(i.e. cumulative return and holding time)
#	to later compare with theoretical results
#first function to compute the predictor

#	d: d is d(m ,r) - presumably the optimal ones found
#	retX: the log returns
#	t: get signal at time t

f <- function(d, retX, t){
	M <- length(d) - 1
	if (t >= M){
		output <- sum(d*rev(retX[(t-M):t]))
		output
	}
	else {print('t is smaller than M')}
}

#next the realized return
#	d: d is d(m ,r) - presumably the optimal ones found
#	retX: the log returns
realized_ret <- function(d, retX){
	re <- c()
	M <- length(d)
	for (t in seq(M + 1,length(retX))){
		re <- c(re, sign(f(d, retX, t-1))*retX[t])
	}
	realized_ret <- sum(re)/length(M:length(retX))
	realized_ret
}

#next is the holding period
holding_period <- function(d, retX){
	num_change <- 0
	M <- length(d)
	for (t in seq(M+1, length(retX))){
		if (sign(f(d, retX,t)) != sign(f(d, retX,t-1))){
			num_change <- num_change + 1
		}
	}
	hold <- length(M:length(retX))/num_change
	hold
}

#function that will take in the array of log returns 
#	and the optimal trading decisions
in_sample_estimate <- function(retX, m, r){
	#number of tickers
	num_stocks <- dim(retX)[1]
	#number of periods/ months
	#to hold in_sample and theoretical values
	djER <- c()
	djH <- c()
	djISR <- c()
	djISH <- c()
	months <- dim(retX)[2]
	for (i in seq(1, num_stocks)){
		ruleReturn_i <- ruleReturn(retX[i, 1:months], m[i], r[i])
		djER <- c(djER, ruleReturn_i$ER)
		djH <- c(djH, ruleReturn_i$H)
		djISR <- c(djISR, 
		realized_ret(d(m[i],r[i]),retX[i, 1:months]))
		djISH <- c(djISH, 
		holding_period(d(m[i],r[i]),retX[i, 1:months]))
	}
	#return estimated and theoretical values 
	list("djER"=djER ,"djH"=djH, "djISR"=djISR, "djISH"=djISH)
}

#Part B Question 1

#	retX: provide the daily log returns only up to the time t
#				(i.e. the time at which the volatility is to be
#				 computed)

ea_volatility <- function(retX){
	#need to solve for delta
	#according to paper delta has been solved for as
	delta <- 0.2
	#number of stocks
	num_stocks <- dim(retX)[1]
	#number of periods provided
	t <- dim(retX)[2]
	#vector to hold the ex ante volatilities
	ea_sd <- c()
	#need to iterate over the stocks
	for (s in seq(1, num_stocks)){
		#compute the exponentially weighted returns
		r_bar <- sum((1-delta)*(delta^c(0:11))*
		rev(retX[s, (t-12):(t-1)]))
		#next the ex-ante volatility (annualized)
		#first the variance
		var <- 12*sum((1-delta)*(delta^c(0:11))*
		(rev(retX[s, (t-12):(t-1)])-r_bar)^2)
		#then add the ex-ante volatility of stock s to list
		ea_sd <- c(ea_sd, sqrt(var))
	}
	
	#return volatilities
	ea_sd
} 


#Question 2: summarizing performance of the equally weighted
#	portfolio
#	retX: 30X(periods) dim matrix of monthly log returns

ew_performance <- function(retX){
	#number of stocks provided
	num_stocks <- dim(retX)[1]
	#number of periods of data provided
	months <- dim(retX)[2]
	#to hold the ruled returns
	re <- c()
	#to hold the cumulative return up to time t
	cre <- c()
	#to hold expected return at time t
	ere <- c()
	#to hold the ex-ante volatility at time t
	ea_sd <- c()
	#to hold sharpe ratio at time t
	sharpes <- c()
	#only go up to as long as window can be fit
	for (i in seq(1, months - 60, 12)){
		print(i)
		#get optimal dma within 60 month window for all stocks
		optimal_dma <- monthlyoptimals_dma(
		retX[1:num_stocks, i:(i+60-1)])
		m <- optimal_dma$monthlyoptimals_m
		r <- optimal_dma$monthlyoptimals_r
		#need to get the ex-ante volatility every 12 months
		ea_sd_i <- ea_volatility(
		retX[1:num_stocks,1:(i+60-1)])
		#Ok, now we use the optimal trading rules to find the 
		#	ruled return for each stock
		#to hold portfolio weights
		#this is an ew portfolio so all weights are the same
		weights <- rep(1/num_stocks, num_stocks)
		
		#this is the portfolio variance
		var_t <- 0
		for (s in seq(1, num_stocks)){
			var_t <- var_t + (weights[s]^2)*ea_sd_i[s]^2
		}
		
		ea_sd <- c(ea_sd, sqrt(var_t))
		
		#going back to old method
		for (t in seq(i + 60, i + 72 -1)){
			if (t <= months){
				re_t <- 0
				for (s in seq(1, num_stocks)){
					s_d <- d(m[s], r[s])
					re_t <- re_t + weights[s]*
					sign(f(s_d, retX[s, 1:months], t-1))*retX[s,t]
				}
				#now add the portfolio ruled return
				re <- c(re, re_t)
				if (length(re) > 1){
					cre <- c(cre, cre[length(cre)] + re_t)
				}else{
					cre <- c(cre, re_t)
				}
				#expected (annual) return at time t
				ere <- c(ere, (12*cre[length(cre)]/length(cre)))
				#as well the sharpe
				sharpe_t <- (ere[length(ere)]-0.02)/
				ea_sd[length(ea_sd)]
				sharpes <- c(sharpes, sharpe_t)
				#note that since this is a EW portfolio, there is 
				#	no rebalancing step as weights remain the same 
			}
		}
		
		
	}
	
	#return results
	list("return"=ere[length(ere)], 
	"volatility"=ea_sd[length(ea_sd)], 
	"sharpe"=sharpes[length(sharpes)], "cumul_return"=cre, 
	"e_return"=ere, "ea_volatilities"=ea_sd, "sharpes"=sharpes)
}

#risk parity function
rp_performance <- function(retX){
	#number of stocks provided
	num_stocks <- dim(retX)[1]
	#number of periods of data provided
	months <- dim(retX)[2]
	#to hold the ruled returns
	re <- c()
	#to hold the cumulative return up to time t
	cre <- c()
	#to hold expected return at time t
	ere <- c()
	#to hold the ex-ante volatility at time t
	ea_sd <- c()
	#to hold sharpe ratio at time t
	sharpes <- c()
	#only go up to as long as window can be fit
	for (i in seq(1, months - 60, 12)){
		print(i)
		#get optimal dma within 60 month window for all stocks
		optimal_dma <- monthlyoptimals_dma(
		retX[1:num_stocks, i:(i+60-1)])
		m <- optimal_dma$monthlyoptimals_m
		r <- optimal_dma$monthlyoptimals_r
		#need to get the ex-ante volatility every 12 months
		ea_sd_i <- ea_volatility(
		retX[1:num_stocks,1:(i+60-1)])
		#Ok, now we use the optimal trading rules to find the 
		#	ruled return for each stock
		#to hold portfolio weights
		#this is an ew portfolio so all weights are the same
		temp <- 1/ea_sd_i
		weights <- temp/sum(temp)
		
		#this is the portfolio variance
		var_t <- 0
		for (s in seq(1, num_stocks)){
			var_t <- var_t + (weights[s]^2)*ea_sd_i[s]^2
		}
		
		ea_sd <- c(ea_sd, sqrt(var_t))
		
		#going back to old method
		for (t in seq(i + 60, i + 72 -1)){
			if (t <= months){
				re_t <- 0
				for (s in seq(1, num_stocks)){
					s_d <- d(m[s], r[s])
					re_t <- re_t + weights[s]*
					sign(f(s_d, retX[s, 1:months], t-1))*retX[s,t]
				}
				#now add the portfolio ruled return
				re <- c(re, re_t)
				if (length(re) > 1){
					cre <- c(cre, cre[length(cre)] + re_t)
				}else{
					cre <- c(cre, re_t)
				}
				#expected (annual) return at time t
				ere <- c(ere, (12*cre[length(cre)]/length(cre)))
				#as well the sharpe
				sharpe_t <- (ere[length(ere)]-0.02)/
				ea_sd[length(ea_sd)]
				sharpes <- c(sharpes, sharpe_t)
				#note that since this is a EW portfolio, there is 
				#	no rebalancing step as weights remain the same 
			}
		}
		
		
	}
	
	#return results
	list("return"=ere[length(ere)], 
	"volatility"=ea_sd[length(ea_sd)], 
	"sharpe"=sharpes[length(sharpes)], "cumul_return"=cre, 
	"e_return"=ere, "ea_volatilities"=ea_sd, "sharpes"=sharpes)
}

#Junk code:
'''
		for (t in seq(i + 60, i + 72 -1)){
			if (t <= months){
				re_t <- 0
				for (s in seq(1, num_stocks)){
					s_d <- d(m[s], r[s])
					re_t <- re_t + (1/num_stocks)*
					sign(f(s_d, retX[s, 1:months], t-1))*retX[s,t]
				}
				#now add the portfolio ruled return
				re <- c(re, re_t)
				#note that since this is a EW portfolio, there is 
				#	no rebalancing step as weights remain the same 
			}
		}
'''

#Part C
#Ques 1: compute the expected h-period holding period return

#	retX: the monthly log returns
#	m and r: the double MA parameters
#	h: the holding periods

e_hperiod_holdperiodreturn <- function(retX, m, r, h){
	#get auto-covariances
	acfs <- acf(retX, plot = F, type = "covariance", lag.max = m)$acf
	#get expected log return
	ER <- mean(retX)
	#get expected holding period return
	ehphpr <- h*sum(d(m,r)*(acfs[2:length(acfs)] + (ER^2))) 
	#return value
	ehphpr
}

#Ques 2: optimal double MA for 30 DJ constituents maximizing 
#	12-period holding period return.
monthlyoptimalEHR_dma <- function(retX){
	#to hold optimal m and r
	optimal_m <- 2
	optimal_r <- 1
	#get ruleReturn for this setting
	currEHR <- e_hperiod_holdperiodreturn(retX, 
	optimal_m, optimal_r, 12)
	#iterate up to 12 max as we are doing monthly
	#loop through r
	for (i in seq(1, 11)){
		for (j in seq(i+1, 12)){
			EHRij <- e_hperiod_holdperiodreturn(retX, 
			j, i, 12)
			if (EHRij > currEHR){
				optimal_m <- j
				optimal_r <- i
				currEHR <- EHRij
			} 
		}
	}
	#return optimal double MA trading rules
	list("monthlyoptimal_m"=optimal_m, "monthlyoptimal_r"=optimal_r)
}

monthlyoptimalsEHR_dma <- function(retsX){
	#to hold the optimal trading rules
	m <- c()
	r <- c()
	#amount of tickers to go through
	numstocks <- dim(retsX)[1]
	#number of periods/ months
	months <- dim(retsX)[2]
	for (i in seq(1, numstocks)){
		optimals <- monthlyoptimalEHR_dma(retsX[i, 1:months])
		#add optimals to list
		m <- c(m, optimals$monthlyoptimal_m)
		r <- c(r, optimals$monthlyoptimal_r)
	} 
	
	#return optimal rules
	list("monthlyoptimals_m"=m, "monthlyoptimals_r"=r)
}

#the version for the getting the optimal rules for each window
windowsmonthlyoptimalsEHR_dma <- function(retX){
	#number of stock
	num_stocks <- dim(retX)[1]
	#the number of periods
	months <- dim(retX)[2]
	#hold optimal m's and r's
	m <- matrix(,nrow = num_stocks, ncol = 1)
	r <- matrix(,nrow = num_stocks, ncol = 1)
	for (i in seq(1, months - 60, 12)){
		#get optimal dma within 60 month window for all stocks
		optimal_dma <- monthlyoptimalsEHR_dma(
		retX[1:num_stocks, i:(i+60-1)])
		if(i == 1){
			m[1:num_stocks, 1] <- optimal_dma$monthlyoptimals_m
			r[1:num_stocks, 1] <- optimal_dma$monthlyoptimals_r	
		}else{
			m <- cbind(m, optimal_dma$monthlyoptimals_m)
			r <- cbind(r, optimal_dma$monthlyoptimals_r)
		}
	}
	
	list("m"=m, "r"=r)
}


#--------------------------------------------------------------
# Test functions here...
#--------------------------------------------------------------
#need to test function,... 
#need to download the csv file data
#code for it:
#change work directory for reading file...
workDirect <- "/Users/amanull9/Documents/STA457/Assignments/STA457-Assn1/Assn 1/DJ Data Monthly"
setwd(workDirect)
#get the data for any one of the constituents - using KO for now
KOdat <- read.csv(file = "KO.csv", header = T, sep=",")
#change back to current directory
workDirect <- "/Users/amanull9/Documents/STA457/Assignments/STA457-Assn1/Assn 1"
setwd(workDirect)

#test getAdC function
dataDir <- "/Users/amanull9/Documents/STA457/Assignments/STA457-Assn1/Assn 1/DJ Data Monthly"
currDir <- "/Users/amanull9/Documents/STA457/Assignments/STA457-Assn1/Assn 1"

monthlydata <- getAdC(dataDir, currDir)
tickers <- monthlydata$tickers
mprices <- monthlydata$prices
mlogrets <- monthlydata$logrets

#aside from getting warning messages, output is as expected to we
# can proceed to part 6

#test monthly dma
KOdma <- monthlyoptimal_dma(mlogrets[15, 1:dim(mlogrets)[2]])
#Seems to be outputting something - dunno if it is correct...
#implement the daily equivalent and see if results match prev assn

#test daily dma
#load SnP data
SnPdat <- read.csv(file = "^GSPC.csv", header = T, sep=",")
#get log daily returns
SnP_AdC <- SnPdat[1:dim(SnPdat)[1], 6]
SnP_logrets <- log(SnP_AdC[2:length(SnP_AdC)]/
SnP_AdC[1:(length(SnP_AdC)-1)])

#test function
SnPdma <- dailyoptimal_dma(SnP_logrets)
#matches... so some sanity check here

#So now will go ahead and finish up ques 1
#get optimals for all 30 constituents
dj_optimals <- monthlyoptimals_dma(mlogrets)

#construct table
djoptimals_table <- t(rbind(tickers, dj_optimals$monthlyoptimals_m,dj_optimals$monthlyoptimals_r))

#column names
colnames(djoptimals_table)[2:3] <- c("m", "r")

#test the windows version
windowsdj_optimals <- windowsmonthlyoptimals_dma(mlogrets)

#construct table
optimal_ms <- cbind(tickers, windowsdj_optimals$m)

#column names
colnames(optimal_ms)[2:dim(optimal_ms)[2]] <- seq(1999, 1999+dim(optimal_ms)[2]-2, 1)

#construct table
optimal_rs <- cbind(tickers, windowsdj_optimals$r)

#column names
colnames(optimal_rs)[2:dim(optimal_rs)[2]] <- seq(1999, 1999+dim(optimal_ms)[2]-2, 1)

#Ok ques 1 done, have optimal dma trading rules for all 30 
#	constituents above!

#test the in sample function next...
#first test the SnP as sanity...
SnP_ISR <- realized_ret(d(SnPdma$dailyoptimal_m, SnPdma$dailyoptimal_r), SnP_logrets)
#this matches with 0.00024 from prev assn

SnP_ISH <- holding_period(d(SnPdma$dailyoptimal_m, SnPdma$dailyoptimal_r), SnP_logrets)
#this also matches with 40.3 from prev assn 

#test the ER stuff
SnP_ruleret <- ruleReturn(SnP_logrets,SnPdma$dailyoptimal_m, SnPdma$dailyoptimal_r)
#theoretical H of 0.00159 is good
#theoretical ret is off (prev assn is 0.00332 and I have 0.000345)
#so are they using something different - ask?

#have test for insample ready
dj_tradestats <- in_sample_estimate(mlogrets,dj_optimals$monthlyoptimals_m,dj_optimals$monthlyoptimals_r)

#Ok outputting something which is good

#code for having this in a table
dj_tradestats_table <- t(rbind(tickers, round(dj_tradestats$djER,4),  round(dj_tradestats$djISR,4), round(dj_tradestats$djH,3), round(dj_tradestats$djISH,3)))

colnames(dj_tradestats_table)[2:5] <- c("theoretical return",  "expected return","theoretical hold", "expected hold")


#Now test the ewperformances
#test ew function
dj_ew_performance <- ew_performance(mlogrets)

#code for having this in a table
dj_ew_performance_table <- t(rbind(round(dj_ew_performance$return,4),  round(dj_ew_performance$volatility,4), round(dj_ew_performance$sharpe,4)))

colnames(dj_ew_performance_table)[1:3] <- c("Annualized expected return",  "Annualized volatility","Annualized sharpe ratio")

#kinda slow but at least have an output...
#also performance does not seem good, am I doing something wrong?

#hmm - check for risk parity and see if it is still bad
dj_rp_performance <- rp_performance(mlogrets)

#code for having this in a table
dj_rp_performance_table <- t(rbind(round(dj_rp_performance$return,4),  round(dj_rp_performance$volatility,4), round(dj_rp_performance$sharpe,4)))

colnames(dj_rp_performance_table)[1:3] <- c("Annualized expected return",  "Annualized volatility","Annualized sharpe ratio")

#code for plots...
par(mfrow = c(2,2))
plot(c(61:dim(mlogrets)[2]), dj_ew_performance$cumul_return, main = "cumulative log returns", type ="l", xlab = "month", ylab="cumulative return", ylim = c(-0.10, 0.30))
lines(c(61:dim(mlogrets)[2]), dj_rp_performance$cumul_return, col="blue")
plot(c(61:dim(mlogrets)[2]), dj_ew_performance$e_return, main = "annualized expected returns", type = 'l', xlab = "month", ylab="annualized expected return")
lines(c(61:dim(mlogrets)[2]), dj_rp_performance$e_return, col="blue")
plot(seq(61,dim(mlogrets)[2],12), dj_ew_performance$ea_volatilities, main = "annualized ex ante volatilities", type = 'l', xlab = "month", ylab="annualized ex ante volatility", ylim = c(0, 0.04))
lines(seq(61,dim(mlogrets)[2],12), dj_rp_performance$ea_volatilities, col="blue")
plot(c(61:dim(mlogrets)[2]), dj_ew_performance$sharpes, main = "sharpe ratio", type = 'l', xlab = "month", ylab="sharpe ratio", ylim = c(-7, 6))
lines(c(61:dim(mlogrets)[2]), dj_rp_performance$sharpes, col="blue")

#test expected h-period holding period ret:
ehphpr_test <- e_hperiod_holdperiodreturn(mlogrets[1, 1:dim(mlogrets)[2]], 3, 2, 12)

#eh, seems to be outputting something so ok...
#So now will go ahead and finish up ques 3
#get optimals for all 30 constituents
dj_optimalsEHR <- monthlyoptimalsEHR_dma(mlogrets)

#construct table
djoptimalsEHR_table <- t(rbind(tickers, dj_optimalsEHR$monthlyoptimals_m,dj_optimalsEHR$monthlyoptimals_r))

#column names
colnames(djoptimalsEHR_table)[2:3] <- c("m", "r")

#test the windows version
windowsdj_optimals <- windowsmonthlyoptimalsEHR_dma(mlogrets)

#construct table
optimal_ms <- cbind(tickers, windowsdj_optimals$m)

#column names
colnames(optimal_ms)[2:dim(optimal_ms)[2]] <- seq(1999, 1999+dim(optimal_ms)[2]-2, 1)

#construct table
optimal_rs <- cbind(tickers, windowsdj_optimals$r)

#column names
colnames(optimal_rs)[2:dim(optimal_rs)[2]] <- seq(1999, 1999+dim(optimal_ms)[2]-2, 1)
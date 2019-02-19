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
		prices <- rbind(prices, price)
		#also the log returns
		logrets <- rbind(logrets, 
		log(price[2:length(price)]/price[1:(length(price)-1)]))
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
ea_volatility <- function(retX_daily){
	#need to solve for delta
	
} 


#Question 2: summarizing performance of the equally weighted
#	portfolio
#	retX: 30X(periods) dim matrix of log returns

ew_performance <- function(retX){
	#number of stocks provided
	num_stocks <- dim(retX)[1]
	#number of periods of data provided
	months <- dim(retX)[2]
	#to hold the ruled returns
	re <- c()
	#print(num_stocks)
	#only go up to as long as window can be fit
	#print(months)
	#print(seq(1, months - 60, 12))
	for (i in seq(1, months - 60, 12)){
		#get optimal dma within 60 month window for all stocks
		print(i)
		optimal_dma <- monthlyoptimals_dma(
		retX[1:num_stocks, i:(i+60-1)])
		m <- optimal_dma$monthlyoptimals_m
		r <- optimal_dma$monthlyoptimals_r
		#print(i)
		#Ok, now use this to get ruled returns for the 12 months
		#	after this rolling window
		#we will iterate through the 12 months following this window
		#	and iterate through the stocks
		for (t in seq(i + 60, i + 72 -1)){
			if (t <= months){
				re_t <- 0
				#print(num_stocks)
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
	}
	#the above will give me the ruled returns - can use this to
	#	get realized return...
	#note that annualized expected return, volatility and sharpe
	#	ratio are needed
	#we can get the monthly expected return as follows:
	realized_ret <- sum(re)/length(re)
	#can convert this to annual as follows:
	realized_ret <- 12*realized_ret
	#question is though how I get the annualized volatility?
	#do I just get the volatility of the assets from the returns
	#	data and then just sum (1/30)^2*volatility, then annualize?
	#do this for now...
	volatility <- 0
	for (s in seq(1, num_stocks)){
		volatility <- volatility + (1/30)^2*(sd(retX[s, 1:months]))^2
	}
	#convert to annualized
	volatility <- sqrt(12)*sqrt(volatility)
	#next is annualized sharpe ratio
	#assuming 0.02 annual-risk free for this assn
	sharpe <- (realized_ret - 0.02)/volatility
	#return results
	list("return"=realized_ret, "volatility"=volatility,
	 "sharpe"=sharpe)
}


#--------------------------------------------------------------
# Test functions here...
#--------------------------------------------------------------
#need to test function,... 
#need to download the csv file data
#code for it:
#change work directory for reading file...
workDirect <- "/Users/amanull9/Documents/STA457/Assignments/Assn 1/DJ Data Monthly"
setwd(workDirect)
#get the data for any one of the constituents - using KO for now
KOdat <- read.csv(file = "KO.csv", header = T, sep=",")
#change back to current directory
workDirect <- "/Users/amanull9/Documents/STA457/Assignments/Assn 1"
setwd(workDirect)

#test getAdC function
dataDir <- "/Users/amanull9/Documents/STA457/Assignments/Assn 1/DJ Data Monthly"
currDir <- "/Users/amanull9/Documents/STA457/Assignments/Assn 1"

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

#test ew function
dj_ew_performance <- ew_performance(mlogrets)

#code for having this in a table
dj_ew_performance_table <- t(rbind(round(dj_ew_performance$return,4),  round(dj_ew_performance$volatility,4), round(dj_ew_performance$sharpe,4)))

colnames(dj_ew_performance_table)[1:3] <- c("Annualized expected return",  "Annualized volatility","Annualized sharpe ratio")

#kinda slow but at least have an output...

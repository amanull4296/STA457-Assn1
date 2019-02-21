# B
# question 2

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

rhoF <- function(d, X){
  #max lag
  M <- length(d)-1
  #get acfs
  acfs <- acf(X, plot = F, type = "covariance", lag.max = M+2)$acf
  #get the covariance
  temp <- d%*%matrix(acfs[as.vector(
    abs(outer(0:M,1:(M+1), "-"))+1)], M+1, M+1)%*%as.vector(d)
  rhoF <- temp/varF(d, X)
  rhoF
}
#Step 4 Intermediate: correlation between X_t and F_{t-1}
#will use a quadratic form
# d: vector of dj coefficients (j = 0, ..., m-2)
# X: log returns

corXF <- function(d, X){
  Mp <- length(d)
  acfs <- acf(X, plot = F, type = "covariance", lag.max = Mp)$acf
  corXF <- sum(d*acfs[-1])/sqrt(acfs[1]*varF(d, X))
  corXF
}
Hold <- function(rho){
  Hold <- pi/acos(rho)
}
muF <- function(d, X){
  muF <- mean(X)*sum(d)
  muF
}

#Step 4: Expected holding period and ruled returns
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

ea_volatility <- function(retX_daily){
  #need to solve for delta
  #according to paper delta has been solved for as
  delta <- 60/61
  #number of stocks
  num_stocks <- dim(retX_daily)[1]
  #number of periods provided
  t <- dim(retX_daily)[2]
  #vector to hold the ex ante volatilities
  ea_sd <- c()
  #need to iterate over the stocks
  for (s in seq(1, num_stocks)){
    #compute the exponentially weighted returns
    r_bar <- sum((1-delta)*(delta^c(0:260))*
                   rev(retX_daily[s, (t-261):(t-1)]))
    #next the ex-ante volatility (annualized)
    #first the variance
    var <- 261*sum((1-delta)*(delta^c(0:260))*
                     (rev(retX_daily[s, (t-261):(t-1)])-r_bar)^2)
    #then add the ex-ante volatility of stock s to list
    ea_sd <- c(ea_sd, sqrt(var))
  }
  
  #return volatilities
  ea_sd
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
# getBestRegression <- function(retX){
#   num_stocks <- dim(retX)[1]
#   t <- dim(retX)[2]
#   best_h_value <- 0
#   
#   for (s in seq(1, num_stocks)){
#     for (h in c(1:228)) {
#       linearMod <- lm()
#     }
#   }
# }

tsMom <- function(retX, retX_daily){
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
    #get optimal dma within 60 month window for all stocks
    optimal_dma <- monthlyoptimals_dma(
      retX[1:num_stocks, i:(i+60-1)])
    m <- optimal_dma$monthlyoptimals_m
    r <- optimal_dma$monthlyoptimals_r
    #need to get the ex-ante volatility every 12 months
    #since primarily using US stocks, will go by 252 trading
    #	days a year - 21 days per month
    ea_sd_i <- ea_volatility(
      retX_daily[1:num_stocks,1:((i+60-1)*21)])
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

dataDir <- "~/GitHub/STA457-Assn1/Assn 1/DJ Data Daily"
currDir <- "~/GitHub/STA457-Assn1/Assn 1"

monthlydata <- getAdC(dataDir, currDir)
tickers <- monthlydata$tickers
mprices <- monthlydata$prices
mlogrets <- monthlydata$logrets
dailydata <- getAdC(dataDir, currDir)
dtickers <- dailydata$tickers
dprices <- dailydata$prices
dlogrets <- dailydata$logrets

eavolatilities <- ea_volatility(dlogrets);
# bestHvalues <- getBestRegression(mlogrets)
val <- tsMom(mlogrets, dlogrets)
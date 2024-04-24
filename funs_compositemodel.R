# A function that returns probability density function of Pareto random variables
dpareto <- function (x, alpha_, theta_) {
    res <- alpha_*theta_^alpha_/(x+theta_)^(alpha_+1)
    return(res)
}

# A function that returns cumulative distribution function of Pareto random variables
ppareto <- function (x, alpha_, theta_, lower.tail = TRUE) {
    if (lower.tail) {
        res <- 1-(theta_/(x+theta_))^alpha_
        return(res)
    } else {
        res <- (theta_/(x+theta_))^alpha_
        return(res)}
}

# A function that returns probability density function of Gamma\Inverse-Gamma composite random variables
dcompositeGInvG <- function(y, parm) {
    alphab <- parm[1]
    thetab <- parm[2]
    alphat <- parm[3]
    thetat <- parm[4]
    
    # provided the parameter values, calculate the threshold and weight parameters first
    theta  <- (alphab+alphat+sqrt((alphab+alphat)^2-4*thetat/thetab))/2*thetab
    phi    <- dgamma(theta, shape=alphab, scale=thetab)/
        pgamma(theta, shape=alphab, scale=thetab)/
        dinvgamma(theta, shape=alphat, scale=thetat)*
        pinvgamma(theta, shape=alphat, scale=thetat, lower.tail=FALSE)
    
    res    <- 1/(1+phi)*dgamma(y,        shape=alphab, scale=thetab)*   
        (y<=theta)/pgamma(theta,    shape=alphab, scale=thetab) +  
        phi/(1+phi)*dinvgamma(y,     shape=alphat, scale=thetat)*
        (y>theta)/pinvgamma(theta, shape=alphat, scale=thetat, lower.tail=FALSE)
    return(res)
}

# A function that returns the log-likelihood of a Gamma\Inverse-Gamma composite random variable
compositeGInvG_loglik <- function(y, parm) {
    #### using exp function to aviod constrained optimization
    parm[1] <- exp(parm[1])
    parm[2] <- exp(parm[2])
    parm[3] <- exp(parm[3])
    parm[4] <- exp(parm[4])
    
    res <- -sum(log(dcompositeGInvG(y, parm)))
    return(res) }

# Gamma\Log-Normal composite density and likelihood
dcompositeGL <- function(y, parm) {
    alphab  <- parm[1]
    thetab  <- parm[2]
    mut     <- parm[3]
    sigmat  <- parm[4]
    
    theta_solve <- function (x) {
        abs(alphab - x/thetab + (log(x)-mut)/sigmat^2)
    }
    theta <- optimize(f = theta_solve, interval = c(0,max(y)))$minimum
    
    phi    <- dgamma(theta, shape=alphab, scale=thetab)/
        pgamma(theta, shape=alphab, scale=thetab)/
        dnorm((log(theta)-mut)/sigmat)*(sigmat*theta)*
        pnorm((log(theta)-mut)/sigmat, lower.tail = FALSE)
    
    res    <- 1/(1+phi)*dgamma(y, shape=alphab, scale=thetab)/
        pgamma(theta, shape=alphab, scale=thetab)*(y<=theta) +  
        phi/(1+phi)*dnorm((log(y)-mut)/sigmat)/(sigmat*y)/
        pnorm((log(theta)-mut)/sigmat, lower.tail = FALSE)*(y>theta)
    return(res)
}

compositeGL_loglik <- function(y, parm) {
    parm[1] <- exp(parm[1])
    parm[2] <- exp(parm[2])
    parm[4] <- exp(parm[4])
    
    res <- -sum(log(dcompositeGL(y, parm)))
    return(res) }

# Gamma\Pareto composite density and likelihood
dcompositeGP <- function(y, parm) {
    alphab <- parm[1]
    thetab <- parm[2]
    alphat <- parm[3]
    thetat <- parm[4]
    
    theta  <- (thetab/2)*((alphab+alphat-thetat/thetab) + 
                              sqrt((alphab+alphat-thetat/thetab)^2+4*(alphab-1)*thetat/thetab))
    
    phi    <- dgamma(theta, shape=alphab, scale=thetab)/
        pgamma(theta, shape=alphab, scale=thetab)/
        dpareto(theta, alpha_ = alphat, theta_ = thetat)*
        ppareto(theta, alpha_ = alphat, theta_ = thetat, lower.tail = FALSE)
    
    res    <- 1/(1+phi)*dgamma(y, shape=alphab, scale=thetab)/
        pgamma(theta, shape=alphab, scale=thetab)*(y<=theta) +  
        phi/(1+phi)*dpareto(y, alpha_ = alphat, theta_ = thetat)/
        ppareto(theta, alpha_ = alphat, theta_ = thetat, lower.tail = FALSE)*(y>theta)
    return(res)
}

compositeGP_loglik <- function(y, parm) {
    parm[1] <- exp(parm[1])
    parm[2] <- exp(parm[2])
    parm[3] <- exp(parm[3])
    parm[4] <- exp(parm[4])
    
    res <- -sum(log(dcompositeGP(y, parm)))
    return(res) }

# Exponential\Log-Normal composite density and likelihood
dcompositeEL <- function(y, parm) {
    thetab  <- parm[1]
    mut     <- parm[2]
    sigmat  <- parm[3]
    
    theta_solve <- function (theta) {
        abs(-theta/thetab+1+(log(theta)-mut)/sigmat^2)
    }
    theta <- optimize(f = theta_solve, interval = c(0,100))$minimum
    
    phi    <- dexp(theta, rate = 1/thetab)/
        pexp(theta, rate = 1/thetab)/
        dnorm((log(theta)-mut)/sigmat)*(sigmat*theta)*
        pnorm((log(theta)-mut)/sigmat, lower.tail = FALSE)
    
    res    <- 1/(1+phi)*dexp(y, rate = 1/thetab)/
        pexp(theta, rate = 1/thetab)*(y<=theta) +  
        phi/(1+phi)*dnorm((log(y)-mut)/sigmat)/(sigmat*y)/
        pnorm((log(theta)-mut)/sigmat, lower.tail = FALSE)*(y>theta)
    return(res)
}

compositeEL_loglik <- function(y, parm) {
    parm[1] <- exp(parm[1])
    parm[3] <- exp(parm[3])
    
    res <- -sum(log(dcompositeEL(y, parm)))
    return(res) }

# Exponential\Inverse-Gamma composite density and likelihood
dcompositeEInvG <- function(y, parm) {
    thetab <- parm[1]
    alphat <- parm[2]
    thetat <- parm[3]
    
    theta  <- ((alphat+1)+sqrt((alphat+1)^2-4*thetat/thetab))*thetab/2
    phi    <- dexp(theta, rate = 1/thetab)/
        pexp(theta, rate = 1/thetab)/
        dinvgamma(theta, shape=alphat, scale=thetat)*
        pinvgamma(theta, shape=alphat, scale=thetat, lower.tail=FALSE)
    
    res    <- 1/(1+phi)*dexp(y, rate = 1/thetab)/
        pexp(theta, rate = 1/thetab)*(y<=theta) +
        phi/(1+phi)*dinvgamma(y, shape=alphat, scale=thetat)/
        pinvgamma(theta, shape=alphat, scale=thetat, lower.tail = FALSE)*(y>theta)
    return(res)
}

compositeEInvG_loglik <- function(y, parm) {
    parm[1] <- exp(parm[1])
    parm[2] <- exp(parm[2])
    parm[3] <- exp(parm[3])
    
    res <- -sum(log(dcompositeEInvG(y, parm)))
    return(res) }

# Exponential\Pareto composite density and likelihood
dcompositeEP <- function(y, parm) {
    thetab <- parm[1]
    alphat <- parm[2]
    thetat <- parm[3]
    
    theta  <- (alphat+1)*thetab - thetat
    phi    <- dexp(theta, rate = 1/thetab)/
        pexp(theta, rate = 1/thetab)/
        dpareto(theta, alpha_=alphat, theta_=thetat)*
        ppareto(theta, alpha_=alphat, theta_=thetat, lower.tail=FALSE)
    
    res    <- 1/(1+phi)*dexp(y, rate = 1/thetab)/
        pexp(theta, rate = 1/thetab)*(y<=theta) +
        phi/(1+phi)*dpareto(y, alpha_=alphat, theta_=thetat)/
        ppareto(theta, alpha_=alphat, theta_=thetat, lower.tail = FALSE)*(y>theta)
    return(res)
}

compositeEP_loglik <- function(y, parm) {
    parm[1] <- exp(parm[1])
    parm[2] <- exp(parm[2])
    parm[3] <- exp(parm[3])
    
    res <- -sum(log(dcompositeEP(y, parm)))
    return(res) }

# The following are the inverse of the CDF of the composite random variables
qcompositeGL <- function(p, parm) {
    alphab      <- parm[1]
    thetab      <- parm[2]
    mut         <- parm[3]
    sigmat      <- parm[4]
    
    # Threshold & weight parameter values
    theta_solve <- function (x) {
        abs(alphab - x/thetab + (log(x)-mut)/sigmat^2)}
    theta       <- optimize(f = theta_solve, interval = c(0,1000))$minimum
    
    phi         <- dgamma(theta, shape=alphab, scale=thetab)/
        pgamma(theta, shape=alphab, scale=thetab)/
        dnorm((log(theta)-mut)/sigmat)*(sigmat*theta)*
        pnorm((log(theta)-mut)/sigmat, lower.tail = FALSE)
    
    # a indicator indicates the input probabilities fall in the head or tail of the distribution
    indicator.p <- 1/(1+phi)
    
    # inverse of the CDF
    f_inv <- function (p) {
        res_inv <- rep(0, length(p))
        res_inv[indicator.p > p]    <- 
            qgamma((1+phi)*p[indicator.p > p]*pgamma(theta, shape=alphab, scale=thetab),
                   shape = alphab, scale = thetab)
        res_inv[!indicator.p > p] <- 
            qlnorm((1+phi)/phi*(1-p[!indicator.p > p])*
                       plnorm(theta, mut, sigmat, lower.tail = F), 
                   meanlog = mut, sdlog = sigmat, lower.tail = F)
        res_inv
    }
    
    res <- rep(0,length(p))
    res[p!=0] <- f_inv(p[p!=0])
    names(res) <- NULL
    
    return(res)
}

qcompositeGP <- function(p, parm) {
    alphab      <- parm[1]
    thetab      <- parm[2]
    alphat      <- parm[3]
    thetat      <- parm[4]
    
    theta  <- (thetab/2)*((alphab+alphat-thetat/thetab) + 
                              sqrt((alphab+alphat-thetat/thetab)^2+4*(alphab-1)*thetat/thetab))
    
    phi    <- dgamma(theta, shape=alphab, scale=thetab)/
        pgamma(theta, shape=alphab, scale=thetab)/
        dpareto(theta, alpha_ = alphat, theta_ = thetat)*
        ppareto(theta, alpha_ = alphat, theta_ = thetat, lower.tail = FALSE)
    
    indicator.p <- 1/(1+phi)
    
    f_inv <- function (p) {
        res_inv <- rep(0, length(p))
        res_inv[indicator.p > p]    <- qgamma((1+phi)*p[indicator.p > p]*
                                                  pgamma(theta, shape=alphab,
                                                         scale=thetab),
                                              shape = alphab, scale = thetab)
        res_inv[!indicator.p > p] <- qpareto((1+phi)/phi*(1-p[!indicator.p > p])*
                                                 ppareto(theta, alphat, thetat, F),
                                             alphat, thetat, F)
        res_inv
    }
    
    res <- rep(0,length(p))
    res[p!=0] <- f_inv(p[p!=0])
    names(res) <- NULL
    return(res)
}



## tau-leap Gillespie algorithm function
tauLeapG <- function(beta, # transmission rate
                     theta, # dispersal scale
                     removal.rate, #removal rate
                     omega, #attractiveness
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputing the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=1, # stoping condition 1: incidence lvl
                     t.end=3650, # stoping condition 2: time after first simulated time step
                     Nvectors=Nvectors, # surface area occupied by one host
                     delta.t=14, # time step
                     ppp, # point pattern as a ppp object, optinally with marks 1/0 for infected/healthy
                     dist.mat=NULL){ # matrix distance if its computation is to be avoided here (for e.g. repeated calls)
  
  ## Initialise host states
  state <- rep("S", ppp$n)
  
  inf.start <- max(1, round(ppp$n * q0))
  
  initial <- sample(ppp$n, inf.start)
  
  state[initial] <- "I"
  ## compute distance matrix if not provided
  if (is.null(dist.mat)){ 
    ## add the kernel computation that can be added directly on the dist matrix to reduce comp time
    dist.mat <- exp(-pairdist(ppp)^b / theta^b)
    diag(dist.mat) <- NA
  }
  attractiveness <- rep(1, length(ppp$x))
  attractiveness[state == "S"] <- 1
  attractiveness[state == "I"] <- omega
  attractiveness[state == "R"] <- 0
  ## function that compute infection event probability, based on the dispersal kernel
  k.norm <- beta * (b/(2*pi*theta^2*gamma(2/b))) # constant part of the exponential power kernel
   p <- attractiveness / sum(attractiveness)
  vectors <- as.vector(rmultinom(
      n = 1,
      size = Nvectors,
      prob = p
    ))
  infection <- function(infected, susceptible, vectors, dist){
    inf <- k.norm * dist[infected, susceptible, drop = FALSE]
    
    for(i in seq_len(nrow(inf))){
      inf[i, ] <- inf[i, ] * vectors[infected][i]
    }
    
    inf[is.na(inf)] <- 0
    
    inf
  }
  ##beginning of the EPIDEMIC
  ## starting time
  time <- 0
  
  ## record infection times
  infection_time <- rep(NA, ppp$n)
  infection_time[state == "I"] <- 0
  
  removal_time <- rep(NA, ppp$n)
  
  #time after which infection has occurred
  infected_hosts <- which(state=="I")
  infection_age <- time - infection_time[infected_hosts]
  
  ## inititate the heavy dataframe that will aggregate all changes
  df.big <- data.frame(
    time = 0,
    who = which(state == "I")
  )
  history <- data.frame(
    time = 0,
    S = sum(state == "S"),
    I = sum(state == "I"),
    R = sum(state == "R")
  )
  ## computation loop
  while (
    any(state == "S") &
    time <= t.end &
    mean(state == "I") < q.end
  ){
    ## infection event probaility
    events <- infection(
      infected    = state=="I",
      susceptible = state=="S",
      vectors     = vectors,
      dist        = dist.mat
    )

    
    ## random proisson draws
    new.infected <- which(state == "S")[
      rpois(
        n = sum(state == "S"),
        lambda = apply(events, 2, sum) * delta.t
      ) > 0
    ]
    ## change marks of newly infected
    state[new.infected] <- "I"
    ##add infection time
    infection_time[new.infected] <- time
    #removal events according to probability argument
    infected_hosts <- which(state == "I")
    remove <- infected_hosts[
      runif(length(infected_hosts)) <
        (1 - exp(-removal.rate * delta.t)) #GAMMA UNDEFINED OR NEGATIVE
    ]
    
    state[remove] <- "R"
    removal_time[remove] <- time
    ## increment time
    time <- time + delta.t
    ## if some infection, increment the big dataframe
    if (length(new.infected) > 0){
      df.big <- rbind(
        df.big,
        data.frame(
          time = time,
          who = new.infected
        )
      )
    }
    
    S <- sum(state == "S")
    I <- sum(state == "I")
    R <- sum(state == "R")
    
    history <- rbind(
      history,
      data.frame(
        time = time,
        S = sum(state == "S"),
        I = sum(state == "I"),
        R = sum(state == "R")
      )
    )
    ## print a dot per new infection
    # cat(paste0(rep('.', length(new.infected)), collapse = '')) ## comment for quiet
  }
  
  ## make compact, time only, version of the big dataframe
  times.i <- unique(df.big[,1])
  times.d <- times.i + sigma
  times <- sort(unique(c(times.i, times.d)))
  infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
  detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
  df.small <- data.frame(time=times, infected=infected, detectable=detectable)
  ## out put the simplified time series, and the big one
  list(df.small[df.small$time <= max(df.big$time),], 
       df.big,
       infection_time,
       removal_time,
       history) 
} 

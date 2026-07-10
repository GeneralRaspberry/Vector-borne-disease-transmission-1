

#runsim <- function(){
history.all <- data.frame()
r.all <- data.frame()   # later for growth rates
for (Nvectors in vector.values){
  
  for (omega in omega.values){

    
    cat("\nRunning Nvectors =", Nvectors,
        " omega =", omega, "\n")
    
    data <- data.frame(
      x = landscape$x,
      y = landscape$y,
      id = 1:hosts
    )
    history.list <- vector("list", iterations)
    
    for(i in 1:iterations){
      
      marks(landscape) <- sample(c(TRUE, rep(FALSE, hosts-1)))
      output <- tauLeapG(
        beta = beta,
        theta = theta,
        b = b,
        omega = omega,
        sigma = sigma,
        delta.t = 14,
        ppp = landscape,
        Nvectors = Nvectors,
        removal.rate = removal.rate
      )
      
      infection_times <- output[[3]]
      removal_times   <- output[[4]]
      history <- output[[5]]
      
      history$sim <- i
      history$Nvectors <- Nvectors
      history$omega <- omega
      
      history.list[[i]] <- history
      
      
      
      data <- cbind(
        data,
        infection_times,
        removal_times
      )
      
      cat(".")
    }
    history.all <- bind_rows(
      history.all,
      bind_rows(history.list)
    )
    cat("\nFinished combination.\n")
  }
}
#}

history.long <-
  history.all %>%
  pivot_longer(
    cols=c(S,I,R),
    names_to="State",
    values_to="Hosts"
  )

history.long$State <- factor(
  history.long$State,
  levels=c("S","I","R")
)

library(ggplot2)

ggplot(history.long,
       aes(time,
           Hosts,
           group=interaction(sim,State),
           colour=State))+
  
  geom_line(alpha=.05,
            linewidth=.35)+
  
  facet_grid(
    Nvectors~omega,
    labeller=label_both
  )+
  
  scale_colour_manual(
    values=c(
      S="forestgreen",
      I="firebrick",
      R="grey40"
    ),
    labels=c(
      S="Susceptible",
      I="Infected",
      R="Removed"
    )
  )+
  
  theme_classic()

##a ten line classic

###add The SIR compartments up over time


## name the additional columns
##colnames(data)[-(1:3)] <- paste0(1:iterations)

## to fit a logistic we need to convert this to number of infected a each time
##data_long <- data%>%pivot_longer(cols = 4:(iterations+3), values_to = 'time', names_to = 'sim')
##times <- sort(unique(data_long$time))

## make a logistic df from this data
data_logistic <- history.all %>% 
    rename(infected = I)
  
## all curves have the same number of points
##ggplot(data_logistic) + geom_line(aes(x=time, y=infected, colour=sim))+
##theme(legend.position = "none")
  ## prepare a logistic function of r to fit
  logis <- function(t, r, K=1, s=0, q0){
    pmin(
      K*q0*exp(r*(t+s)) / (K + q0*(exp(r*(t+s)) - 1)),
      K) # numerical errors can happen for high r and sigma
  }
  
  eval <- function(r, df){
    sum((logis(r=r, t=df$time, K=1000, q0=1) - df$infected)^2) ## sum of square errors between predictions and observations
  }
  
  for(Nv in unique(history.all$Nvectors)){
    
    for(om in unique(history.all$omega)){
      
      temp <- subset(history.all,
                     Nvectors == Nv &
                       omega == om)
      
      for(i in unique(temp$sim)){
        
        epidemic <- subset(temp, sim == i)
        
        epidemic <- subset(epidemic,
                           I > 0 & I < 250)
        
        if(nrow(epidemic) > 5){
          
          r.calc <- optimize(
            eval,
            interval = c(0,1),
            df = data.frame(
              time = epidemic$time,
              infected = epidemic$I
            )
          )$minimum
          
          r.all <- rbind(
            r.all,
            data.frame(
              Nvectors = Nv,
              omega = om,
              sim = i,
              r = r.calc
            )
          )
          
        }
        
      }
      
    }
    
  }
  r.summary <- r.all %>%
    group_by(Nvectors, omega) %>%
    summarise(
      mean.r = mean(r),
      sd.r = sd(r),
      .groups = "drop"
    )

 #plot 
  ggplot(r.all, aes(x = r)) +
    
    geom_density(
      fill = "#4F81BD",
      colour = "#1F4E79",
      alpha = 0.8,
      linewidth = 0.8,
      adjust = 1.1,
      na.rm = TRUE
    ) +
    
    facet_grid(
      Nvectors ~ omega,
      labeller = labeller(
        omega = label_both,
        Nvectors = label_both
      )
    ) +
    
    labs(
      x = expression("Estimated growth rate, "*italic(r)),
      y = "Density"
    ) +
    
    theme_bw(base_size = 14) +
    
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )+
    xlim(c(0,.02))
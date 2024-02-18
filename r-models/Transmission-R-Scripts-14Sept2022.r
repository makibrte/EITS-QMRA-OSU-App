
# ---- LOAD PACKAGES ----

  if("deSolve" %in% installed.packages()==FALSE){install.packages("deSolve")}


# ---- CONSTANTS ----
  # -- General constants --
    N = 1000         # Fixed population
    scenarios <- c('Air', 'Frequent-Fomite', 'Infrequent-Fomite')

  # -- Constants for Air --
    g.air <- 0.2         # gamma, recovery rate 1/day
    p.air <- 0.0517      # pi, infectivity
    m.air <- 8.64        # mu, elimination rate 1/day
    a.air <- 693         # alpha, deposit rate, pathogens/infected/day
    r.air <- 0.0000877   # ro pick up 1/person/day

  # -- Constants for Frequently Touched Fomites --
    g.freq.fomite <- 0.2          # gamma, recovery rate 1/day
    p.freq.fomite <- 0.0000693    # pi, infectivity
    m.freq.fomite <- 2.88         # mu, elimination rate 1/day
    a.freq.fomite <- 5244         # alpha, deposit rate, pathogens/infected/day
    r.freq.fomite <- 0.297        # ro pick up 1/person/day

  # -- Constants for Infrequently Touched Fomites --
    g.infreq.fomite <- 0.2         # gamma, recovery rate 1/day
    p.infreq.fomite <- 0.0000693;  # pi, infectivity
    m.infreq.fomite <- 2.88        # mu, elimination rate 1/day
    a.infreq.fomite <- 1040177     # alpha, deposit rate, pathogens/infected/day
    r.infreq.fomite <- 0.0000145   # ro pick up 1/person/day


# ---- FUNCTIONS ----
  # -- R0 --
    R_0 <- function(Pars,N){(a/g)*(r*N/(m+r*N))*p}

  # -- EITS --

    EITSModel <- function(time, y, params){
    S <- y[1]
    I <- y[2]
    R <- y[3]
    E <- y[4]

    g <- params['g']
    p <- params['p']
    m <- params['m']
    a <- params['a']
    r <- params['r']

    dSdt <- (-r)*p*S*E
    dIdt <- r*p*S*E-(g*I)
    dRdt <- g*I
    dEdt <- a*I-E*(1000*r+m)

    dydt <- c(dSdt,dIdt,dRdt,dEdt)

    list(dydt)
    }

# ---- ODE Solution ----

  # -- Set Initial coditions, define parameters --

    ini.cond <- c(S=(N-1), I=1, R=0,E=0)

    times <- seq(0,150,0.5)
    para.air <- c(g = g.air, p = p.air, m = m.air, a = a.air, r = r.air)
    para.freq.fomite <- c(g = g.freq.fomite, p = p.freq.fomite,
      m = m.freq.fomite, a = a.freq.fomite, r = r.freq.fomite)
    para.infreq.fomite <- c(g = g.infreq.fomite, p = p.infreq.fomite,
      m = m.infreq.fomite, a = a.infreq.fomite, r = r.infreq.fomite)

  # -- deSolve package called for ode solution --

    EITS.Air <- deSolve::ode(ini.cond,times,EITSModel,para.air)
    EITS.freq.fomite <- deSolve::ode(ini.cond,times,EITSModel,para.freq.fomite)
    EITS.infreq.fomite <- deSolve::ode(ini.cond,times,EITSModel,para.infreq.fomite)


# ---- OUTPUTS

  # -- Summary Stats --

    Air.summary <- summary(EITS.Air)
      colnames(Air.summary) <- c("Susceptible", "Infected", "Recovered", "Exposed")
    Freq.fomite.summary <- summary(EITS.freq.fomite)
      colnames(Freq.fomite.summary) <- c("Susceptible", "Infected", "Recovered", "Exposed")
    Infreq.fomite.summary <- summary <- summary(EITS.infreq.fomite)
      colnames(Infreq.fomite.summary) <- c("Susceptible", "Infected", "Recovered", "Exposed")

    # -- Write stats to disk --
      write.csv(Air.summary, file = "Air-Summary-Stats.csv")
      write.csv(Freq.fomite.summary, file = "Ferquently-Touched-Fomite-Summary-Stats.csv")
      write.csv(Infreq.fomite.summary, file = "Inferquently-Touched-Fomite-Summary-Stats.csv")

  # -- Plots --

      # Method for plotting better than lots of cutting and pasting ?sprintf if you want to learn how
    for(i in 1:length(scenarios)){
    I <- scenarios[i]
    if(i==1){Z <- EITS.Air}; if(i==2){Z <- EITS.freq.fomite}; if(i==3){Z <- EITS.infreq.fomite}
    png(sprintf('%s-Incidence.png',I),500,500)
           plot(Z[,1], Z[,3], xlab = expression(bold(Time~(days))),
        ylab = expression(bold(Prevalence)), main = expression(bold(Prevalence~Time~Series)),
        type = 'l', col = 'blue', lwd=3)
    dev.off()

    png(sprintf("%s-Cumulative-Incidence.png",I),500,500)
      plot(Z[,1], (Z[,3]+Z[,4]), xlab = expression(bold(Time~(days))),
        ylab = expression(bold(Cumulatiive~Incidence)), main = expression(bold(Cumulative~Incidence~Time~Series)),type = 'l', col = 'blue', lwd=3  )
    dev.off()

          leg.text <- c('Susceptible','Cumulative Incidence','Incidence')
    png(sprintf("%s-Combined-Plot.png",I),500,500)
      plot(Z[,1], Z[,2], ylim=c(min(Z[,3]), max(Z[,2])), xlab = expression(bold(Time~(days))),
        ylab = expression(bold(Number~of~People)), main = expression(bold(Combined~Incidence~Time~Series)),type = 'l', col = 'blue', lwd=3  )
          lines(Z[,1], (Z[,3]+Z[,4]), col = 'red', lwd = 3)
          lines(Z[,1], Z[,3], col = 'green', lwd =3)
          legend('topright',leg.text, col = c('blue','red','green'), lty = c(1,1,1))
    dev.off()
    }


# END OF CODE
	# You have reached the end of your journey. Take a rest, consider the concept of Memento Mori.
	# Memento Mori is a concept from Stoic philosophy, where it forces you to consider the vital and changable
	# It is not a morbid one as many shortsighted people may think. Rather it is a way of contemplating
	# and then ensure that you engage in those things that truly matter. Family first and always is part of
	# my Memento Mori, I work to live, not live to work. Go get a run, play with your children, workout,
	# or do something else that is good for you, not coding and reading this.

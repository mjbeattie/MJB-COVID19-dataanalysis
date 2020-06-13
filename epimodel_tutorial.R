# EpiModel Tutorial
# May 23, 2020
#
# Taken from Jeness, "Disease over Networks", Journal of Statistical Software, 2018

library("EpiModel")

set.seed(12345)

# SIS Model Example: closed population, some assortive mixing
# 1000 persons, two equally sized risk groups.  Independent model.

# Estimate network structure
nw <- network::network.initialize(n=1000,directed=FALSE)
nw <- network::set.vertex.attribute(nw, "risk", rep(0:1, each=500))

formation <- ~ edges + nodefactor("risk") + nodematch("risk") + concurrent

# Data is passed to the formation formula via target statistics
# These stats are expected values for each time step of the network
target.stats <- c(250, 375, 225, 100)

# Set dissolution parameter.  Done here by specifying the mean
# duration of the edges in number of time steps.  Offset means coef is explicitly set.
coef.diss <- dissolution_coefs(dissolution= ~ offset(edges), duration=80)

# Estimate coefficients for formulation and dissolution from tergm call
est1 <- netest(nw, formation, target.stats, coef.diss)

# Diagnose the network to ensure it fits with observed statistics
dx <- netdx(est1, nsims = 10, nsteps = 1000)
dx

# Plot the diagnostics
plot(dx)
par(mfrow = c(1, 2)) 
plot(dx, type = "duration") 
plot(dx, type = "dissolution")

# RUN THE EPIDEMIC SIMULATION ON THE CONSTRUCTED NETWORK
# Set initial number of infected nodes
init <- init.net(i.num = 50)

# Base SIS model takes three parameters, infection probability (inf.prob),
# act rate (act.rate), recovery rate (rec.rate).  These are setup using
# the param.net helper function.  Mean recovery rate is the inverse of the
# time spent infected.
param <- param.net(inf.prob=0.1, act.rate=5, rec.rate=0.02)


# Control settings include epidemic type, etc.
control <- control.net(type="SIS", nsteps=500, nsims=10, epi.by="risk")

# Pass fitted model and control parameters to netsim
sim1 <- netsim(est1, param, init, control)

# Print model
sim1

# Print simulation for a specific time step
summary(sim1, at=500)

# Default plot
plot(sim1)

# Plot incidence of infection and recovery
plot(sim1, y=c("si.flow","is.flow"), legend=TRUE)

# Plot results by risk group
plot(sim1, y=c("i.num.risk0", "i.num.risk1"), legend=TRUE)

# Plot network at specific time steps
par(mfrow=c(1,2), mar=c(0,0,1,0))
plot(sim1, type="network", at=1, sims="mean", col.status=TRUE,
     main="Prevalence at t1")
plot(sim1, type="network", at=500, sims="mean", col.status=TRUE,
     main="Prevalence at t500")


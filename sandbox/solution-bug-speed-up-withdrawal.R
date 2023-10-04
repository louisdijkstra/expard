nrow = 5

m1 <- matrix(rep(0, 10), nrow = nrow)
m1[2,2] <- 1

i = which(m1 == 1)

m2 <- matrix(rep(0, 10), nrow = nrow)
m2[2,2] <- "hello"

m2[i]

tm1 <- t(m1)

# determine the time_steps_ago for every patient
time_steps_ago <- unlist(
  lapply(1:nrow, function(i) { 
    tm1[,i]
  }))

j = which(time_steps_ago == 1)

m2[j]


l <- lapply(1:nrow, function(i) { 
  m1[i,]
})

as.matrix(l)
## THIS IS THE TRICK!
which(do.call(rbind, l) == 1)


if (FALSE){
  library(ggplot2)
  n <- 100 
  set.seed(10)
  y <- abs(rnorm(n))
  P <- abs(rnorm(n))
  S <- abs(cbind(rnorm(n), rnorm(n), rnorm(n), rnorm(n), rnorm(n)))
  dimnames(S)[[2]] <- c("a", "b", "c", "d", "e")
  
  data = as.data.frame(cbind(y, P, S))
  loss <- "y"
  score <- c("a", "b", "c")
  base = c("P", "c")
  
  
  x1 <- gini(loss = "y", score = c("a", "b", "c", "d", "e"), base = "P", data = data)
  x1
  plot(x1)
  plot(x1, overlay = FALSE)
  
  
  x3 <- gini(loss = "y", score = c("b", "c", "d", "e"), base = NULL, data = data)                 
  x3
  plot(x3)
  plot(x3, overlay = FALSE)
  
  plot(x)                 
  plot(x, overlay = FALSE)
  
  x = ans
  digits = 3
  plot(x)
}
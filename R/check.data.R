check.data <-
function(dat)
{
  # test if data set contains continuous variables
  cont = apply(dat, 2, function(x18) sum(abs(x18-round(x18))))
  if (any(cont!=0)){
    stop("data set contains continuous variables: ", paste0("V",which(cont!=0),collapse=','))
  }
  
  # test if data set contains constant variables
  Vari = apply(dat, 2, var)
  if (any(Vari==0)){
    stop("data set contains variables with zero variance: ", paste0("V",which(Vari==0),collapse=','))
    }

  cat("All clear!",'\n')
}


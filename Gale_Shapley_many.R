#' This is a modified version of Gale-Shapley stable matching algorithm that can do 
#' one-to-one, one-to-many, and many-to-many matching. Outside option (self-matched)
#' is allowed.
#' 
#' @param U The utility matrix for the women's side. Each row is a woman, each column is a man.
#' The matrix entry (i,j) is the utility that woman \code{i} gains from pairing with man \code{j}. 
#' In other words, the utility is computed from woman \code{i}'s perspective.
#' @param V The utility matrix for the men's side. Each column is a man, each row is a woman.
#' The matrix entry (i,j) is the utility that man \code{j} gains from pairing with woman \code{i}. 
#' In other words, the utility is computed from man \code{j}'s perspective.
#' @param U_slots The number of choice situations each woman has.
#' @param V_slots The number of choice situations each man has.
#' @return This function returns the following matrix: 
#' \item{mu}{The matching matrix, where 1 represents a pairing, 0 otherwise. 
#' Each row is a woman, each column is a man. The order of the rows is the same as the 
#' rows in \code{U}. The order of the columns is the same as the columns in \code{V}.}
#' @seealso fitrpm_R_CP
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples
#' 

# return indices with the n largest elements greater than 0
maxNIdx = function(x, n) {
  res = sort(x, method = "radix", decreasing = TRUE, index.return = TRUE)
  # idx = na.omit(res$ix[res$x > 0][1:n])
  idx = res$ix[res$x > 0][1:n]
  idx
}

# return the nth largest value greater than 0 (or return 0)
maxNVal = function(x, n) {
  res = sort(x, method = "radix", decreasing = TRUE, index.return = TRUE)
  val = na.omit(res$x[res$x > 0][1:n])
  if (length(val)>0) {
    val = min(val)
  } else {
    val = 0
  }
  val
}

Gale_Shapley_many <- function(U,V, U_slots, V_slots){

# first argument proposing side
# proposing side: individuals correspond to rows

nw <- nrow(U)
nm <- ncol(U)
U_temp <- U
nmax <- 10*nw*nm


for (i in 1:nmax){

 # check dimensions!
 # no need to worry about ties
    # Prop <- sweep(U_temp, 2, apply(rbind(U_temp,0),2,max),"==") # incorrect swap of first and second dimensions
    # Prop <- sweep(U_temp, 1, apply(cbind(U_temp,0),1,max),"==")
    Prop = matrix(FALSE, nrow=nw, ncol=nm)
    colIdx = apply(U_temp, 1, maxNIdx, U_slots)
    if (U_slots==1) { 
      Prop[cbind(1:nw, colIdx)] = TRUE
    } else {
      for (jj in 1:U_slots) {
        Prop[cbind(1:nw, colIdx[jj,])] = TRUE
      }
    }
    
     
    # Rej <- (Prop*V < sweep(Prop,1,apply(cbind(Prop*V,0),1,max),"*")) # incorrect swap of first and second dimensions
    # Rej <- (Prop*V < sweep(Prop,2,apply(rbind(Prop*V,0),2,max),"*"))
    Rej <- (Prop*V < sweep(Prop,2,apply(Prop*V,2,maxNVal, V_slots),"*"))
    
    
    U_temp[Rej&Prop] = -1

    if (all(Rej == FALSE)) {
      break
    }
    # if (sum(sum(Rej))==0){
    #     break
    # }
}

# mu <- sweep(U_temp, 2, apply(rbind(U_temp,0),2,max),"==") # incorrect swap of first and second dimensions
mu <- sweep(U_temp, 1, apply(cbind(U_temp,0),1,max),"==")
if (U_slots==1 & V_slots == 1) {
  mu <- sweep(U_temp, 1, apply(cbind(U_temp,0),1,max),"==")
} else {
  mu = Prop
}

#mu <- cbind(1:nrow(mu),apply(mu,1,function(x){match(1,x,nomatch=0)}))
#mu[mu[,2]>0,]
1*mu
}

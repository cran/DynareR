#' Write a new `mod` file.
#'
#' Use `write_mod(model,code)` if you want the `Dynare` file to live in the current working directory.
#' Use `write_mod(model,code,path)` if you want the `Dynare` file to live in the path different from the current working directory.
#'
#' @inheritParams run_dynare
#' @return Set of \code{Dynare} (open-source software for DSGE modelling) outputs
#' @examples library(DynareR)
#' \dontrun{
#' DynareCodes='var y, c, k, a, h, b;
#' varexo e, u;
#' parameters beta, rho, alpha, delta, theta, psi, tau;
#' alpha = 0.36;
#' rho   = 0.95;
#' tau   = 0.025;
#' beta  = 0.99;
#' delta = 0.025;
#' psi   = 0;
#' theta = 2.95;
#' phi   = 0.1;
#' model;
#' c*theta*h^(1+psi)=(1-alpha)*y;
#' k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
#'           *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
#' y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
#' k = exp(b)*(y-c)+(1-delta)*k(-1);
#' a = rho*a(-1)+tau*b(-1) + e;
#' b = tau*a(-1)+rho*b(-1) + u;
#' end;
#' initval;
#' y = 1.08068253095672;
#' c = 0.80359242014163;
#' h = 0.29175631001732;
#' k = 11.08360443260358;
#' a = 0;
#' b = 0;
#' e = 0;
#' u = 0;
#' end;
#'
#' shocks;
#' var e; stderr 0.009;
#' var u; stderr 0.009;
#' var e, u = phi*0.009*0.009;
#' end;
#'
#' stoch_simul;'
#' model<-"example1" # This is "example1" of the `Dynare` example files
#' code<-DynareCodes
#' write_mod(model,code)
#'}
#' @seealso write_dynare eng_dynare run_model run_dynare
#' @keywords documentation
#' @export
write_mod <- function(model,code,path="") {
  #model=paste0("DynareR_",model)
if(path==""){
  file.create(paste0(model, '.', "mod"))
  f <-paste0(model, '.', "mod")
  }else{
    file.create(paste0(path,"/",model, '.', "mod"))
  f <-paste0(path,"/",model, '.', "mod")
  }
  writeLines(code, f)
}




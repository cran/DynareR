## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval=F,
  echo = T
  
)


library(magrittr)
library(DynareR)

## ----installation,eval=FALSE--------------------------------------------------
#  install.packages("DynareR")
#  
#            OR
#  
#  devtools::install_github('sagirumati/DynareR')

## ----IRF,out.extra='trim={0cm 7cm 0cm 7cm},clip',fig.cap="Another of figure generated from Dynare software"----
#  
#  include_IRF(model="bkk",IRF = "E_H2")
#  
#  # Alternatively, use the path argument
#  
#  include_IRF(path="bkk/bkk/graphs/bkk_IRF_E_H2.pdf")
#  

## ----exampleDyn---------------------------------------------------------------
#  dynareCodes='var y, c, k, a, h, b;
#  varexo e, u;
#  parameters beta, rho, alpha, delta, theta, psi, tau;
#  alpha = 0.36;
#  rho   = 0.95;
#  tau   = 0.025;
#  beta  = 0.99;
#  delta = 0.025;
#  psi   = 0;
#  theta = 2.95;
#  phi   = 0.1;
#  model;
#  c*theta*h^(1+psi)=(1-alpha)*y;
#  k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
#            *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
#  y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
#  k = exp(b)*(y-c)+(1-delta)*k(-1);
#  a = rho*a(-1)+tau*b(-1) + e;
#  b = tau*a(-1)+rho*b(-1) + u;
#  end;
#  initval;
#  y = 1.08068253095672;
#  c = 0.80359242014163;
#  h = 0.29175631001732;
#  k = 11.08360443260358;
#  a = 0;
#  b = 0;
#  e = 0;
#  u = 0;
#  end;
#  
#  shocks;
#  var e; stderr 0.009;
#  var u; stderr 0.009;
#  var e, u = phi*0.009*0.009;
#  end;
#  
#  stoch_simul;'
#  
#  
#  write_dyn(code=dynareCodes, model="example1")
#  
#  write_dyn(code=dynareCodes,model="DynareR/write_dyn/example1")

## ----examplemod---------------------------------------------------------------
#  DynareCodes='var y, c, k, a, h, b;
#  varexo e, u;
#  parameters beta, rho, alpha, delta, theta, psi, tau;
#  alpha = 0.36;
#  rho   = 0.95;
#  tau   = 0.025;
#  beta  = 0.99;
#  delta = 0.025;
#  psi   = 0;
#  theta = 2.95;
#  phi   = 0.1;
#  model;
#  c*theta*h^(1+psi)=(1-alpha)*y;
#  k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
#            *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
#  y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
#  k = exp(b)*(y-c)+(1-delta)*k(-1);
#  a = rho*a(-1)+tau*b(-1) + e;
#  b = tau*a(-1)+rho*b(-1) + u;
#  end;
#  initval;
#  y = 1.08068253095672;
#  c = 0.80359242014163;
#  h = 0.29175631001732;
#  k = 11.08360443260358;
#  a = 0;
#  b = 0;
#  e = 0;
#  u = 0;
#  end;
#  
#  shocks;
#  var e; stderr 0.009;
#  var u; stderr 0.009;
#  var e, u = phi*0.009*0.009;
#  end;
#  
#  stoch_simul;'
#  
#  
#  write_mod(model="example1",code=dynareCodes)
#  
#  write_mod(code=dynareCodes,model="DynareR/write_mod/example1")
#  

## ----exampleRunDynare---------------------------------------------------------
#  DynareCodes='var y, c, k, a, h, b;
#  varexo e, u;
#  parameters beta, rho, alpha, delta, theta, psi, tau;
#  alpha = 0.36;
#  rho   = 0.95;
#  tau   = 0.025;
#  beta  = 0.99;
#  delta = 0.025;
#  psi   = 0;
#  theta = 2.95;
#  phi   = 0.1;
#  model;
#  c*theta*h^(1+psi)=(1-alpha)*y;
#  k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
#            *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
#  y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
#  k = exp(b)*(y-c)+(1-delta)*k(-1);
#  a = rho*a(-1)+tau*b(-1) + e;
#  b = tau*a(-1)+rho*b(-1) + u;
#  end;
#  initval;
#  y = 1.08068253095672;
#  c = 0.80359242014163;
#  h = 0.29175631001732;
#  k = 11.08360443260358;
#  a = 0;
#  b = 0;
#  e = 0;
#  u = 0;
#  end;
#  
#  shocks;
#  var e; stderr 0.009;
#  var u; stderr 0.009;
#  var e, u = phi*0.009*0.009;
#  end;
#  
#  stoch_simul;'
#  
#  run_dynare(code=DynareCodes,model="example1",import_log = T)
#  run_dynare(code=DynareCodes,model="DynareR/run_dynare/example1")
#  

## ----exampleRunModels---------------------------------------------------------
#  
#  demo(agtrend)
#  demo(bkk)
#  demo(example1)
#  
#  # Provide the list of the `Dynare` files in a vector
#  # Ensure that "agtrend.mod", "bkk.mod" and "example1.mod"
#  # live in the current working directory
#  
#  # Copy the dynare files to the current working directory
#  
#  lapply(c("agtrend","bkk","example1"),\(x) file.copy(paste0(x,"/",x,".mod"),"."))
#  
#  run_models(c("agtrend","bkk","example1")) # Run the models in the vector.
#  
#  

## ----runAllDynare-------------------------------------------------------------
#  run_models() # Run all models in Current Working Directory.

## ----runAllDynare1------------------------------------------------------------
#  
#  # Copy the dynare files to the 'DynareR/run_dynare' directory
#  
#  lapply(c("agtrend","bkk","example1"),\(x) file.copy(paste0(x,".mod"),"DynareR/run_dynare"))
#  
#  run_models(model = 'DynareR/run_dynare*') # notice the * at the end
#  

## ----importLog----------------------------------------------------------------
#  import_log(model="bkk")
#  
#  import_log(path="bkk/bkk.log")
#  
#  knitr::kable(dynare$bkk$autocorrelation)  %>% kableExtra::kable_styling(latex_options = c("basic","hold_position","scale_down")) %>%
#   kableExtra::footnote(general="Some footnote with equation $\\alpha x^2+\\beta x+c=0$", general_title = "*",footnote_as_chunk=T,threeparttable=T,escape=F) %>%
#  kableExtra::row_spec(0,bold=T)

## ----dynareVersion,eval=F-----------------------------------------------------
#  set_dynare_version("6-unstable-2022-04-03-0800-700a0e3a")
#  

## ----OctavePath---------------------------------------------------------------
#  set_octave_path('C:/Program Files/GNU Octave/Octave-6.4.0/mingw64/bin/octave20.exe')

## ----addPath,eval=FALSE-------------------------------------------------------
#  
#  add_path('/usr/lib/dynare/matlab')#  Default for Linux
#  
#  add_path('c:/dynare/5.1/matlab') # Default for Windows, but 5.1 can change if later version of
#  # `Dynare` is installed.
#  
#  add_path('/usr/lib/dynare/matlab') # Default for macOS
#  

## ----demo---------------------------------------------------------------------
#  demo(run_dynare)
#  demo(run_models)
#  demo(import_log)


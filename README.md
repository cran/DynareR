DynareR: Bringing the Power of Dynare to R, R Markdown, and Quarto
================

#  <img src="inst/Dynare/DynareR.png" align="right" width="120" />

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/DynareR)](https://cran.r-project.org/package=DynareR)
[![CRAN_Status_Badge](https://cranlogs.r-pkg.org/badges/grand-total/DynareR?color=49C31B)](https://cranlogs.r-pkg.org/badges/grand-total/DynareR?color=49C31B)
[![](https://cranlogs.r-pkg.org/badges/DynareR?color=49C31B)](https://cranlogs.r-pkg.org/badges/DynareR?color=49C31B)

<!-- badges: end -->

# About the Author

The author of this package, **Sagiru Mati**, obtained his PhD in
Economics from the Near East University, North Cyprus. He works at the
Department of Economics, Yusuf Maitama Sule (Northwest) University,
Kano, Nigeria. Please visit his [website](https://smati.com.ng) for more
details.

Please follow his publications on [**ORCID: 0000-0003-1413-3974**](https://orcid.org/0000-0003-1413-3974)

# 1 About DynareR

DynareR is an R package that can run `Dynare` program from R Markdown.

# 2 Requirements

Users need the following in order to knit this document:

1.  Dynare 4.6.1 or above

2.  Octave 5.2.0 or above

3.  Dynare is installed in the standard location as follows:

-   `/usr/lib/dynare/matlab` for `Linux`

-   `/usr/lib/dynare/matlab` for `macOS`

-   `c:/dynare/x.y/matlab` for `Windows`, where `x.y` is `Dynare`
    version number.

If `dynare` and `Octave` are installed in standard location, `DynareR`
package will take care of the configurations, which include adding
`matlab` directory to path, using the latest installed `dynare` and so
on. Otherwise, users have to specify the `matlab` folder using
`add_path` function, set the `Octave` path using the `set_octave_path`
function, or set `dynare` version using the `set_dynare_version`
function.

# 3 Installation

DynareR can be installed using the following commands in R.

``` r
install.packages("DynareR")

          OR
          
devtools::install_github('sagirumati/DynareR')
```

# 4 Usage

Please load the DynareR package as follows:

    ```{r DynareR}                                                             
    library(DynareR)
    ```

Then create a chunk for `dynare` (adopted from Dynare example file
`bkk`) as shown below:

    ```{dynare bkk,eval=T} 
    /*
     * This file implements the multi-country RBC model with time to build,
     * described in Backus, Kehoe and Kydland (1992): "International Real Business
     * Cycles", Journal of Political Economy, 100(4), 745-775.
     *
     * The notation for the variable names are the same in this file than in the paper.
     * However the timing convention is different: we had to taken into account the
     * fact that in Dynare, if a variable is denoted at the current period, then
     * this variable must be also decided at the current period.
     * Concretely, here are the differences between the paper and the model file:
     * - z_t in the model file is equal to z_{t+1} in the paper
     * - k_t in the model file is equal to k_{t+J} in the paper
     * - s_t in the model file is equal to s_{J,t}=s_{J-1,t+1}=...=s_{1,t+J-1} in the paper
     *
     * The macroprocessor is used in this file to create a loop over countries.
     * Only two countries are used here (as in the paper), but it is easy to add
     * new countries in the corresponding macro-variable and completing the
     * calibration.
     *
     * The calibration is the same than in the paper. The results in terms of
     * moments of variables are very close to that of the paper (but not equal
     * since the authors a different solution method).
     *
     * This implementation was written by Sebastien Villemot. Please note that the
     * following copyright notice only applies to this Dynare implementation of the
     * model.
     */

    /*
     * Copyright (C) 2010 Dynare Team
     *
     * This file is part of Dynare.
     *
     * Dynare is free software: you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation, either version 3 of the License, or
     * (at your option) any later version.
     *
     * Dynare is distributed in the hope that it will be useful,
     * but WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     * GNU General Public License for more details.
     *
     * You should have received a copy of the GNU General Public License
     * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
     */

    @#define countries = [ "H", "F" ]
    @#define J = 4

    @#for co in countries
    var C_@{co} L_@{co} N_@{co} A_@{co} K_@{co} Z_@{co} X_@{co} LAMBDA_@{co} S_@{co} NX_@{co} Y_@{co};

    varexo E_@{co};

    parameters beta_@{co} alpha_@{co} eta_@{co} mu_@{co} gamma_@{co} theta_@{co} nu_@{co} sigma_@{co} delta_@{co} phi_@{co} psi_@{co} rho_@{co}_@{co};
    @#endfor

    // Lagrange multiplier of aggregate constraint
    var LGM;

    parameters rho_@{countries[1]}_@{countries[2]} rho_@{countries[2]}_@{countries[1]};

    model;
    @#for co in countries

    Y_@{co} = ((LAMBDA_@{co}*K_@{co}(-@{J})^theta_@{co}*N_@{co}^(1-theta_@{co}))^(-nu_@{co}) + sigma_@{co}*Z_@{co}(-1)^(-nu_@{co}))^(-1/nu_@{co});
    K_@{co} = (1-delta_@{co})*K_@{co}(-1) + S_@{co};
    X_@{co} =
    @# for lag in (-J+1):0
              + phi_@{co}*S_@{co}(@{lag})
    @# endfor
    ;

    A_@{co} = (1-eta_@{co})*A_@{co}(-1) + N_@{co};
    L_@{co} = 1 - alpha_@{co}*N_@{co} - (1-alpha_@{co})*eta_@{co}*A_@{co}(-1);

    // Utility multiplied by gamma
    # U_@{co} = (C_@{co}^mu_@{co}*L_@{co}^(1-mu_@{co}))^gamma_@{co};

    // FOC with respect to consumption
    psi_@{co}*mu_@{co}/C_@{co}*U_@{co} = LGM;

    // FOC with respect to labor
    // NOTE: this condition is only valid for alpha = 1
    psi_@{co}*(1-mu_@{co})/L_@{co}*U_@{co}*(-alpha_@{co}) = - LGM * (1-theta_@{co})/N_@{co}*(LAMBDA_@{co}*K_@{co}(-@{J})^theta_@{co}*N_@{co}^(1-theta_@{co}))^(-nu_@{co})*Y_@{co}^(1+nu_@{co});

    // FOC with respect to capital
    @# for lag in 0:(J-1)
     +beta_@{co}^@{lag}*LGM(+@{lag})*phi_@{co}
    @# endfor
    @# for lag in 1:J
     -beta_@{co}^@{lag}*LGM(+@{lag})*phi_@{co}*(1-delta_@{co})
    @# endfor
     = beta_@{co}^@{J}*LGM(+@{J})*theta_@{co}/K_@{co}*(LAMBDA_@{co}(+@{J})*K_@{co}^theta_@{co}*N_@{co}(+@{J})^(1-theta_@{co}))^(-nu_@{co})*Y_@{co}(+@{J})^(1+nu_@{co});

    // FOC with respect to stock of inventories
     LGM=beta_@{co}*LGM(+1)*(1+sigma_@{co}*Z_@{co}^(-nu_@{co}-1)*Y_@{co}(+1)^(1+nu_@{co}));

    // Shock process
    @# if co == countries[1]
    @#  define alt_co = countries[2]
    @# else
    @#  define alt_co = countries[1]
    @# endif
     (LAMBDA_@{co}-1) = rho_@{co}_@{co}*(LAMBDA_@{co}(-1)-1) + rho_@{co}_@{alt_co}*(LAMBDA_@{alt_co}(-1)-1) + E_@{co};


    NX_@{co} = (Y_@{co} - (C_@{co} + X_@{co} + Z_@{co} - Z_@{co}(-1)))/Y_@{co};

    @#endfor

    // World ressource constraint
    @#for co in countries
      +C_@{co} + X_@{co} + Z_@{co} - Z_@{co}(-1)
    @#endfor
        =
    @#for co in countries
      +Y_@{co}
    @#endfor
        ;

    end;

    @#for co in countries
    beta_@{co} = 0.99;
    mu_@{co} = 0.34;
    gamma_@{co} = -1.0;
    alpha_@{co} = 1;
    eta_@{co} = 0.5; // Irrelevant when alpha=1
    theta_@{co} = 0.36;
    nu_@{co} = 3;
    sigma_@{co} = 0.01;
    delta_@{co} = 0.025;
    phi_@{co} = 1/@{J};
    psi_@{co} = 0.5;
    @#endfor

    rho_H_H = 0.906;
    rho_F_F = 0.906;
    rho_H_F = 0.088;
    rho_F_H = 0.088;

    initval;
    @#for co in countries
    LAMBDA_@{co} = 1;
    NX_@{co} = 0;
    Z_@{co} = 1;
    A_@{co} = 1;
    L_@{co} = 0.5;
    N_@{co} = 0.5;
    Y_@{co} = 1;
    K_@{co} = 1;
    C_@{co} = 1;
    S_@{co} = 1;
    X_@{co} = 1;

    E_@{co} = 0;
    @#endfor

    LGM = 1;
    end;

    shocks;
    var E_H; stderr 0.00852;
    var E_F; stderr 0.00852;
    corr E_H, E_F = 0.258;
    end;

    steady;
    check;

    stoch_simul(order=1, hp_filter=1600);
    ```  

The above chunk creates a Dynare program with the chunk’s content, then
automatically run Dynare, which will save Dynare outputs in the current
directory.

Please note that DynareR uses the chunk name as the model name. So, the
outpus of Dynare are saved in a folder with its respective chunk name.
Thus a new folder `bkk/` will be created in your current working
directory.

By default, `dynare` chunk imports log output as a list of dataframes,
which can be accessed via `dynare$modelName`. Therefore to access the
outputs of the `bkk` model produced by the `dynare` chunk, use
`dynare$bkk`.

Use inline code `` `r dynare$bkk$moments[2,3]` `` to access the value of
second row and third column of the `moments`, which is 0.0024.

# 5 Plotting the IRF

The Impulse Response Function (IRF) is saved by default in
`bkk/bkk/graphs/` folder with the IRF’s name `bkk_IRF_E_H2.pdf`, where
`bkk` is the Dynare model’s name. Therefore, you need to add
`stoch_simul(graph_format = (pdf))` to change the default saving
behaviour of `Dynare` from `eps` to `pdf`.

# 6 DynareR functions for base R

The DynareR package is also designed to work with base R. The following
functions show how to work with DynareR outside the R Markdown or Quarto
documents.

## 6.1 The include_IRF function

Use this function to embed the graphs Impulse Response Function (IRF) in
R Markdown or Quarto document.

The Impulse Response Function (IRF) of the `bkk` model can be fetched
using the following R chunk. Note that only the last part of the IRF’s
name (`E_H2`) is needed, that is `bkk_IRF_` is excluded. Also note that
`out.extra='trim={0cm 7cm 0cm 7cm},clip'` is used to trim the white
space above and below the IRF.

    ```{r IRF,out.extra='trim={0cm 7cm 0cm 7cm},clip',fig.cap="Another of figure generated from Dynare software"} 
    include_IRF("bkk","E_H2")

    # Alternatively, use the path argument 

    ```

``` r
include_IRF(model="bkk",IRF = "E_H2")

# Alternatively, use the path argument 

include_IRF(path="bkk/bkk/graphs/bkk_IRF_E_H2.pdf")
```

However, Dynare figure can only be dynamically included if the output
format is pdf as Dynare produces pdf and eps graphs only.

## 6.2 The write_dyn function

This function writes a new `dyn` file.

Use `write_dyn(code="code",model="someModel")` if you want the `Dynare`
file to live in the current working directory. Use
`write_dyn(code="code",model="path/to/someDirectory/someModel")` if you
want the Dynare file to live in the path different from the current
working directory.

``` r
dynareCodes='var y, c, k, a, h, b;
varexo e, u;
parameters beta, rho, alpha, delta, theta, psi, tau;
alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;
phi   = 0.1;
model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
          *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;
initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul;'


write_dyn(code=dynareCodes, model="example1")

write_dyn(code=dynareCodes,model="DynareR/write_dyn/example1")
```

## 6.3 The write_mod function

This function writes a new `mod` file.

Use `write_mod(code="code",model="someModel")` if you want the `Dynare`
file to live in the current working directory. Use
`write_mod(code="code",model="path/to/someDirectory/someModel")` if you
want the Dynare file to live in the path different from the current
working directory.

``` r
DynareCodes='var y, c, k, a, h, b;
varexo e, u;
parameters beta, rho, alpha, delta, theta, psi, tau;
alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;
phi   = 0.1;
model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
          *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;
initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul;'


write_mod(model="example1",code=dynareCodes)

write_mod(code=dynareCodes,model="DynareR/write_mod/example1")
```

## 6.4 The run_dynare function

Create and run Dynare `mod` file

Use this function to create and run Dynare mod file. Use
`run_dynare(code="code",model="someModel")` if you want the Dynare files
to live in the current working directory. Use
`run_dynare(code="code",model="path/to/someDirectory/someModel")` if you
want the Dynare files to live in the path different from the current
working directory. Use `import_log=T` argument to return the `dynare`
log file as list of dataframes in an environment `dynare`, which can be
accessed via `dynare$modelName`.

``` r
DynareCodes='var y, c, k, a, h, b;
varexo e, u;
parameters beta, rho, alpha, delta, theta, psi, tau;
alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;
phi   = 0.1;
model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
          *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;
initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul;'

run_dynare(code=DynareCodes,model="example1",import_log = T)
run_dynare(code=DynareCodes,model="DynareR/run_dynare/example1")
```

## 6.5 The run_models function

Run multiple existing `mod` or `dyn` files.

Use this function to execute multiple existing Dynare files. Use
`run_models(model="someModel")` if the Dynare files live in the current
working directory. Use
`run_models(model="path/to/someDirectory/someModel")` if the Dynare
files live in the path different from the current working directory. Use
`run_models()` to exectute all the `dynare` models in the current
working directory. Use `run_models("path/to/someDirectory*)` to run all
the `dynare` models in `path/to/someDirectory`.

Where `agtrend.mod`, `bkk.mod` and `example1.mod` are the Dynare model
files (with `mod` or `dyn` extension), which live in the current working
directory.

``` r
demo(agtrend)
demo(bkk)
demo(example1)

# Provide the list of the `Dynare` files in a vector
# Ensure that "agtrend.mod", "bkk.mod" and "example1.mod"
# live in the current working directory

# Copy the dynare files to the current working directory

lapply(c("agtrend","bkk","example1"),\(x) file.copy(paste0(x,"/",x,".mod"),"."))

run_models(c("agtrend","bkk","example1")) # Run the models in the vector.
```

To run all `Dynare` models that live in the current working directory,
use the following:

``` r
run_models() # Run all models in Current Working Directory.
```

To run all `Dynare` models that live in particular path (for example
‘DynareR/run_dynare/’ folder), use the following:

``` r
# Copy the dynare files to the 'DynareR/run_dynare' directory

lapply(c("agtrend","bkk","example1"),\(x) file.copy(paste0(x,".mod"),"DynareR/run_dynare"))

run_models(model = 'DynareR/run_dynare*') # notice the * at the end
```

# 7 import_log function

This function returns the `dynare` log output as a list of dataframes,
which include `summary`, `shocks`, `policy`, `moments`, `decomposition`,
`correlation` and `autocorrelation`. The list is accessible via
`dynare$modelName`. if the model name is `bkk`, the policy variables can
be obtained via `dynare$bkk$policy` as a dataframe.

``` r
import_log(model="bkk")

import_log(path="bkk/bkk.log")

knitr::kable(dynare$bkk$autocorrelation) 
```

# 8 set_dynare_version function

On Windows, you can set the version of dynare you want to use. By
default, `DynareR` package does this for you if the dynare version
ranges from 4.6.1 to 9.9. However, if you are using the development
version of `dynare`, for example version
`6-unstable-2022-04-03-0800-700a0e3a`, you can override the default as
follows

``` r
set_dynare_version("6-unstable-2022-04-03-0800-700a0e3a")
```

# 9 set_octave_path function

You can use this function if `Octave` is not installed in the standard
location

``` r
set_octave_path('C:/Program Files/GNU Octave/Octave-6.4.0/mingw64/bin/octave20.exe')
```

# 10 add_path function

This function is a wrapper of `addpath` in `Octave`. If `dynare` is not
installed in the standard location, use this function to add the
`matlab` subdirectory. By default, `DynareR` does this for if `dynare`
is installed in the standard location.

``` r
add_path('/usr/lib/dynare/matlab')#  Default for Linux

add_path('c:/dynare/5.1/matlab') # Default for Windows, but 5.1 can change if later version of
# `Dynare` is installed.

add_path('/usr/lib/dynare/matlab') # Default for macOS
```

# 11 Demo

The demo files are included and can be accessed via
demo(package=“DynareR”)

``` r
demo(run_dynare)
demo(run_models)
demo(import_log)
```

# 12 Template

Template for R Markdown is created. Go to
`file->New File->R Markdown-> From Template->DynareR`.

# Similar Packages

Similar packages include
[EviewsR](https://github.com/sagirumati/EviewsR) (Mati 2022b, 2020b,Mati,
Civcir, and Abba 2023),
[gretlR](https://github.com/sagirumati/gretlR) (Mati 2020c, 2022c), and
[URooTab](https://github.com/sagirumati/URooTab) (Mati 2023b, 2023a)

For further details, consult Mati (2020a) and Mati (2022a).

<br><br><br><br>


Please download the example files from
[Github](https://github.com/sagirumati/DynareR/tree/master/inst/examples/).


# References

Mati, Sagiru. 2020a. “DynareR: Bringing the Power of Dynare to
<span class="nocase">R, R Markdown, and Quarto</span>.” *CRAN*.
<https://CRAN.R-project.org/package=DynareR>.

———. 2020b. *EviewsR: A Seamless Integration of EViews and R*.
<https://CRAN.R-project.org/package=EviewsR>.

———. 2020c. *gretlR: A Seamless Integration of Gretl and R*.
<https://CRAN.R-project.org/package=gretlR>.

———. 2021. “Do as Your Neighbours Do? Assessing the Impact of Lockdown
and Reopening on the Active COVID-19 Cases in Nigeria.” *Social Science
&Amp; Medicine* 270 (February): 113645.
<https://doi.org/10.1016/j.socscimed.2020.113645>.

———. 2022a. “Package ‘DynareR’.”
<https://CRAN.R-project.org/package=DynareR/DynareR.pdf>.

———. 2022b. “Package ‘EviewsR’.”
<https://CRAN.R-project.org/package=EviewsR/EviewsR.pdf>.

———. 2022c. “Package ‘gretlR’.”
<https://CRAN.R-project.org/package=gretlR/gretlR.pdf>.

———. 2023a. “Package ‘URooTab’.”
<https://CRAN.R-project.org/package=URooTab/URooTab.pdf>.

———. 2023b. *URooTab: Tabular Reporting of EViews Unit Root Tests*.
<https://github.com/sagirumati/URooTab>.

Mati, Sagiru, Irfan Civcir, and S. I. Abba. 2023. “EviewsR: An r Package
for Dynamic and Reproducible Research Using EViews, r, r Markdown and
Quarto.” *The R Journal* 15 (2): 169–205.
<https://doi.org/10.32614/rj-2023-045>.

Mati, Sagiru, Irfan Civcir, and Hüseyin Ozdeser. 2019. “ECOWAS COMMON
CURRENCY: HOW PREPARED ARE ITS MEMBERS?” *Investigación Económica* 78
(308): 89. <https://doi.org/10.22201/fe.01851667p.2019.308.69625>.

Mati, Sagiru, Irfan Civcir, and Hüseyin Özdeşer. 2023. “ECOWAS Common
Currency, a Mirage or Possibility?” *Panoeconomicus* 70 (2): 239–60.
<https://doi.org/10.2298/pan191119015m>.

Mati, Sagiru, Magdalena Radulescu, Najia Saqib, Ahmed Samour, Goran
Yousif Ismael, and Nazifi Aliyu. 2023. “Incorporating Russo-Ukrainian
War in Brent Crude Oil Price Forecasting: A Comparative Analysis of
ARIMA, TARMA and ENNReg Models.” *Heliyon* 9 (11): e21439.
<https://doi.org/10.1016/j.heliyon.2023.e21439>.

# Dynare 0.1.5

* Added support for Matlab. The package now works with both Octave and Matlab

* Added ``set_matlab_path()` function, which sets path to the Matlab executable to use.


* include_IRF() function gained a new argument `crop=TRUE`. This crops the white space around IRF graphs

* The package now automatically converts IRF graphs to `png` if the output format is not latex

# Dynare 0.1.4

What is new? 

* New function `input_tex()` has been added. It can be used to include raw LaTeX file into R Markdown or Quarto.

* Added new funciton `add_matlab_path()` to replace ``add_path()` function, which is still available for backward compatibility.

* The dependency on kableExtra has been removed

* Updated Octave search path

* The package now chooses the `Octave` version  compatible with installed Dynare, or with the dynare version chosen via `set_dynare_version()` function

* Updated `set_dynare_version()` function

* Author's email changed to <sagirumati@gmail.com>

* Package Title changed to "Bringing the Power of Dynare to R, R Markdown, and Quarto"

* Startup message with citation information


# Dynare 0.1.3

What is new? 

* Bug fixes: updated `dynare` configurations for all platforms

* New function `import_log`, `set_dynare_version`, `set_octave_path` and `add_path`.

* `set_dynare_version`, `set_octave_path` and `add_path` functions are only used if `Dynare` and `Octave` are not installed in the standard location.

* `run_models` and `run_dynare` functions gained new argument `import_log=F`. Whether to import `Dynare` log file as a list of dataframes (`FALSE` by default)

* `path` argument has been dropped in `write_dyn`, `write_mod`, `run_models` and `run_dynare` functions. You can now include the path in the `model` argument. For example `model='someModel'` for model in the current working directory, `model='some/path/to/the/somemodel'` for model in `some/path/to/the/` directory.

* The default directory of `dynare` chunk is the current working directory.

* Updated demo

* Updated Vignettes

* Updated R Markdown template

# Dynare 0.1.2

What is new? 

* DynareR is now platform independent. It works on the major operating systems


# Dynare 0.1.1

What's new?

* Template for R Markdown is created. Go to `file->New File->R Markdown-> From Template->DynareR`.

* DynareR can be used with both base R and R Markdown

* New functions `include_IRF`, `run_dynare`, `run_models`, `write_dyn` and `write_mod` are created

* Embed the graph Inpulse Response Function in R Markdown with `include_IRF`

* Demo files are accessible via demo(package="DynareR")


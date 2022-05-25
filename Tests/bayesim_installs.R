
# only install packages once
required_packages <- c("tidyr", "roxygen2", "posterior", "latex2exp", "ggthemes", "dplyr", "stats", "matrixStats", "ggplot2")
print("To run Bayesim, you need a couple of packages. Those are:")
print(required_packages)
print("If you like to proceed, type yes / no")
response <- readline()
if(grepl("yes", tolower(response), fixed = TRUE)) {
  print("Install required packages...")
  for(package in required_packages) {
    if(!requireNamespace(package)) {
      install.packages(package, dependencies = TRUE)
    }
  }
} else {
  print("Packages will not be installed.")
  print("Given, without those packages, this cannot be run, the script stops now.")
  print("Please Re-Run the script, if you wish to still install the packages!")
  stop("END INSTALL SCRIPT PREMATURELY")
}

suggested_packages <- c("cmdstanr", "foreach", "parallel", "doParallel")
print("For Bayesim, some packages are suggested, but not required. Those are:")
print(suggested_packages)
print("If you like to proceed, type yes / no")
response <- readline()
if(grepl("yes", tolower(response), fixed = TRUE)) {
  print("Install suggested packages...")
  for(package in suggested_packages) {
    if(!requireNamespace(package)) {
      install.packages(package, dependencies = TRUE)
    }
  }
} else {
  print("Packages will not be installed.")
}

# ask before installing test-packages and only install it once
required_test_packages <- c("extraDistr", "testthat")
print("Do you want to install the packages required to run Unit-tests? Those are:")
print(required_test_packages)
print("If so, type: yes / no")
response <- readline()
if(grepl("yes", tolower(response), fixed = TRUE)) {
  print("Install test packages...")
  for(package in required_test_packages) {
    if(!requireNamespace(package)) {
      install.packages(package, dependencies = TRUE)
    }
  }
} else {
  print("Packages will not be installed. Re-run script, if packages are required.")
}

# ask, before installing (given remotes is not required for normal operation)
# are only installed once, if a new version was available
required_github_package <- c("stan-dev/cmdstanr", "paul-buerkner/brms")
print("Do you want to install the reuqired github-packages? Those are:")
print(required_github_package)
print("For this, the remotes package will be installed as well. If that is ok, type: yes / no")
response <- readline()
if(grepl("yes", tolower(response), fixed = TRUE)) {
  print("Install github packages...")
  if(!requireNamespace("remotes")) {
    install.packages("remotes", dependencies = TRUE)
  }
  for(package in required_github_package) {
    remotes::install_github(package)
  }
} else {
  print("Packages will not be installed. Re-run script, if packages are required.")
}

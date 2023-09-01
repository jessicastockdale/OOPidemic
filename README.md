
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OOPidemic

<!-- badges: start -->

<!-- badges: end -->

The goal of OOPidemic is to simulate an disease epidemic using an object
oriented approach.

## Installation

You can install the development version of OOPidemic like so:

``` r
install.packages("devtools")
library(devtools)
# This repro is private so generate a personal access token in github to allow you to install it
# read the auth_token parameter details in the install_github docs for more details
?install_github 
install_github("snailvet/OOPidemic", auth_token = [Your personal access token])
```

## Example

This is a basic example which shows you how to run a simple simluation:

``` r
library(OOPidemic)

# set up a reference strain with a randomised genome
ref_strain <- ReferenceStrain$new(
  name = "ref_strain",
  g_len = 1000
)

# set up a group 
group <- Group$new(
  id = 1,
  ref_strain = ref_strain
)

# set up a lab to take samples 
lab <- Lab$new()

# run simluation
group$run_simulation(lab)
#> Group: 
#>    Id:                 1
#>    Time:               10
#>    Group Size:         100
#>    Susceptible Hosts:  26
#>    Exposed Hosts:      29
#>    Infectious Hosts:   37
#>    Recovered Hosts:    8
#> Group: 
#>    Id:                 1
#>    Time:               20
#>    Group Size:         100
#>    Susceptible Hosts:  0
#>    Exposed Hosts:      0
#>    Infectious Hosts:   11
#>    Recovered Hosts:    89
#> Group: 
#>    Id:                 1
#>    Time:               25
#>    Group Size:         100
#>    Susceptible Hosts:  0
#>    Exposed Hosts:      0
#>    Infectious Hosts:   0
#>    Recovered Hosts:    100
#> Lab: 
#>    WG Sequences:   100
```

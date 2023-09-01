
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OOPidemic

<!-- badges: start -->

<!-- badges: end -->

The goal of OOPidemic is to simulate an disease epidemic using an object
oriented approach.

## Installation

You can install the development version of OOPidemic like so:

``` r
# install.packages("devtools")
devtools::install_github("snailvet/OOPidemic")
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
#>    Susceptible Hosts:  29
#>    Exposed Hosts:      30
#>    Infectious Hosts:   37
#>    Recovered Hosts:    4
#> Group: 
#>    Id:                 1
#>    Time:               20
#>    Group Size:         100
#>    Susceptible Hosts:  0
#>    Exposed Hosts:      0
#>    Infectious Hosts:   12
#>    Recovered Hosts:    88
#> Group: 
#>    Id:                 1
#>    Time:               27
#>    Group Size:         100
#>    Susceptible Hosts:  0
#>    Exposed Hosts:      0
#>    Infectious Hosts:   0
#>    Recovered Hosts:    100
#> Lab: 
#>    WG Sequences:   100
```

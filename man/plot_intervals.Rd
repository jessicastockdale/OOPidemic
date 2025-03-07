% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-intervals.R
\name{plot_intervals}
\alias{plot_intervals}
\title{Plot a histogram of the desired interval in an population}
\usage{
plot_intervals(
  x,
  interval = "transmission",
  show_distribution = FALSE,
  combine = FALSE,
  group_labels = NULL
)
}
\arguments{
\item{x}{A Population or Group object that has had it's outbreak simulated}

\item{interval}{The interval that you want plotted (default is the transmision interval)}

\item{show_distribution}{Logical. If \code{TRUE} and parameters were provided for the interval being plotted then the expected distribution will be overlaid on top of the histogram.}

\item{combine}{A boolean indicating if bins for different groups should be combined (only has an effect if a Population object with multiple groups is provided). If combine is \code{TRUE} and \code{x} is a Population object then \code{show_distribution} will be overridden and set to \code{FALSE}.}

\item{group_labels}{A vector of strings used to label multiple groups. Should be in the same order as the groups in the Population object. (similarly to combine, this only has an effect if a population with multiple groups is provided)}
}
\value{
A ggplot object
}
\description{
Plot a histogram of the desired interval in an population
}
\examples{
# An SEIR model outbreak in a single homogenous group using the 
# serial interval method to find infectors
ref_strain <- ReferenceStrain$new("ref_strain")
group <- Group$new(
     1, 
     ref_strain, 
     inf_rate = 0.75, 
     find_infector_method = "serial", 
     si_shape = 6, si_rate = 1
)
group$run_simulation(Lab$new())
figure <- plot_intervals(group, interval = "serial", show_distribution = TRUE)
plot(figure)

# set.seed(1)
# An SEIR model outbreak in a population with two groups using the tranmission method 
# for finding infectors.
# The two groups will have different transmission interval shapes
population <- Population$new(
    c(
        Group$new(1, ref_strain, init_inf = 5, init_sus = 95),
        Group$new(2, ref_strain, init_inf = 0, init_sus = 100, trans_int_shape = 5)
    ),
    matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
    Lab$new()
)
population$run_simulation()
fig_1 <- plot_intervals(
     population, 
     show_distribution = TRUE, group_labels = c("Canine", "Feline")
)
plot(fig_1)

fig_2 <- plot_intervals(population, combine = TRUE, group_labels = c("Canine", "Feline"))
plot(fig_2)

fig_3 <- plot_intervals(
     population, 
     show_distribution = TRUE, interval = "generation", group_labels = c("Canine", "Feline")
)
plot(fig_3)

}

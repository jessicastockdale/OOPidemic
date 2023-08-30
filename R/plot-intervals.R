###############################################################################
# interval functions

incubation_time <- function(host) {
    # calculate the incubation time for this host 

    if (is.na(host$exposure_time)) {
        # susceptible host
        return(NA)
    } else { 
        return(host$infectious_time - host$exposure_time)
    } 
}

infectiousness_duration <- function(host) {
    # calculate the time this host was infectious 

    if (is.na(host$exposure_time)) {
        # susceptible host
        return(NA)
    } else { 
        return(host$recovery_time - host$infectious_time)
    } 
}

generation_time <- function(host) {
    # calculate the generation time for this host and their infector

    if (is.na(host$exposure_time)) {
        # susceptible host
        return(NA)
    } else if (host$is_index) {
         # ignore index cases as we can't know their generation time
        return(NA)
    } else {
        return(host$exposure_time - host$infector$exposure_time)
    }
}

serial_interval <- function(host) {
    # calculate the serial interval for this host and their infector

    if (is.na(host$exposure_time)) {
        # susceptible host
        return(NA)
    } else if (host$is_index) {
         # ignore index cases as we can't know their serial interval
        return(NA)
    } else {
        return(host$infectious_time - host$infector$infectious_time)
    }
}

transmission_interval <- function(host) {
    # calculate the serial interval for this host and their infector

    if (is.na(host$exposure_time)) {
        # susceptible host
        return(NA)
    } else if (host$is_index) {
         # ignore index cases as we can't know their transmission interval
        return(NA)
    } else {
        return(host$exposure_time - host$infector$infectious_time)
    }
}

###############################################################################

add_expected_distribution <- function(population, figure, interval) {

    shape <- NULL

    # get gamma parameters
    if (interval == "incubation") {
        shape <- population$inc_shape
        rate <- population$inc_rate

    } else if (interval == "infectious") {
        shape <- population$rec_shape
        rate <- population$rec_rate

    } else if (all(
        interval == "serial",
        population$find_infector_method == "serial"
    )) {
        shape <- population$si_shape
        rate <- population$si_rate
        
    } else if (all(
        interval == "transmission",
        population$find_infector_method == "transmission"
    )) {
        shape <- population$trans_int_shape
        rate <- population$trans_int_rate

    } else if (all(
        interval == "generation",
        population$find_infector_method == "generation"
    )) {
        shape <- population$gen_shape
        rate <- population$gen_rate
        
    }

    if (!is.null(shape)) {
        figure <- figure + ggplot2::stat_function(
            fun = stats::dgamma, 
            args = list(shape = shape, rate = rate),
            geom = "line",
            ggplot2::aes(linetype = as.character(population$id))
        )

    }
    
    return(figure)
}

###############################################################################
# plot functions 

#' Plot a histogram of the desired interval in an outbreak
#'
#' @param x An Outbreak or Population object that has had it's outbreak simulated
#' @param interval The interval that you want plotted (default is the transmision interval)
#' @param show_distribution Logical. If `TRUE` and parameters were provided for the interval being plotted then the expected distribution will be overlaid on top of the histogram. 
#' @param combine A boolean indicating if bins for different populations should be combined (only has an effect if an Outbreak object with multiple populations is provided). If combine is `TRUE` and `x` is an Outbreak object then `show_distribution` will be overridden and set to `FALSE`.
#' @param population_labels A vector of strings used to label multiple populations. Should be in the same order as the populations in the Outbreak object. (similarly to combine, this only has an effect if an outbreak with multiple populations is provided)
#' 
#' @importFrom stats density
#' @importFrom dplyr count
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # An SEIR model outbreak in a single homogenous population using the 
#' # serial interval method to find infectors
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' population <- Population$new(
#'      1, 
#'      ref_strain, 
#'      inf_rate = 0.75, 
#'      find_infector_method = "serial", 
#'      si_shape = 6, si_rate = 1
#' )
#' population$run_simulation(Lab$new())
#' figure <- plot_intervals(population, interval = "serial", show_distribution = TRUE)
#' plot(figure)
#' 
#' # set.seed(1)
#' # An SEIR model outbreak with two populations using the tranmission method 
#' # for finding infectors.
#' # The two populations will have different transmission interval shapes
#' outbreak <- Outbreak$new(
#'     c(
#'         Population$new(1, ref_strain, init_inf = 5, init_sus = 95),
#'         Population$new(2, ref_strain, init_inf = 0, init_sus = 100, trans_int_shape = 5)
#'     ),
#'     matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
#'     Lab$new()
#' )
#' outbreak$run_simulation()
#' fig_1 <- plot_intervals(
#'      outbreak, 
#'      show_distribution = TRUE, population_labels = c("Canine", "Feline")
#' )
#' plot(fig_1)
#' 
#' fig_2 <- plot_intervals(outbreak, combine = TRUE, population_labels = c("Canine", "Feline"))
#' plot(fig_2)
#' 
#' fig_3 <- plot_intervals(
#'      outbreak, 
#'      show_distribution = TRUE, interval = "generation", population_labels = c("Canine", "Feline")
#' )
#' plot(fig_3)
#' 
plot_intervals <- function(
    x, 
    interval = "transmission", 
    show_distribution = FALSE,
    combine = FALSE,
    population_labels = NULL
) {

    # validations for x
    stopifnot(R6::is.R6(x))
    if ("Population" %in% class(x)) {
        populations <- c(x)
    } else if ("Outbreak" %in% class(x)) {
        populations <- x$populations
        if (all(!is.null(population_labels), length(populations) != length(population_labels))) {
            stop("The length of `population_labels` must equal the number of populations in x")
        }

    } else {
        stop("`x` must be an Outbreak or Population object")
    }

    # setup interval function and title
    if (interval == "serial") {
        interval_function <- serial_interval
        title <- "Serial Interval"
    } else if (interval == "generation") {
        interval_function <- generation_time
        title <- "Generation Time"
    } else if (interval == "incubation") {
        interval_function <- incubation_time
        title <- "Duration of Incubation"
    } else if (interval == "infectious") {
        interval_function <- infectiousness_duration
        title <- "Duration of Infectiousness"
    } else if (interval == "transmission") {
        interval_function <- transmission_interval
        title <- "Transmission Interval"
    } else {
        # interval validation
        stop(paste0("interval = ", interval, " not a valid option. \n Please use one of: 'serial', 'generation', 'incubation', 'infectious', or transmission"))
    }

    # other validations
    stopifnot(is.logical(show_distribution), length(show_distribution) == 1)
    stopifnot(is.logical(combine), length(combine) == 1)
    stopifnot(any(is.null(population_labels), is.character(population_labels)))

    # get a vector with all hosts
    hosts <- unlist(
        lapply(populations, function(pop) pop$hosts),
        recursive = FALSE
    )

    # get intervals for hosts (susceptible hosts and index hosts will return an NA)
    data <- data.frame(
        pop_id = unlist(lapply(populations, function(pop) rep(as.character(pop$id), pop$size)), recursive = FALSE),
        value = vapply(hosts, interval_function, numeric(1L))
    )

    # remove NAs
    data <- data[!is.na(data$value), ]

    # if we are combining populations then it's easiest to set pop_id to the same thing
    if (combine) {data$pop_id <- 1}

    # # plot histogram
    figure <- ggplot2::ggplot(
        data, 
        ggplot2::aes(x = data$value, fill = data$pop_id)
    ) + 
        ggplot2::geom_histogram(
            ggplot2::aes(y = ggplot2::after_stat(density)), 
            binwidth = 1, 
            position = "dodge"
        )
    
    # add bin sizes 
    figure <- figure + ggplot2::stat_bin(
        ggplot2::aes(y = ggplot2::after_stat(density), label = ggplot2::after_stat(count)), 
        binwidth = 1,
        geom = "text",
        position = "dodge", 
        vjust = - 1
    )  

    # set theme, title, and axis labels
    figure <- figure + ggplot2::theme_minimal() +
        ggplot2::ggtitle(title, paste0("n = ", length(data$value))) +
        ggplot2::labs(x = "Time", y = "Density") 

    # show distributions if wanted and if we aren't combining an Outbreaks's populations
    if (show_distribution) {
        for (population in populations) {
            figure <- add_expected_distribution(population, figure, interval)
        }
    }

    # Do we need to display a legend for populations?
    if (length(unique(data$pop_id)) == 1) {
        # legend isn't needed if there is only one population or we are combining 
        figure <- figure + ggplot2::theme(legend.position = "none")

    } else if (!is.null(population_labels)) {
        # not combining and we have population labels
        figure <- figure + ggplot2::scale_fill_discrete(
            name = "Population",
            labels = population_labels
        )
        figure <- figure + ggplot2::scale_linetype_discrete(
            name = "Expected\n Distributions",
            labels = population_labels
        )

    } else {
        # no population labels
        figure <- figure + ggplot2::scale_fill_discrete(name = "Population")
        figure <- figure + ggplot2::scale_linetype_discrete(name = "Population")
    }


    return(figure)
}

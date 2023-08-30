###############################################################################
# interval functions
# Get intervals from a host

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

#' Add a gamma distribution to a figure
#'
#' @param group A Group object
#' @param figure The figure to add the distribution plot to
#' @param interval The interval being plotted
#' 
#' @return A ggplot object
#'
#' @noRd 
add_expected_distribution <- function(group, figure, interval) {

    shape <- NULL

    # get gamma parameters
    if (interval == "incubation") {
        shape <- group$inc_shape
        rate <- group$inc_rate

    } else if (interval == "infectious") {
        shape <- group$rec_shape
        rate <- group$rec_rate

    } else if (all(
        interval == "serial",
        group$find_infector_method == "serial"
    )) {
        shape <- group$si_shape
        rate <- group$si_rate
        
    } else if (all(
        interval == "transmission",
        group$find_infector_method == "transmission"
    )) {
        shape <- group$trans_int_shape
        rate <- group$trans_int_rate

    } else if (all(
        interval == "generation",
        group$find_infector_method == "generation"
    )) {
        shape <- group$gen_shape
        rate <- group$gen_rate
        
    }

    if (!is.null(shape)) {
        figure <- figure + ggplot2::stat_function(
            fun = stats::dgamma, 
            args = list(shape = shape, rate = rate),
            geom = "line",
            ggplot2::aes(linetype = as.character(group$id))
        )

    }
    
    return(figure)
}

###############################################################################
# plot functions 

#' Plot a histogram of the desired interval in an population
#'
#' @param x A Population or Group object that has had it's outbreak simulated
#' @param interval The interval that you want plotted (default is the transmision interval)
#' @param show_distribution Logical. If `TRUE` and parameters were provided for the interval being plotted then the expected distribution will be overlaid on top of the histogram. 
#' @param combine A boolean indicating if bins for different groups should be combined (only has an effect if a Population object with multiple groups is provided). If combine is `TRUE` and `x` is a Population object then `show_distribution` will be overridden and set to `FALSE`.
#' @param group_labels A vector of strings used to label multiple groups. Should be in the same order as the groups in the Population object. (similarly to combine, this only has an effect if a population with multiple groups is provided)
#' 
#' @importFrom stats density
#' @importFrom dplyr count
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # An SEIR model outbreak in a single homogenous group using the 
#' # serial interval method to find infectors
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' group <- Group$new(
#'      1, 
#'      ref_strain, 
#'      inf_rate = 0.75, 
#'      find_infector_method = "serial", 
#'      si_shape = 6, si_rate = 1
#' )
#' group$run_simulation(Lab$new())
#' figure <- plot_intervals(group, interval = "serial", show_distribution = TRUE)
#' plot(figure)
#' 
#' # set.seed(1)
#' # An SEIR model outbreak in a population with two groups using the tranmission method 
#' # for finding infectors.
#' # The two groups will have different transmission interval shapes
#' population <- Population$new(
#'     c(
#'         Group$new(1, ref_strain, init_inf = 5, init_sus = 95),
#'         Group$new(2, ref_strain, init_inf = 0, init_sus = 100, trans_int_shape = 5)
#'     ),
#'     matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
#'     Lab$new()
#' )
#' population$run_simulation()
#' fig_1 <- plot_intervals(
#'      population, 
#'      show_distribution = TRUE, group_labels = c("Canine", "Feline")
#' )
#' plot(fig_1)
#' 
#' fig_2 <- plot_intervals(population, combine = TRUE, group_labels = c("Canine", "Feline"))
#' plot(fig_2)
#' 
#' fig_3 <- plot_intervals(
#'      population, 
#'      show_distribution = TRUE, interval = "generation", group_labels = c("Canine", "Feline")
#' )
#' plot(fig_3)
#' 
plot_intervals <- function(
    x, 
    interval = "transmission", 
    show_distribution = FALSE,
    combine = FALSE,
    group_labels = NULL
) {

    # validations for x
    stopifnot(R6::is.R6(x))
    if ("Group" %in% class(x)) {
        groups <- c(x)
    } else if ("Population" %in% class(x)) {
        groups <- x$groups
        if (all(!is.null(group_labels), length(groups) != length(group_labels))) {
            stop("The length of `group_labels` must equal the number of groups in x")
        }

    } else {
        stop("`x` must be a Population or Group object")
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
    stopifnot(any(is.null(group_labels), is.character(group_labels)))

    # get a vector with all hosts
    hosts <- unlist(
        lapply(groups, function(group) group$hosts),
        recursive = FALSE
    )

    # get intervals for hosts (susceptible hosts and index hosts will return an NA)
    data <- data.frame(
        group_id = unlist(lapply(groups, function(grp) rep(as.character(grp$id), grp$size)), recursive = FALSE),
        value = vapply(hosts, interval_function, numeric(1L))
    )

    # remove NAs
    data <- data[!is.na(data$value), ]

    # if we are combining groups then it's easiest to set group_id to the same thing
    if (combine) {data$group_id <- 1}

    # # plot histogram
    figure <- ggplot2::ggplot(
        data, 
        ggplot2::aes(x = data$value, fill = data$group_id)
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

    # show distributions if wanted and if we aren't combining an Populations's groups
    if (show_distribution) {
        for (group in groups) {
            figure <- add_expected_distribution(group, figure, interval)
        }
    }

    # Do we need to display a legend for groups?
    if (length(unique(data$group_id)) == 1) {
        # legend isn't needed if there is only one group or we are combining 
        figure <- figure + ggplot2::theme(legend.position = "none")

    } else if (!is.null(group_labels)) {
        # not combining and we have group labels
        figure <- figure + ggplot2::scale_fill_discrete(
            name = "Group",
            labels = group_labels
        )
        figure <- figure + ggplot2::scale_linetype_discrete(
            name = "Expected\n Distributions",
            labels = group_labels
        )

    } else {
        # no group labels
        figure <- figure + ggplot2::scale_fill_discrete(name = "Group")
        figure <- figure + ggplot2::scale_linetype_discrete(name = "Group")
    }


    return(figure)
}

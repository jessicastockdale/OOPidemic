#' Plot S(E)IR compartment dynamics for disease outbreak
#'
#' @param x An Outbreak or Population object that has had it's outbreak simulated
#' @param frequency Integer indicating how frequently the compartment sizes should be measured
#' @param combine A boolean indicating if compartment sizes for different populations should be combined (only has an effect if an Outbreak object with multiple populations is provided) 
#' @param population_labels A vector of strings used to label multiple populations. Should be in the same order as the populations in the Outbreak object. (similarly to combine, this only has an effect if an outbreak with multiple populations is provided)
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # An SIR model outbreak in a single homogenous population
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' population <- Population$new(1, ref_strain, inf_rate = 0.75, inc_shape = 0)
#' lab_pop <- Lab$new()
#' population$run_simulation(lab_pop)
#' plot_compartment_dynamics(population)
#' 
#' # An SEIR model outbreak with two populations
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' outbreak <- Outbreak$new(
#'     c(
#'         Population$new(1, ref_strain, init_inf = 5, init_sus = 95),
#'         Population$new(2, ref_strain, init_inf = 0, init_sus = 100)
#'     ),
#'     matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
#'     Lab$new()
#' )
#' outbreak$run_simulation()
#' # plot each population separately
#' plot_compartment_dynamics(outbreak, population_labels = c("Senior", "Adult"))
#' # plot populations combined
#' plot_compartment_dynamics(outbreak, combine = TRUE)
#' 
plot_compartment_dynamics <- function(
    x,
    frequency = 1,
    combine = FALSE,
    population_labels = NULL
) {

    # R6 validations
    stopifnot(R6::is.R6(x))
    # correct class of x is checked later
    stopifnot(is.numeric(frequency), length(frequency) == 1, frequency %% 1 == 0)
    stopifnot(is.logical(combine), length(combine) == 1)
    stopifnot(any(is.null(population_labels), is.character(population_labels)))

    x_classes <- class(x)
    if ("Population" %in% x_classes) {
        populations <- c(x)

    } else if ("Outbreak" %in% x_classes) {
        populations <- x$populations
        if (all(!is.null(population_labels), length(populations) != length(population_labels))) {
            stop("The length of `population_labels` must equal the number of populations in x")
        }
    } else {
        stop("`x` must be an Outbreak or Population object")
    }

    times <- base::seq.int(0, x$time, frequency)
    data <- data.frame(
        t = numeric(0),
        id = character(0),
        S = numeric(0),
        E = numeric(0),
        I = numeric(0),
        R = numeric(0)
    )
    
    for (time in times) {
        for (population in populations) {
            data <- rbind(data, compartment_sizes(population, time))
        }
    }    

    # check if there is any exposed compartment
    title <- "SEIR Outbreak"
    compartment_labels <- c("Exposed", "Infectious", "Recovered", "Susceptible")
    if (all(data$E == 0)) {
        # there is no exposed compartment
        data$E <- NULL
        title <- "SIR Outbreak"
        compartment_labels <- compartment_labels[2:length(compartment_labels)]
    }

    # combine populations if needed / delete the id column if only a single population
    start_col <- 3
    if (any(combine, length(unique(data$id)) == 1)) {
        # combine data from each population
        data$id <- NULL
        start_col <- start_col - 1
        data <- stats::aggregate(. ~ t, data, FUN = sum)
    }

    # pivot data longer for plotting
    data_long <- tidyr::pivot_longer(data, cols = c(seq(start_col, ncol(data))))

    figure <- generate_figure(data_long, title, compartment_labels, population_labels)

    plot(figure)

    return(figure)
}

#' Get the SEIR compartment sizes for a population at a given time
#'
#' @param population A Population object
#' @param time The time that you want the compartyment sizes for
#'
#' @return A list with: time, population id, and SEIR sizes
#'
#' @noRd 
compartment_sizes <- function(population, time) {

    result <- list(
        t = time,
        id = as.character(population$id),
        S = length(population$susceptible_hosts(time)),
        E = length(population$exposed_hosts(time)),
        I = length(population$infectious_hosts(time)),
        R = length(population$recovered_hosts(time))
    )

    return(result)
}

#' Build the ggplot figure for displaying compartment sizes
#'
#' @param data A dataframe with the following columns: 
#'  * t: time
#'  * id: the population id (optional)
#'  * name: the name of the compartment. Options are S, E, I, and R
#'  * value: the size of the the compartment
#' @param title the title of the figure
#' @param compartment_labels the labels of the compartments (lexicographic order)
#' @param population_labels alternative labels for populations if there is more then one of them. If there is multiple and this is left as NULL then the values in the id column of `data` will be used.
#'
#' @return A ggplot figure
#'
#' @noRd
generate_figure <- function(data, title, compartment_labels, population_labels) {

    if ("id" %in% names(data)) {
        # if there still is an id column then we have multiple populations to plot

        figure <- ggplot2::ggplot(
            data, 
            ggplot2::aes(
                x = t, y = data$value, 
                colour = data$name, 
                linetype = data$id,
                group = interaction(data$name, data$id)
            )
        ) 
        figure <- figure + ggplot2::geom_line(lwd = 1)

    } else {
        figure <- ggplot2::ggplot(
            data, 
            ggplot2::aes(x = t, y = data$value, group = data$name, colour = data$name)
        ) + ggplot2::geom_line(lwd = 2)

    }

    figure <- figure + ggplot2::theme_minimal() 
    figure <- figure + ggplot2::xlab("Time")
    figure <- figure + ggplot2::ylab("Number of individuals")
    figure <- figure + ggplot2::guides(linetype = ggplot2::guide_legend(title = "Population"))
    figure <- figure + ggplot2::scale_colour_discrete(
        name = "Compartment",
        labels = compartment_labels
    ) 
    figure <- figure + ggplot2::ggtitle(title)

    # add population labels if they've been supplied
    if (!is.null(population_labels)) {
        figure <- figure + ggplot2::scale_linetype_discrete(
            name = "Population", 
            labels = population_labels
        )
    }

    return(figure)
}


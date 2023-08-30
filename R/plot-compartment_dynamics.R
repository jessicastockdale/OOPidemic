#' Plot S(E)IR compartment dynamics for disease outbreak
#'
#' @param x A Population or Group object that has had it's outbreak simulated
#' @param frequency Integer indicating how frequently the compartment sizes should be measured
#' @param combine A boolean indicating if compartment sizes for different groups should be combined (only has an effect if an Population object with multiple groups is provided) 
#' @param group_labels A vector of strings used to label multiple groups. Should be in the same order as the groups in the Population object. (similarly to combine, this only has an effect if a population with multiple groups is provided)
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' # An SIR model outbreak in a single homogenous group
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' group <- Group$new(1, ref_strain, inf_rate = 0.75, inc_shape = 0)
#' lab_group <- Lab$new()
#' group$run_simulation(lab_group)
#' plot_compartment_dynamics(group)
#' 
#' # An SEIR model outbreak ina population with two groups
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' population <- Population$new(
#'     c(
#'         Group$new(1, ref_strain, init_inf = 5, init_sus = 95),
#'         Group$new(2, ref_strain, init_inf = 0, init_sus = 100)
#'     ),
#'     matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
#'     Lab$new()
#' )
#' population$run_simulation()
#' # plot each group separately
#' plot_compartment_dynamics(population, group_labels = c("Senior", "Adult"))
#' # plot groups combined
#' plot_compartment_dynamics(population, combine = TRUE)
#' 
plot_compartment_dynamics <- function(
    x,
    frequency = 1,
    combine = FALSE,
    group_labels = NULL
) {

    # R6 validations
    stopifnot(R6::is.R6(x))
    # correct class of x is checked later
    stopifnot(is.numeric(frequency), length(frequency) == 1, frequency %% 1 == 0)
    stopifnot(is.logical(combine), length(combine) == 1)
    stopifnot(any(is.null(group_labels), is.character(group_labels)))

    x_classes <- class(x)
    if ("Group" %in% x_classes) {
        groups <- c(x)

    } else if ("Population" %in% x_classes) {
        groups <- x$groups
        if (all(!is.null(group_labels), length(groups) != length(group_labels))) {
            stop("The length of `group_labels` must equal the number of groups in x")
        }
    } else {
        stop("`x` must be a Population or Group object")
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
        for (group in groups) {
            data <- rbind(data, compartment_sizes(group, time))
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

    # combine groups if needed / delete the id column if only a single group
    start_col <- 3
    if (any(combine, length(unique(data$id)) == 1)) {
        # combine data from each group
        data$id <- NULL
        start_col <- start_col - 1
        data <- stats::aggregate(. ~ t, data, FUN = sum)
    }

    # pivot data longer for plotting
    data_long <- tidyr::pivot_longer(data, cols = c(seq(start_col, ncol(data))))

    figure <- generate_figure(data_long, title, compartment_labels, group_labels)

    plot(figure)

    return(figure)
}

#' Get the SEIR compartment sizes for a group at a given time
#'
#' @param group A Group object
#' @param time The time that you want the compartyment sizes for
#'
#' @return A list with: time, group id, and SEIR sizes
#'
#' @noRd 
compartment_sizes <- function(group, time) {

    result <- list(
        t = time,
        id = as.character(group$id),
        S = length(group$susceptible_hosts(time)),
        E = length(group$exposed_hosts(time)),
        I = length(group$infectious_hosts(time)),
        R = length(group$recovered_hosts(time))
    )

    return(result)
}

#' Build the ggplot figure for displaying compartment sizes
#'
#' @param data A dataframe with the following columns: 
#'  * t: time
#'  * id: the group id (optional)
#'  * name: the name of the compartment. Options are S, E, I, and R
#'  * value: the size of the the compartment
#' @param title the title of the figure
#' @param compartment_labels the labels of the compartments (lexicographic order)
#' @param group_labels alternative labels for groups if there is more then one of them. If there is multiple and this is left as NULL then the values in the id column of `data` will be used.
#'
#' @return A ggplot figure
#'
#' @noRd
generate_figure <- function(data, title, compartment_labels, group_labels) {

    if ("id" %in% names(data)) {
        # if there still is an id column then we have multiple groups to plot

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
    figure <- figure + ggplot2::guides(linetype = ggplot2::guide_legend(title = "Group"))
    figure <- figure + ggplot2::scale_colour_discrete(
        name = "Compartment",
        labels = compartment_labels
    ) 
    figure <- figure + ggplot2::ggtitle(title)

    # add group labels if they've been supplied
    if (!is.null(group_labels)) {
        figure <- figure + ggplot2::scale_linetype_discrete(
            name = "Group", 
            labels = group_labels
        )
    }

    return(figure)
}


###############################################################################
# support functions

#' Get the level of a host's node in a transmission tree
#'
#' @param host A Host object
#' 
#' @return Integer
#'
#' @noRd 
node_level <- function(host) {

    # base case
    if (host$is_index) return(1L)

    # recursive
    return(node_level(host$infector) + 1L)
}

#' Get the shape of a host's node in a transmission tree
#' 
#' @description A index case will have a square node and all other hosts will have a circular node
#'
#' @param host A Host object
#' 
#' @return character
#'
#' @noRd 
node_shape <- function(host) {
    if (host$is_index) "square"
    else "circle"
}

#' Get the title of a host's node in a transmission tree
#'
#' @param host A Host object
#' 
#' @return character
#'
#' @noRd 
node_title <- function(host) {
    paste0(
        "<p>",
            "<b>ID:</b>         ", host$id, "<br>",
            "<b>Population:</b> ", host$population$id, "<br>",
            "<b>Exposure Time:</b> ", host$exposure_time, "<br>",
            "<b>Infectious Time:</b> ", host$infectious_time, "<br>",
            "<b>Recovery Time:</b> ", host$recovery_time, "<br>",
            "<b>Infector:</b> ", host$infector_id, "<br>",
            "<b>Infectees:</b> ", paste0(vapply(host$infectees, function(i) i$id, numeric(1L)), collapse = ", "), "<br>",
        "</p>"
    )
}

#' Get an interval for a host
#'
#' @param host A Host object (infected hosts only!)
#' @param interval The interval wanted: "serial", "generation", or "transmission" (default)
#' 
#' @return An integer
#'
#' @noRd 
interval <- function(host, interval) {
    if (host$is_index) NA
    else {
        if (interval == "serial") {
            # edge length is proportional to serial interval 
            host$infectious_time - host$infector$infectious_time
        } else if (interval == "generation") {
            # edge length is proportional to generation time 
            host$exposure_time - host$infector$exposure_time
        } else {
            # edge length is proportional to transmission time
            host$exposure_time - host$infector$infectious_time
        }
    }
}

#' Build a dataframe of nodes for a population
#'
#' @param population A Population object
#' 
#' @return A dataframe with the following columns
#'  * `group` The population that the host is in
#'  * `id` A node id with the form <host$population$id>-<host$id>
#'  * `level` The level of the node in the tree
#'  * `label` The label of the node. Host id is used here
#'  * `shape` The shape of the node
#'  * `title` The title displayed when a visNetwork node is hovered over
#'
#' @noRd 
population_nodes <- function(population) {
    
    infected_hosts <- population$hosts[vapply(
        population$hosts, 
        function(host) !host$is_susceptible(population$time),
        logical(1L)
    )]
    
    id <- vapply(infected_hosts, function(host) host$id, integer(1L))
    nodes <- data.frame(
        group = rep(population$id, length(infected_hosts)),
        id = paste0(population$id, "-", as.character(id)),
        level = vapply(infected_hosts, node_level, integer(1L)),
        label = as.character(id),
        shape = vapply(infected_hosts, node_shape, character(1L)),
        title = vapply(infected_hosts, node_title, character(1L))
    )
    return(nodes)
}

#' Build a dataframe of edges for a population
#'
#' @param population A Population object
#' 
#' @return A dataframe with the following columns
#'  * `from`  id of the node that the edge comes from
#'  * `to`  id of the node that the edge goes to
#'  * `serial_interval` The serial interval for this infector / infectee pair (NA if the infectee is an index case)
#'  * `generation_time` The generation time for this infector / infectee pair (NA if the infectee is an index case)
#'  * `transmission_interval` The transmission interval for this infector / infectee pair (NA if the infectee is an index case)
#' 
#' @noRd 
population_edges <- function(population) {

    infected_hosts <- population$hosts[vapply(
        population$hosts, 
        function(host) !host$is_susceptible(population$time),
        logical(1L)
    )]

    to_id <- vapply(infected_hosts, function(host) host$id, integer(1L))
    edges <- data.frame(
        from = vapply(
            infected_hosts, 
            function(host) {
                if (!host$is_index) {
                    paste0(host$infector$population$id, "-", host$infector$id)
                } else {
                    ""
                }
            },
            character(1L)
        ),
        to = paste0(population$id, "-", as.character(to_id)),
        serial_interval = vapply(infected_hosts, interval, numeric(1L), interval = "serial"),
        generation_time = vapply(infected_hosts, interval, numeric(1L), interval = "generation"),
        transmission_interval = vapply(infected_hosts, interval, numeric(1L), interval = "transmission")
    )

    return(edges)
}

#' Build a dataframe of nodes for an outbreak
#'
#' @param population An Outbreak object
#' 
#' @return A dataframe with the same columns as population_nodes
#' 
#' @noRd 
outbreak_nodes <- function(outbreak) {

    nodes <- do.call(rbind, lapply(outbreak$populations, population_nodes))

    return(nodes)
}

#' Build a dataframe of edges for an outbreak
#'
#' @param population An Outbreak object
#' 
#' @return A dataframe with the same columns as population_edges
#' 
#' @noRd 
outbreak_edges <- function(outbreak) {

    edges <- do.call(rbind, lapply(outbreak$populations, population_edges))

    return(edges)
}

###############################################################################
# exported functions
###############################################################################

#' Get a dataframe of all nodes in a outbreak transmission tree
#' 
#' @param x An Outbreak or Population object that has had it's outbreak simulated
#'
#' @return A dataframe
#' @export
#'
#' @examples 
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' population <- Population$new(
#'      1, 
#'      ref_strain, 
#'      inf_rate = 0.75, 
#'      find_infector_method = "serial", 
#'      si_shape = 6, si_rate = 1
#' )
#' population$run_simulation(Lab$new())
#' nodes <- transmission_tree_nodes(population)
#' 
transmission_tree_nodes <- function(x) {

    # validations for x
    stopifnot(R6::is.R6(x))
    if ("Population" %in% class(x)) {
        nodes <- population_nodes(x)
    } else if ("Outbreak" %in% class(x)) {
        nodes <- outbreak_nodes(x)
    } else {
        stop("`x` must be an Outbreak or Population object")
    }

    return(nodes)
}

#' Get a dataframe of all edges in a outbreak transmission tree
#' 
#' @param x An Outbreak or Population object that has had it's outbreak simulated
#'
#' @return A dataframe
#' @export
#'
#' @examples 
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' population <- Population$new(
#'      1, 
#'      ref_strain, 
#'      inf_rate = 0.75, 
#'      find_infector_method = "serial", 
#'      si_shape = 6, si_rate = 1
#' )
#' population$run_simulation(Lab$new())
#' edges <- transmission_tree_edges(population)
#' 
transmission_tree_edges <- function(x) {

    # validations for x
    stopifnot(R6::is.R6(x))
    if ("Population" %in% class(x)) {
        edges <- population_edges(x)
    } else if ("Outbreak" %in% class(x)) {
        edges <- outbreak_edges(x)
    } else {
        stop("`x` must be an Outbreak or Population object")
    }

    return(edges)
}

#' Plot a transmission tree for a disease outbreak
#' 
#' @param x An Outbreak or Population object that has had it's outbreak simulated
#' @param rooted Logical. Should the tree be rooted at the index cases 
#' @param edge_length Which interval the edge lengths should be proportional to
#'  * "serial": Use serial interval
#'  * "generation": Use the generation time
#'  * "transmission": Use the transmission interval (default)
#' 
#' @importFrom magrittr %>%
#'
#' @return A visNetwork plot
#' @export
#' 
#' @examples
#' # An SEIR model outbreak in a single homogenous population using the 
#' # serial interval method to find infectors
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' population <- Population$new(
#'      1, 
#'      ref_strain, 
#'      inf_rate = 0.75, 
#'      find_infector_method = "serial", 
#'      si_shape = 6, si_rate = 1
#' )
#' population$run_simulation(Lab$new())
#' figure <- plot_transmission_tree(population, rooted = FALSE, edge_length = "serial")
#' print(figure)
#' 
#' # An SEIR model outbreak with two populations using the tranmission method 
#' # for finding infectors.
#' # The two populations will have different transmission interval shapes
#' set.seed(1)
#' outbreak <- Outbreak$new(
#'     c(
#'         Population$new(1, ref_strain, init_inf = 5, init_sus = 95),
#'         Population$new(2, ref_strain, init_inf = 0, init_sus = 100)
#'     ),
#'     matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
#'     Lab$new()
#' )
#' outbreak$run_simulation(0)
#' 
#' figure <- plot_transmission_tree(outbreak)
#' print(figure)
#' 
plot_transmission_tree <- function(x, rooted = TRUE, edge_length = "transmission") {

    stopifnot(edge_length %in% c("serial", "generation", "transmission"))

    nodes <- transmission_tree_nodes(x)
    edges <- transmission_tree_edges(x)

    if (edge_length == "serial") {
        # edge length is proportional to serial edge_length 
        edges$length <- edges$serial_interval
    } else if (edge_length == "generation") {
        # edge length is proportional to generation time 
        edges$length <- edges$generation_time
    } else {
        # edge length is proportional to transmission time
        edges$length <- edges$transmission_interval
    }
    
    print(nodes$id)
    print(edges)
    graph <- visNetwork::visNetwork(nodes, edges, width = "100%") %>%
        visNetwork::visEdges(arrows = "to") %>%
        visNetwork::visInteraction(navigationButtons = TRUE) %>%
        visNetwork::visPhysics(stabilization = FALSE)

    if (rooted) {
        graph <- graph %>% visNetwork::visHierarchicalLayout()
    }

    return(graph)
}

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
            "<b>Group:</b> ", host$group$id, "<br>",
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

#' Build a dataframe of nodes for a group
#'
#' @param group A Group object
#' 
#' @return A dataframe with the following columns
#'  * `group` The group that the host is in
#'  * `id` A node id with the form <host$group$id>-<host$id>
#'  * `level` The level of the node in the tree
#'  * `label` The label of the node. Host id is used here
#'  * `shape` The shape of the node
#'  * `title` The title displayed when a visNetwork node is hovered over
#'
#' @noRd 
group_nodes <- function(group) {
    
    infected_hosts <- group$hosts[vapply(
        group$hosts, 
        function(host) !host$is_susceptible(group$time),
        logical(1L)
    )]
    
    id <- vapply(infected_hosts, function(host) host$id, integer(1L))
    nodes <- data.frame(
        group = rep(group$id, length(infected_hosts)),
        id = paste0(group$id, "-", as.character(id)),
        level = vapply(infected_hosts, node_level, integer(1L)),
        label = as.character(id),
        shape = vapply(infected_hosts, node_shape, character(1L)),
        title = vapply(infected_hosts, node_title, character(1L))
    )
    return(nodes)
}

#' Build a dataframe of edges for a group
#'
#' @param group A Group object
#' 
#' @return A dataframe with the following columns
#'  * `from`  id of the node that the edge comes from
#'  * `to`  id of the node that the edge goes to
#'  * `serial_interval` The serial interval for this infector / infectee pair (NA if the infectee is an index case)
#'  * `generation_time` The generation time for this infector / infectee pair (NA if the infectee is an index case)
#'  * `transmission_interval` The transmission interval for this infector / infectee pair (NA if the infectee is an index case)
#' 
#' @noRd 
group_edges <- function(group) {

    infected_hosts <- group$hosts[vapply(
        group$hosts, 
        function(host) !host$is_susceptible(group$time),
        logical(1L)
    )]

    to_id <- vapply(infected_hosts, function(host) host$id, integer(1L))
    edges <- data.frame(
        from = vapply(
            infected_hosts, 
            function(host) {
                if (!host$is_index) {
                    paste0(host$infector$group$id, "-", host$infector$id)
                } else {
                    ""
                }
            },
            character(1L)
        ),
        to = paste0(group$id, "-", as.character(to_id)),
        serial_interval = vapply(infected_hosts, interval, numeric(1L), interval = "serial"),
        generation_time = vapply(infected_hosts, interval, numeric(1L), interval = "generation"),
        transmission_interval = vapply(infected_hosts, interval, numeric(1L), interval = "transmission")
    )

    return(edges)
}

#' Build a dataframe of nodes for an population
#'
#' @param group A Population object
#' 
#' @return A dataframe with the same columns as `group_nodes()`
#' 
#' @noRd 
population_nodes <- function(population) {

    nodes <- do.call(rbind, lapply(population$groups, group_nodes))

    return(nodes)
}

#' Build a dataframe of edges for an population
#'
#' @param group A Population object
#' 
#' @return A dataframe with the same columns as `group_edges()`
#' 
#' @noRd 
population_edges <- function(population) {

    edges <- do.call(rbind, lapply(population$groups, group_edges))

    return(edges)
}

###############################################################################
# exported functions
###############################################################################

#' Get a dataframe of all nodes in an outbreak transmission tree
#' 
#' @param x A Population or Group object that has had it's outbreak simulated
#'
#' @return A dataframe
#' @export
#'
#' @examples 
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' group <- Group$new(
#'      1, 
#'      ref_strain, 
#'      inf_rate = 0.75, 
#'      find_infector_method = "serial", 
#'      si_shape = 6, si_rate = 1
#' )
#' group$run_simulation(Lab$new())
#' nodes <- transmission_tree_nodes(group)
#' 
transmission_tree_nodes <- function(x) {

    # validations for x
    stopifnot(R6::is.R6(x))
    if ("Group" %in% class(x)) {
        nodes <- group_nodes(x)
    } else if ("Population" %in% class(x)) {
        nodes <- population_nodes(x)
    } else {
        stop("`x` must be an Population or Group object")
    }

    return(nodes)
}

#' Get a dataframe of all edges in an outbreak transmission tree
#' 
#' @param x A Population or Group object that has had it's outbreak simulated
#'
#' @return A dataframe
#' @export
#'
#' @examples 
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' group <- Group$new(
#'      1, 
#'      ref_strain, 
#'      inf_rate = 0.75, 
#'      find_infector_method = "serial", 
#'      si_shape = 6, si_rate = 1
#' )
#' group$run_simulation(Lab$new())
#' edges <- transmission_tree_edges(group)
#' 
transmission_tree_edges <- function(x) {

    # validations for x
    stopifnot(R6::is.R6(x))
    if ("Group" %in% class(x)) {
        edges <- group_edges(x)
    } else if ("Population" %in% class(x)) {
        edges <- population_edges(x)
    } else {
        stop("`x` must be an Population or Group object")
    }

    return(edges)
}

#' Plot a transmission tree for a disease outbreak
#' 
#' @param x A Population or Group object that has had it's outbreak simulated
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
#' # An SEIR model outbreak in a single homogenous group using the 
#' # serial interval method to find infectors
#' set.seed(1)
#' ref_strain <- ReferenceStrain$new("ref_strain")
#' group <- Group$new(
#'      1, 
#'      ref_strain, 
#'      inf_rate = 0.75, 
#'      find_infector_method = "serial", 
#'      si_shape = 6, si_rate = 1
#' )
#' group$run_simulation(Lab$new())
#' figure <- plot_transmission_tree(group, rooted = FALSE, edge_length = "serial")
#' print(figure)
#' 
#' # An SEIR model outbreak in a population with two groups using the tranmission method 
#' # for finding infectors.
#' set.seed(1)
#' population <- Population$new(
#'     c(
#'         Group$new(1, ref_strain, init_inf = 5, init_sus = 95),
#'         Group$new(2, ref_strain, init_inf = 0, init_sus = 100)
#'     ),
#'     matrix(c(0.75, 0.5, 0.25, 0.1), ncol = 2),
#'     Lab$new()
#' )
#' population$run_simulation(0)
#' 
#' figure <- plot_transmission_tree(population)
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

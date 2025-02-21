% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-group.R
\name{Group}
\alias{Group}
\title{Group}
\description{
An R6 class representing a group of individuals. The Default model for an disease outbreak is SEIR but can be switched to an SIR model by setting \code{inc_shape} to zero.
}
\examples{
ref_strain <- ReferenceStrain$new(name = "ref_strain")
Group$new(id = 1, ref_strain = ref_strain)
Group$new(
 id = 1, 
 ref_strain = ref_strain, 
 find_infector_method = "serial",
 si_shape = 6,
 si_rate = 2 
)
Group$new(
 id = 1, 
 ref_strain = ref_strain, 
 find_infector_method = "random",
 max_init_dist = 20
)

## ------------------------------------------------
## Method `Group$get_infectees`
## ------------------------------------------------

ref_strain <- ReferenceStrain$new(name = "ref_strain")
group <- Group$new(id = 1, ref_strain = ref_strain)
group$get_infectees(3)


## ------------------------------------------------
## Method `Group$get_infectors`
## ------------------------------------------------

ref_strain <- ReferenceStrain$new(name = "ref_strain")
group <- Group$new(id = 1, ref_strain = ref_strain, init_inf = 10)
group$get_infectors(c(1, 2))
group <- Group$new(
 id = 1, ref_strain = ref_strain, 
 init_inf = 10, find_infector_method = "random"
)
group$get_infectors(c(0, 0, 0))

}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{id}}{The id of this group ** this must be unique if multiple groups are being used in an \code{Population} object}

\item{\code{init_sus}}{The number of initial susceptible individuals}

\item{\code{init_inf}}{The number of initial infected individuals}

\item{\code{inf_rate}}{The infectious pressure exerted by one infectious \code{host} on the rest of this \code{group}}

\item{\code{find_infector_method}}{The method with which an infectees infector is chosen (see definition of new for valid choices)}

\item{\code{si_shape, si_rate}}{The Gamma shape and rate parameters for the serial interval}

\item{\code{gt_shape, gt_rate}}{The Gamma shape and rate parameters for the generation time}

\item{\code{trans_int_shape, trans_int_rate}}{The Gamma shape and rate parameters for the transmission itnerval}

\item{\code{inc_shape, inc_rate}}{The Gamma shape and rate parameters for the incubation time}

\item{\code{rec_shape, rec_rate}}{The Gamma shape and rate parameters for the recovery time.}

\item{\code{ref_strain}}{The \code{ReferenceStrain} object used to create \code{Strain} objects for the index cases if there was any}

\item{\code{min_init_dist, max_init_dist}}{The minimum and maximum distances that index \code{Strain} objects should be from the \code{ref_strain}}

\item{\code{sample_schedule, sample_freq}}{The method used to schedule sampling times and the time interval between sampling (see new definition for valid choices):}

\item{\code{time}}{The current time for this group}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{si_shape, si_rate}}{The Gamma shape and rate parameters for the serial interval}

\item{\code{gt_shape, gt_rate}}{The Gamma shape and rate parameters for the generation time}

\item{\code{trans_int_shape, trans_int_rate}}{The Gamma shape and rate parameters for the transmission itnerval}

\item{\code{inc_shape, inc_rate}}{The Gamma shape and rate parameters for the incubation time}

\item{\code{rec_shape, rec_rate}}{The Gamma shape and rate parameters for the recovery time.}

\item{\code{min_init_dist, max_init_dist}}{The minimum and maximum distances that index \code{Strain} objects should be from the \code{ref_strain}}

\item{\code{sample_schedule, sample_freq}}{The method used to schedule sampling times and the time interval between sampling (see new definition for valid choices):}

\item{\code{exposed_size}}{The number of hosts in the exposed compartment at the groups current time}

\item{\code{hosts}}{All host objects in a vector}

\item{\code{hosts_due_for_sampling}}{A vector of hosts that are due for sampling at the group's current time}

\item{\code{is_outbreak_active}}{Logical indicating if there are still exposed or infectious individuals at the group's current time}

\item{\code{infectious_size}}{The number of hosts in the infectious compartment at the groups current time}

\item{\code{group_interval_stack}}{Pop a set of intervals from the interval stack. The interval stack is a preprepared set of at least incubation periods and the interval being used to find an infector (serial interval and generation time are also present if that method is being used)}

\item{\code{recovered_size}}{The number of hosts in the recovered compartment at the groups current time}

\item{\code{susceptible_size}}{The number of hosts in the susceptible compartment at the groups current time}

\item{\code{outbreak_size}}{The number of infected so far}

\item{\code{size}}{The group size}

\item{\code{time_since_infectious_hosts_infected}}{A vector containing the time elapsed wince all currently infectious hosts were infected}

\item{\code{time_since_infectious_hosts_infectious}}{A vector containing the time elapsed wince all currently infectious hosts became infectious}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Group-new}{\code{Group$new()}}
\item \href{#method-Group-exposed_hosts}{\code{Group$exposed_hosts()}}
\item \href{#method-Group-get_infectees}{\code{Group$get_infectees()}}
\item \href{#method-Group-get_infectors}{\code{Group$get_infectors()}}
\item \href{#method-Group-get_num_infectees}{\code{Group$get_num_infectees()}}
\item \href{#method-Group-infect}{\code{Group$infect()}}
\item \href{#method-Group-infectious_hosts}{\code{Group$infectious_hosts()}}
\item \href{#method-Group-print}{\code{Group$print()}}
\item \href{#method-Group-recovered_hosts}{\code{Group$recovered_hosts()}}
\item \href{#method-Group-run_simulation}{\code{Group$run_simulation()}}
\item \href{#method-Group-susceptible_hosts}{\code{Group$susceptible_hosts()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-new"></a>}}
\if{latex}{\out{\hypertarget{method-Group-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{group} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$new(
  id,
  ref_strain = NULL,
  init_sus = 99,
  init_inf = 1,
  inf_rate = 1,
  find_infector_method = "transmission",
  si_shape = NA,
  si_rate = NA,
  gt_shape = NA,
  gt_rate = NA,
  trans_int_shape = 3,
  trans_int_rate = 2,
  inc_shape = 3,
  inc_rate = 2,
  rec_shape = 5,
  rec_rate = 1,
  max_init_dist = 0,
  min_init_dist = 0,
  sample_schedule = "random",
  sample_freq = NA,
  host_class = Host,
  strain_class = Strain
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{The id of this group ** this must be unique if multiple groups are being used in an \code{Population} object}

\item{\code{ref_strain}}{The time that this strain is 'spawned'}

\item{\code{init_sus}}{The number of initial susceptible individuals (default 99)}

\item{\code{init_inf}}{The number of initial infected individuals (default 1)}

\item{\code{inf_rate}}{The infectious pressure exerted by one infectious \code{host} on the rest of this \code{group}}

\item{\code{find_infector_method}}{The method with which an infectees infector is chosen (see definition of new for valid choices)
\itemize{
\item 'transmission' (Default) Uses the transition interval to pick the infector. A transition interval is the time elapsed between an infector becoming infectious and the infectee being exposed
\item 'serial': Uses the serial interval to pick the infector. A serial interval is the time elapsed between an infector and their infectee becoming infectious
\item 'generation': Uses the generation time to pick the infector. A generation time is the time elapsed between an infector and their infectee being exposed
\item 'random': An infector is randomly selected from all currently infectious hosts
}}

\item{\code{si_shape, si_rate}}{The Gamma shape and rate parameters for the serial interval (Defaults are \code{6} and \code{2})}

\item{\code{gt_shape, gt_rate}}{The Gamma shape and rate parameters for the generation time}

\item{\code{trans_int_shape, trans_int_rate}}{The Gamma shape and rate parameters for the transmission interval}

\item{\code{inc_shape, inc_rate}}{The Gamma shape and rate parameters for the incubation time. If running an SIR simulation then set \code{inc_shape} to zero so that a host spends 0 time in the exposed compartment. (Defaults are \code{3} and \code{2} respectively)}

\item{\code{rec_shape, rec_rate}}{The Gamma shape and rate parameters for the recovery time. (Defaults are \code{5} and \code{1} respectively)}

\item{\code{min_init_dist, max_init_dist}}{The minimum and maximum distances that index \code{Strain} objects should be from the \code{ref_strain} (Default for both is \code{0})}

\item{\code{sample_schedule, sample_freq}}{The method used to schedule sampling times and the time interval between sampling
\itemize{
\item 'random': (Default) A random time is chosen between a \code{host}'s infectious time and their recovery time
\item 'calendar': Samples are taken from every infectious host at \code{sample_freq} intervals (starting from time = \code{sample_freq})
\item 'individual': A \code{host} is sampled every \code{sample_freq} time intervals up until they recover
}}

\item{\code{host_class}}{The class used to create hosts. Default is \code{Host} a subclassed class can be used here instead}

\item{\code{strain_class}}{The class used to create strains. Default is \code{Strain} a subclassed class can be used here instead}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new 'Group' object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-exposed_hosts"></a>}}
\if{latex}{\out{\hypertarget{method-Group-exposed_hosts}{}}}
\subsection{Method \code{exposed_hosts()}}{
Get all hosts that are in the exposed compartment at \code{time}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$exposed_hosts(time = self$time)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{time}}{The time being querried (Defaults to this group's time)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A vector of \code{Host} objects
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-get_infectees"></a>}}
\if{latex}{\out{\hypertarget{method-Group-get_infectees}{}}}
\subsection{Method \code{get_infectees()}}{
Select a \code{num_infectees} \code{host}s from the susceptible hosts in this \code{group}. All susceptible hosts are equally likely to be selected
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$get_infectees(num_infectees)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{num_infectees}}{An integer indicating the number of susceptible hosts to select}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A vector of \code{Host} objects
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{ref_strain <- ReferenceStrain$new(name = "ref_strain")
group <- Group$new(id = 1, ref_strain = ref_strain)
group$get_infectees(3)

}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-get_infectors"></a>}}
\if{latex}{\out{\hypertarget{method-Group-get_infectors}{}}}
\subsection{Method \code{get_infectors()}}{
Select infectors. \code{find_infector_method} determines how \code{infector_intervals} is used to find the infectors. The length of \code{infector_intervals} determines how many infectors are chosen.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$get_infectors(infector_intervals)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{infector_intervals}}{A vector of integers indicating how long ago an infector was infected / infectious (depends on \code{find_infector_method}). If \code{find_infector_method="random"} then pass a vector of zeros.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A vector of \code{Host} objects
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{ref_strain <- ReferenceStrain$new(name = "ref_strain")
group <- Group$new(id = 1, ref_strain = ref_strain, init_inf = 10)
group$get_infectors(c(1, 2))
group <- Group$new(
 id = 1, ref_strain = ref_strain, 
 init_inf = 10, find_infector_method = "random"
)
group$get_infectors(c(0, 0, 0))

}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-get_num_infectees"></a>}}
\if{latex}{\out{\hypertarget{method-Group-get_num_infectees}{}}}
\subsection{Method \code{get_num_infectees()}}{
Determine number of infectees
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$get_num_infectees(foi)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{foi}}{The force of infection exerted on this group}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An integer
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-infect"></a>}}
\if{latex}{\out{\hypertarget{method-Group-infect}{}}}
\subsection{Method \code{infect()}}{
Trigger an infection cycle in this \code{group}. All new infectees (if any) are selected, their infectors are selected and infections of hosts are carried out.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$infect()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
This group object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-infectious_hosts"></a>}}
\if{latex}{\out{\hypertarget{method-Group-infectious_hosts}{}}}
\subsection{Method \code{infectious_hosts()}}{
Get all hosts that are in the infectious compartment at \code{time}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$infectious_hosts(time = self$time)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{time}}{The time being querried (Defaults to this group's time)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A vector of \code{Host} objects
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-print"></a>}}
\if{latex}{\out{\hypertarget{method-Group-print}{}}}
\subsection{Method \code{print()}}{
Print a description of the \code{Group} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{arguments will be ignored}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
No return value. Description is printed to the console
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-recovered_hosts"></a>}}
\if{latex}{\out{\hypertarget{method-Group-recovered_hosts}{}}}
\subsection{Method \code{recovered_hosts()}}{
Get all hosts that are in the recovered compartment at \code{time}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$recovered_hosts(time = self$time)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{time}}{The time being querried (Defaults to this group's time)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A vector of \code{Host} objects
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-run_simulation"></a>}}
\if{latex}{\out{\hypertarget{method-Group-run_simulation}{}}}
\subsection{Method \code{run_simulation()}}{
Simulate an outbreak in this group. This is intended to be used if you are siming an outbreak in a single homogenous group. It's recommended that you use an \code{Population} class if your simulation requires a heterogenous population with more then one groups of hosts
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$run_simulation(lab, feedback = 10)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{lab}}{A lab object to collect samples from this group}

\item{\code{feedback}}{How frequently to print group summary to the console (turn off by setting to zero)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
This \code{group} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Group-susceptible_hosts"></a>}}
\if{latex}{\out{\hypertarget{method-Group-susceptible_hosts}{}}}
\subsection{Method \code{susceptible_hosts()}}{
Get all hosts that are in the susceptible compartment at \code{time}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Group$susceptible_hosts(time = self$time)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{time}}{The time being querried (Defaults to this group's time)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A vector of \code{Host} objects
}
}
}

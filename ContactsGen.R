# ContactsGen.R is part of Food INdustry CoViD Control Tool
# (FInd CoV Control), version 3.0.
# Copyright (C) 2020-2024 Cornell University.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


# Basic model: Employees are organized into crews, each with one foreman and
# one or more (ordinary) workers. Crews are organized into teams, each of which
# is under a supervisor. Supervisors, in turn, are all under a single manager.
# Crews may be of different sizes, and teams may include different numbers of crews.
# (In the current version of the complete model, all crews are the same size, and all
# teams include the same number of the crews, but this is not required by the
# functions in this file. If this requirement is not relaxed in future versions,
# this file may be simplifed by taking explicit account of it.


library(Matrix)

# The following are helper functions used in Primary_Matrices, not meant to be
# called directly by other functions

# Generates a single vector, for a single crew
Crew_Worker = function(size) {
    c(0, rep(1, size)) # 0 for the foreman, 1 for the workers
}

# Generates a single vector, for a single crew
Crew_Foreman = function(size) {
    c(1, rep(0, size)) # 1 for the foreman, 0 for the workers
}

# Prepends 0 to the first vector in a list, to represent the supervisor, for a
# list of crews, or the manager, for a list of teams. This is a little kludgy
# (the supervisor is not actually a member of the first crew on their team;
# the manager is not actually a member of the first team on the farm), but it
# generates the correct block matrix when fed to bdiag.
Prepend_Zero = function(list_) {
    list_[[1]] = c(0, list_[[1]])
    list_
}

# Generates a list of vectors for use in creating the crew_worker matrix.
# (And, after flattening into a single vector, for use in creating the
# team_worker matrix.)
# One call generates vectors for all crews in one team.
Team_Crew_Worker = function(sizes) {
    Prepend_Zero(lapply(sizes, Crew_Worker)) # Prepending 0 for the supervisor
}

# Generates a list of vectors for use in creating the crew_foreman matrix.
# (And, after flattening into a single vector, for use in creating the
# team_foreman matrix.)
# One call generates vectors for all crews in one team.
Team_Crew_Foreman = function(sizes) {
    Prepend_Zero(lapply(sizes, Crew_Foreman)) # Prepending 0 for the supervisor
}

# Generates a single vector, for a single team.
Team_Worker = function(sizes) {
    unlist(Team_Crew_Worker(sizes), use.names = FALSE)
}

# Generates a single vector, for a single team.
Team_Foreman = function(sizes) {
    unlist(Team_Crew_Foreman(sizes), use.names = FALSE)
}

# Generates a single vector, for a single team.
Team_Supervisor = function(sizes) {
    c(1, rep(0, sum(sizes) + length(sizes))) # 1 for the supervisor,
                                             # 0 for the foreman and workers.
}

#Prepends 0 to represent the manager
List_to_Matrix = function(list_) {
    as.matrix(bdiag(Prepend_Zero(list_)))
}


# Generates "primary" matrices for each combination of "level" ("crew,"
# "team," or "all") and employee type (worker, foreman, supervisor,
# manager).
#
# All have N rows, one for each individual in the population.
# The number of columns varies:
#   "crew" level matrices have one column per crew,
#   "team" level matrices have one column per team, and
#   "all" level matrices have only a single column.
#
# In all cases, a value of 1 indicates membership, and a value of 0 indicates
# lack of membership. Example: crew_worker[i,j] = 1 if and only if employee 1
# is an (ordinary) worker who is part of crew number j.
#
# For convenience, the order of employees reflects the structure of the
# workforce, i.e., employee #1 is the manager, each supervisor is followed by
# all the employees (foremen and workers) under that supervisor, and each
# foreman is followed by all the workers under that foreman.
#
# Parameters:
# crews_by_team: a vector of how many crews are in each team
# crew_sizes: a vector of how many workers (not including the foreman) are in
# each crew (ordered by team)
Primary_Matrices = function(crews_by_team, crew_sizes) {

    # Turns the single vector crew_sizes into a list
    # The latter could simply have been given as a parameter, rather than
    # using the two current parameters. This may or may not be changed
    # in the future.
    crew_sizes_split = split(crew_sizes, rep(1:length(crews_by_team),
                                             crews_by_team))

    crew_worker = List_to_Matrix(unlist(lapply(crew_sizes_split,
                                               Team_Crew_Worker),
                                        recursive = FALSE))
    crew_foreman = List_to_Matrix(unlist(lapply(crew_sizes_split,
                                                Team_Crew_Foreman),
                                         recursive = FALSE))

    team_worker = List_to_Matrix(lapply(crew_sizes_split, Team_Worker))
    team_foreman = List_to_Matrix(lapply(crew_sizes_split, Team_Foreman))
    team_supervisor = List_to_Matrix(lapply(crew_sizes_split, Team_Supervisor))

    all_worker = as.matrix(rowSums(team_worker))
    all_foreman = as.matrix(rowSums(team_foreman))
    all_supervisor = as.matrix(rowSums(team_supervisor))

    # 1 for the manager, 0 for workers, foremen, and supervisors
    all_manager = as.matrix(c(1, rep(0, sum(crew_sizes) + length(crew_sizes) +
                                     length(crews_by_team))))

    list(crew_worker = crew_worker,
         crew_foreman = crew_foreman,
         team_worker = team_worker,
         team_foreman = team_foreman,
         team_supervisor = team_supervisor,
         all_worker = all_worker,
         all_foreman = all_foreman,
         all_supervisor = all_supervisor,
         all_manager = all_manager)
}

# Generates "secondary" matrices for each combination of "level" ("same crew,"
# "same team [but different crew]" or "other") and two employee types (worker,
# foreman, supervisor, manager; can be the same or different).
#
# All have N rows and N columns, one for each individual in the population.
#
# For ease of future use, I've adopted a more compact notation than in
# Primary_Matrices: level is indicated by "c," "t," or "o," and employee types
# by "w," "f," "s," or "m," with an underscore separating level from types.
# Because these matrices are treated as symmetrical, only one is stored for,
# e.g., "t_wf" and "t_fw"; in all cases, the lower-level position is named first
# (e.g., "t_wf," not "t_fw").
#
# In all cases, a value of 1 indicates a pair of employees that have the
# types and share the level in question, and a value of 0 indicates a pair of
# employees that do not.

# Example: t_wf[i,j] = 1 if and only if one of employee i and employee j is a
# foreman and the other is a worker, and they are on the same team but not on
# the same crew.
Secondary_Matrices = function(primary_matrices) {
    cw = primary_matrices[['crew_worker']]
    tw = primary_matrices[['team_worker']]
    aw = primary_matrices[['all_worker']]

    cf = primary_matrices[['crew_foreman']]
    tf = primary_matrices[['team_foreman']]
    af = primary_matrices[['all_foreman']]

    ts = primary_matrices[['team_supervisor']]
    as_ = primary_matrices[['all_supervisor']] #underscore to get around builtin

    am = primary_matrices[['all_manager']]

    N = dim(cw)[1]

    c_ww = cw %*% t(cw)
    t_ww = tw %*% t(tw) - c_ww
    o_ww = aw %*% t(aw) - (c_ww + t_ww)
    diag(c_ww) = 0 # we do not count pairs of an individual with themself

    t_ff = tf %*% t(tf)
    o_ff = af %*% t(af) - t_ff
    diag(t_ff) = 0 # we do not count pairs of an individual with themself

    c_wf = cw %*% t(cf) + cf %*% t(cw) 
    t_wf = tw %*% t(tf) + tf %*% t(tw) - c_wf
    o_wf = aw %*% t(af) + af %*% t(aw) - (c_wf + t_wf)
    #diag inherently 0

    o_ss = as_ %*% t(as_)
    diag(o_ss) = 0 # we do not count pairs of an individual with themself

    t_ws = tw %*% t(ts) + ts %*% t(tw)
    o_ws = aw %*% t(as_) + as_ %*% t(aw) - t_ws
    #diag inherently 0

    t_fs = tf %*% t(ts) + ts %*% t(tf)
    o_fs = af %*% t(as_) + as_ %*% t(af) - t_fs
    #diag inherently 0

    #o_mm is inherently 0, given that we do not count pairs of an individual
    #with themself, as there is only one manager.

    o_wm = aw %*% t(am) + am %*% t(aw)
    o_fm = af %*% t(am) + am %*% t(af)
    o_sm = as_ %*% t(am) + am %*% t(as_)
    #diag inherently 0 for all three

    list(c_ww = c_ww,
         t_ww = t_ww,
         o_ww = o_ww,
         t_ff = t_ff,
         o_ff = o_ff,
         c_wf = c_wf,
         t_wf = t_wf,
         o_wf = o_wf,
         o_ss = o_ss,
         t_ws = t_ws,
         o_ws = o_ws,
         t_fs = t_fs,
         o_fs = o_fs,
         o_wm = o_wm,
         o_fm = o_fm,
         o_sm = o_sm
    )
}

# Calculate work contact rates (not binary yes/no) from a list of secondary
# matrices, a list of rates for each type of pair, and an average (daily)
# contact rate _per person_ (not per pair).
# NB: As currently coded, this will give wrong results if the order of the
# matrices and the order of the rates are not the same. Care is required.
Work_Contacts = function(secondary_matrices, rates, average) {
    raw = Reduce('+', mapply('*', secondary_matrices, rates, SIMPLIFY = FALSE))
    
    raw * average / (mean(raw) * dim(raw)[1]) #note: could equally well use
                                              #dim(raw)[2]; since the matrices
                                              #are square, these are identical
}


#Since home contacts are currently being modeled as homogeneous mixing, not
#writing code to generate matrices for that.

#sample_rates here is a bunch of provisional assumptions, since we need
#*something* to use here. A future version of the model may change these
#values; it may also allow end users to specify their own.
example_rates = list(c_ww = 1,      #by definition (reference value);
                                    #these are relative rates,
                                    #so making the rate at which a worker
                                    #contacts another worker on the same crew 1
                                    #may be a natural way to scale it
                    t_ww = 0.1,
                    o_ww = 0.1,
                    t_ff = 0.2,
                    o_ff = 0.1,
                    c_wf = 1,
                    t_wf = 0.1,
                    o_wf = 0.1,
                    o_ss = 0.2,
                    t_ws = 0.3,
                    o_ws = 0.1,
                    t_fs = 1,
                    o_fs = 0.1,
                    o_wm = 0.01,
                    o_fm = 0.1,
                    o_sm = 1
    )

#Okay, now for a single master function:
ContactsGen <- function(crews_by_team, crew_sizes, contact_rates, average) {
    primary_matrices = Primary_Matrices(crews_by_team, crew_sizes)
    secondary_matrices = Secondary_Matrices(primary_matrices)
    work_contacts = Work_Contacts(secondary_matrices, contact_rates, average)

    work_contacts
}

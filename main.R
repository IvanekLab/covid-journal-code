# main.R is part of Food INdustry CoViD Control Tool
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


#ABM model
#individual based
#each person/agent has properties = parameters, attributes behaviours
#model each person individually
main_produce_farm_fn = function(kConstants) { #goal: get more meaningful debug data
#library(Rlab) #TBD renable this somewhere sensible
 
source("AgentGen.R")
source("ContactsGen.R")
source("ABM.R")

#General note: foo = get('bar', baz) is similar to foo = baz[['bar']], *except*
#that it will throw an error if baz has no element named 'bar', instead of
#setting foo = NULL
#This helps to make errors due to typos, or to incomplete updating following
#code refactoring, happen where the mistake was made, instead of propagating
#silently. Consequently, I'm working on shifting to this idiom (even if it is a
#bit less "R-like") whereever there isn't a reason not to. "Fail early, fail
#loudly!"


#These shouldn't vary on a scenario-by-scenario or facility-by-facility basis,
#with the possible exception of "what if a more human-to-human transmissible
#strain evolves and becomes dominant"--type scenarios.
    #2021-07-18 Oh hey, foreshadowing!
#Otherwise, modification should occur only in response to updates in knowledge
#about viral properties.
#Default values from Moghadas et al. 2020 (gives range of .0575 to .0698 for
#p_trans_IP)
p_trans_IP = get('p_trans_IP', kConstants)
p_trans_IA = get('p_trans_IA', kConstants)
p_trans_IM = get('p_trans_IM', kConstants)

virus_parameters = list(p_trans_IP = p_trans_IP,
                        p_trans_IA = p_trans_IA,
                        p_trans_IM = p_trans_IM
                        )

#some of these are redundant steps that should be eliminated
ScenarioParameters = function(work_R0, dormitory_R0, days, housing_dormitory,
                              work_testing_rate,
                              home_vaccination_rate, lambda = 0, crews_by_team,
                              crew_sizes, virus_params) {
    p_trans_IP = get('p_trans_IP', virus_params)
    p_trans_IA = get('p_trans_IA', virus_params)
    p_trans_IM = get('p_trans_IM', virus_params)

    #marginal values based on default age distribution
    #should probably be made more flexible later, but this should do for now

    p_symptomatic = 0.671
    duration_IA = get('duration_IA_mean', kConstants)
    duration_IP = get('duration_IP_mean', kConstants)
    duration_IM = get('duration_IM_mean', kConstants)
    #p_trans_IS = .89 * .0575, # from Moghadas et al. 2020, but irrelevant
                               # (for now, at least) since we are assuming all
                               # severe cases are hospitalized

    #R0 per contact per day, = 0.2349583
    r0pcpd = p_symptomatic * (duration_IP * p_trans_IP +
                              duration_IM * p_trans_IM) +
            (1 - p_symptomatic) * duration_IA * p_trans_IA

    work_contacts = work_R0 / r0pcpd
    dormitory_contacts = dormitory_R0 / r0pcpd
    par1 = list(average = work_contacts,
                lambda = lambda, #force of infection "from the community" per
                                 #home shift, if in community housing
                nTime1 = days,  #duration of one simulation iteration, in days
                crews_by_team = crews_by_team,
                crew_sizes = crew_sizes,
                dormitory_intensity = dormitory_contacts,
                housing_dormitory = housing_dormitory,
                home_vaccination_rate = home_vaccination_rate,
                work_testing_rate = work_testing_rate
    )
    par1
}

##########################################
###     set scenario and output files 
############################################
scenario_parameters = ScenarioParameters(work_R0 = net_work_R0,
                                         dormitory_R0 = dormitory_R0,
                                         days = days,
                                         housing_dormitory = TRUE, 
                                         work_testing_rate = work_testing_rate,
                                         home_vaccination_rate = vaccination_rate,
                                         lambda = community_foi,
                                         crews_by_team = crews_by_team, 
                                         crew_sizes = crew_sizes,
                                         virus_params = virus_parameters)

###################################
##################################

#create a nPop1 x nPop1 matrix for contacts in the population and extract list
#of individuals met by each Agent
####
#edits for facility here
####
lambda_home = scenario_parameters$lambda

if(farm_or_facility == 'farm') {
    work_contacts <- ContactsGen(scenario_parameters$crews_by_team,
                                 scenario_parameters$crew_sizes,
                                 example_rates,
                                 scenario_parameters$average)

    production_shift_1 = work_contacts
    production_shift_2 = matrix(0, N, N)
    cleaning_shift_full = matrix(0, N, N)
    shift_sum = work_contacts

    on_ps_1 = rep(1, N)
    on_ps_2 = rep(0, N)
    on_cs = rep(0, N)

} else { #only alternative that we allow is facility

    source('custom-contacts-gen-general.R')

    contacts_matrices = facility_contacts_gen(
        workers_per_line = workers_per_crew,
        #keeping old names for analogous parameters in full_run
        n_lines = crews_per_supervisor,
        n_production_shifts = supervisors,
        n_shift_floaters = n_shift_floaters,
        n_cleaners = n_cleaners,
        n_all_floaters = n_all_floaters
    )

    production_shift_1 = contacts_matrices[['production_shift_1']]
    production_shift_2 = contacts_matrices[['production_shift_2']]
    cleaning_shift_full = contacts_matrices[['cleaning_shift_full']]
    shift_sum =  contacts_matrices[['shift_sum']]


    ####
    #TBD (eventually): Move this to the contacts generation file
    ####
    psX_only_size = 1 + workers_per_crew * crews_per_supervisor + n_shift_floaters
    if(supervisors > 1) {
        on_ps_1 = c(1/3, rep(1, psX_only_size), rep(0, psX_only_size),
                    rep(0, n_cleaners), rep(1/3, n_all_floaters))
        on_ps_2 = c(1/3, rep(0, psX_only_size), rep(1, psX_only_size),
                    rep(0, n_cleaners), rep(1/3, n_all_floaters))
        on_cs = c(1/3, rep(0, 2 * psX_only_size), rep(1, n_cleaners),
                rep(1/3, n_all_floaters))
    } else {
        on_ps_1 = c(1/2, rep(1, psX_only_size), rep(0, n_cleaners),
                    rep(1/2, n_all_floaters))
        on_ps_2 = rep(0, 1 + psX_only_size + n_cleaners + n_all_floaters)
        on_cs = c(1/2, rep(0, psX_only_size), rep(1, n_cleaners),
                rep(1/2, n_all_floaters))
    }
}



if(any(on_ps_1 + on_ps_2 + on_cs != rep(1, N))) {
    stop('Some presences do not add up to 1.')
}

###
#working on proper dormitory_contacts parameters
#treating dormitory as the shift after work (or work shift and shift after work
#on the weekend)
#Actually, this really belongs in ContactsGen or its replacement, because it's
#actually non-trivial.
#TBD (eventually): Move this crap to ContactsGen and its facility analogue.
###

quantitative_presence_list = list(ps_1 = on_ps_1,
                           ps_2 = on_ps_2,
                           cs =   on_cs,
                           weekend_ps_1 = 0,
                           weekend_ps_2 = 0,
                           weekend_cs = 0)

lambda_list = list(ps_1 = on_cs * lambda_home,
                   ps_2 = on_ps_1 * lambda_home,
                   cs = on_ps_2 * lambda_home,
                   weekend_ps_1 = (on_cs + on_ps_1) * lambda_home,
                   weekend_ps_2 = (on_ps_1 + on_ps_2) * lambda_home,
                   weekend_cs = (on_ps_2 + on_cs) * lambda_home)

make_contact_matrix = function(v) {
    v = matrix(v)
    M = v %*% t(v)
    diag(M) = 0
    M
}

raw_home_contacts_ps_1 = make_contact_matrix(on_cs)
raw_home_contacts_ps_2 = make_contact_matrix(on_ps_1)
raw_home_contacts_cs = make_contact_matrix(on_ps_2)
raw_home_contacts_weekend_ps_1 = make_contact_matrix(on_cs + on_ps_1)
raw_home_contacts_weekend_ps_2 = make_contact_matrix(on_ps_1 + on_ps_2)
raw_home_contacts_weekend_cs = make_contact_matrix(on_ps_2 + on_cs)

average_raw_home_contacts_per_day = (5 * (raw_home_contacts_ps_1 +
                                          raw_home_contacts_ps_2 +
                                          raw_home_contacts_cs) +
                                     2 * (raw_home_contacts_weekend_ps_1 +
                                          raw_home_contacts_weekend_ps_2 +
                                          raw_home_contacts_weekend_cs)
                                     ) / 7

#a quick sanity check suggested the dominant eigenvalue is actually pretty close
#to the mean of the rowSums, so we'll use the latter for now, for consistency
#with work_R0

home_scaling_factor = ifelse(scenario_parameters$dormitory_intensity == 0,
    0,
    scenario_parameters$dormitory_intensity /
            mean(rowSums(average_raw_home_contacts_per_day))
)
#no need for 7/9, since this is directly calculated from the weekly sum
#but ifelse *is* needed to avoid a 0/0 issue if housing is in the community
work_scaling_factor = ifelse(scenario_parameters$average == 0,
    0,
    scenario_parameters[['average']] / (sum(shift_sum) / N) * 7/5
)    #this could (and perhaps should) be done with a weekly sum as well

contacts_list = list(ps_1 = production_shift_1 * work_scaling_factor +
                            raw_home_contacts_ps_1 * home_scaling_factor,
                     ps_2 = production_shift_2 * work_scaling_factor +
                            raw_home_contacts_ps_2 * home_scaling_factor,
                     cs = cleaning_shift_full * work_scaling_factor +
                          raw_home_contacts_cs * home_scaling_factor,
                     weekend_ps_1 = raw_home_contacts_weekend_ps_1 *
                            home_scaling_factor,
                     weekend_ps_2 = raw_home_contacts_weekend_ps_2 *
                         home_scaling_factor,
                     weekend_cs = raw_home_contacts_weekend_cs *
                         home_scaling_factor)


vaccination_rate_list = list(
        ps_1 = on_cs * scenario_parameters$home_vaccination_rate,
        ps_2 = on_ps_1 * scenario_parameters$home_vaccination_rate,
        cs   = on_ps_2 * scenario_parameters$home_vaccination_rate
)
boosting_rate_list = list(
        ps_1 = on_cs * boosting_rate,
        ps_2 = on_ps_1 * boosting_rate,
        cs   = on_ps_2 * boosting_rate
)
vaccination_rate_list[['weekend_ps_1']] = vaccination_rate_list[['ps_1']]
vaccination_rate_list[['weekend_ps_2']] = vaccination_rate_list[['ps_2']]
vaccination_rate_list[['weekend_cs']] = vaccination_rate_list[['cs']]

boosting_rate_list[['weekend_ps_1']] = boosting_rate_list[['ps_1']]
boosting_rate_list[['weekend_ps_2']] = boosting_rate_list[['ps_2']]
boosting_rate_list[['weekend_cs']] = boosting_rate_list[['cs']]

vaccination_interval = get('second_shot_interval', kConstants)

workday = c('ps_1', 'ps_2', 'cs')
day_off = c('weekend_ps_1', 'weekend_ps_2', 'weekend_cs')
step_length_list = list(ps_1 = 1/3, ps_2 = 1/3, cs = 1/3, weekend_ps_1 = 1/3,
                        weekend_ps_2 = 1/3, weekend_cs = 1/3)
testing_rate_list = list(ps_1 = get('work_testing_rate', scenario_parameters),
                         ps_2 = get('work_testing_rate', scenario_parameters),
                         cs =   get('work_testing_rate', scenario_parameters),
                         weekend_ps_1 = 0,
                         weekend_ps_2 = 0,
                         weekend_cs = 0)
steps = scenario_parameters$nTime1 * 3
step_index = (1:steps) * (1/3) #step_length

###### code to run simulation with num_sims iterations
source('safe-random-functions.R')
if(!exists('FIXED_SEED') || FIXED_SEED == TRUE) {
    safe_set_seed(-778276078) #random 32-bit signed integer generated using
                              #atmospheric noise for reproducible output
}

sys_time_start = Sys.time()
for (i in 1:num_sims) {
    start_day = sample(1:7, 1)
    #week = c(rep(workday, 5), rep(day_off, 2))
    if(start_day %in% 1:5) {
        week = c(rep(workday, 6 - start_day),
                 rep(day_off, 2),
                 rep(workday, start_day - 1))
    } else {
        week = c(rep(day_off, 8 - start_day),
                 rep(workday, 5),
                 rep(day_off, start_day - 6))
    }
    schedule = rep(week, ceiling(days/7))[1:(3 * days)]

    floater_randomizers = runif(N, 0, 1)
    #print(paste0('Test floater_randomizers: ', floater_randomizers[1], '\t', floater_randomizers[N]))
    if(i %in% 1:2) {
        cat('intervention:', index_i, 'run i:', i, 'After floater_randomizers; Test value:', runif(1, 0, 1), '\n')
    }
    on_ps_1_randomized = floater_randomizers <= on_ps_1
    on_ps_2_randomized = on_ps_1 < floater_randomizers &
                         floater_randomizers <= on_ps_1 + on_ps_2
    on_cs_randomized = on_ps_1 + on_ps_2 < floater_randomizers 
    #on_ps_1 + on_ps_2 + on_cs was already checked to == rep(1, N)
    if(any(on_ps_1_randomized + on_ps_2_randomized + on_cs_randomized != 1)) {
        stop('Shift randomization for vaccination and testing somehow failed.')
    }

    agent_presence_list = list(ps_1 = on_ps_1_randomized,
                           ps_2 = on_ps_2_randomized,
                           cs = on_cs_randomized,
                           weekend_ps_1 = FALSE,
                           weekend_ps_2 = FALSE,
                           weekend_cs = FALSE)

    agents <- AgentGen(N, E0 = n_exposed, IA0 = 0, IP0 = 0, IM0 = n_mild,
                       initial_recovered = initial_recovered,
                       initial_V1 = initial_V1, initial_V2 = initial_V2,
                       ffv_last_five_months = ffv_last_five_months,
                       fraction_boosted_ever = fraction_boosted_ever,
                       fraction_boosted_last_five_months =
                           fraction_boosted_last_five_months,
                       protection_functions = protection_functions,
                       kConstants = kConstants)
    if(i %in% 1:2) {
        cat('intervention:', index_i, 'run i:', i, 'After AgentGen; Test value:', runif(1, 0, 1), '\n')
    }
    model <- ABM(agents, contacts_list = contacts_list,
                 lambda_list = lambda_list, schedule = schedule,
                 virus_parameters, testing_parameters,
                 vaccination_interval,
                 scenario_parameters,
                 steps = steps, step_length_list = step_length_list,
                 testing_rate_list = testing_rate_list,
                 vaccination_rate_list = vaccination_rate_list,
                 agent_presence_list = agent_presence_list,
                 quantitative_presence_list = quantitative_presence_list,
                 boosting_rate_list = boosting_rate_list,
                 protection_functions = protection_functions,
                 kConstants
    )
    agents = model$agents
    output = model$Out1

    if(i == 1) {
        full_output = array(0, c(steps, dim(as.matrix(output))[2], num_sims))
        #doing this this way guarantees it gets created with the right number of
        #return variables
        start_days = numeric(num_sims)
    }
    full_output[,,i] = as.matrix(output) #this works; for whatever reason,
                                         #as.array does not
    start_days[i] = start_day
} # for (i in 1:num_sims)
#print_rand_state(paste('intervention:', index_i, 'printing state'))

cat('intervention:', index_i, 'All runs completed; Test value:', runif(1, 0, 1), '\n')

sys_time_end = Sys.time()
cat(sys_time_end - sys_time_start, 'for', row_name,'\n')

colnames(full_output) = colnames(output)
#cat('\n\nfull_output_save_name:', full_output_save_name, '\n\n')
saveRDS(full_output, full_output_save_name)
#cat('through the first hurdle\n')
fragments = unlist(strsplit(full_output_save_name, '/'))
#cat(fragments[1],'\n',
        #start_days = readRDS(paste0(fragments[1], '/start_days--', fragments[2]))
start_days_save_name = paste0(fragments[1], '/start_days--', fragments[2])
#cat('\n\nstart_days_save_name:', start_days_save_name, '\n\n')
saveRDS(start_days, start_days_save_name)
} #main_produce_farm_fn 

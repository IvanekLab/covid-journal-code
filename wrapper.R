# wrapper.R is part of Food INdustry CoViD Control Tool
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


# behavior of wrapper.R is controlled by other files using several tests for
# all-caps variables. Production code has not run wrapper.R directly for some,
# time now; if it is still possible to run it directly without further
# modification, then the following should evaluate to FALSE:
# exists('ANALYZE') && ANALYZE == TRUE
# exists('DOUBLE_WRAPPED') && DOUBLE_WRAPPED == TRUE

# If using double-wrapped.R to run a batch of scenarios, a number of
# parameters may be overwritten by the "double_wrapped" variables; as a reminder
# of this, they are marked with comments of the form:
# batch mode: [double_wrapped variable name]
# These comments can be ignored entirely if not intending to ever use batch mode
# (assuming that not using batch mode still works).

# Number of individuals initially in various states other than completely
# susceptible
wrapper_fn = function(index_i, kConstants) {

initial_recovered = 0 # batch mode: double_wrap_initial_recovered
initial_V1 = 0        # batch mode: double_wrap_initial_V1
initial_V2 = 0        # batch mode: double_wrap_initial_V2

###########################################
# Flags controlling various interventions #
###########################################

# Social distancing interventions
theoretical_social_distancing_R0_reduction = 0

#surveillance -- currently, these cannot both be true,
#but that limitation may at some point be lifted
temperature_screening = FALSE
temperature_threshold = 38      # batch mode: double_wrap_temp_test
                                # values currently recognized: 38
                                # combining other values with
                                # temperature_testing = TRUE will cause an error
viral_testing_rate = 0.3        # batch mode: double_wrap_viral_test_rate
rational_testing = TRUE

#vaccination
vaccination_rate = 0            # batch mode: double_wrap_vax_rate
rational_vaccination = FALSE    # TRUE is not currently implemented, but may be
                                # at some point in the future


##############################
# Program control parameters #
##############################
num_sims = 100 # number of simulations; 100 takes 1-2 minutes
               # (less on a faster machine); for production use,
               # most likely want 1000+
               # batch mode: double_wrap_num_sims


################################################################
# Code to allow batch running, by overriding parameters        #
# set in this file with those set in the invoking file.        #
# If not using batch mode, this section may safely be ignored. #
################################################################

if(exists('DOUBLE_WRAPPED') && DOUBLE_WRAPPED == TRUE) {
    #TBD (eventually): Really, this whole setup should be handled by parameters
    #passed into the function either from double-wrapped.R or elsewhere
    initial_recovered = double_wrap_initial_recovered
    theoretical_social_distancing_R0_reduction = double_wrap_reduction
    temperature_screening = (double_wrap_temp_test != FALSE)
    temperature_threshold = double_wrap_temp_test
    viral_testing_rate = double_wrap_viral_test_rate
    vaccination_rate = double_wrap_vax_rate
    community_foi = double_wrap_community_foi
    num_sims = double_wrap_num_sims
    baseline_work_R0 = double_wrap_baseline_work_R0
    initial_V1 = double_wrap_initial_V1
    initial_V2 = double_wrap_initial_V2
    rational_testing = double_wrap_rational_testing
}

#############################################################################
# Derived values; simple use will not require modification below this point #
#############################################################################

social_distancing_R0_factor = (
                               (1 - theoretical_social_distancing_R0_reduction)
                              )

total_R0_factor = social_distancing_R0_factor

net_work_R0 = total_R0_factor * baseline_work_R0
#for now, assuming that R0 reductions only apply at work

if(temperature_screening) {
    if(viral_testing_rate != 0) {
        stop() #both at once is not yet implemented
    }

    if (temperature_threshold == 38) {
        sensitivity = .052
        specificity = 1.00
    } else {
        stop('Invalid screening temperature')
    }
    
    testing_parameters = list(
        ### the following are parameters for virus testing
        asymptomatic_FNR = specificity, #presumably, we'll catch a few by chance
        presymptomatic_FNR = specificity,
        mild_FNR = 1 - sensitivity,
        FPR = 1 - specificity,
        rational_testing = rational_testing
    )

    work_testing_rate = 1
} else {
    testing_parameters = list(
        ### the following are parameters for virus testing
        asymptomatic_FNR = 1 - 0.9,
        presymptomatic_FNR = 1 - 0.9,
        mild_FNR = 1 - 0.9,
        FPR = 1 - 0.9995,
        rational_testing = rational_testing
    )
    work_testing_rate = viral_testing_rate
}

#########################
# filename construction #
#########################

filename_core = paste(subdirectory, unique_id, '_community-', community_foi,
                      ',work_R0-', baseline_work_R0, sep = '')

R0_factor_string = paste()
if(total_R0_factor != 1) {
   filename_core = paste(filename_core , 'x(1-', 1 - total_R0_factor, ')',
                         sep = '')
}

if(dormitory_R0 > 0) {
    filename_core = paste(filename_core, ',dormitory_R0-', dormitory_R0,
                          sep = '')
}

filename_core = paste(filename_core, ',E0-', n_exposed, sep = '')

if(n_mild > 0) {
    filename_core = paste(filename_core, ',IM0-', n_mild, sep = '')
}

if(temperature_screening) {
    filename_core = paste(filename_core, ',T.test-', temperature_threshold,
                          sep = '')
}

if(viral_testing_rate != 0) {
    filename_core = paste(filename_core, ',v.test-', viral_testing_rate,
                          sep = '')
    if(rational_testing) {
        filename_core = paste(filename_core, '-rational', sep = '')
    }
}

if(vaccination_rate != 0) {
    filename_core = paste(filename_core, ',vax-rate', vaccination_rate,
                          sep = '')
}


if(initial_recovered > 0) {
    filename_core = paste(filename_core, ',initial_recovered-',
                          initial_recovered, sep = '')
}

if(initial_V2 > 0) {
    filename_core = paste(filename_core, ',initial_V2-', initial_V2, sep = '')
}

if(initial_V1 > 0) {
    filename_core = paste(filename_core, ',initial_V1-', initial_V1, sep = '')
}

filename_core = paste(filename_core, ',n_sims-', num_sims, 'index_i-', index_i,
                      sep = '')
full_output_save_name = paste(filename_core, '_full-output.rds', sep = '')

#and now let's run the model
if(!(exists('ANALYZE') && ANALYZE == TRUE)) {
    source('main.R', local = TRUE)
    main_produce_farm_fn(kConstants)
}

full_output_save_name
} #wrapper_fn

# AgentGen.R is part of Food INdustry CoViD Control Tool
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

# create database of all agents and their attributes
#The default age_probabilities vector is based on finding the MLE 4-parameter
#beta distribution that best fits age category data on agricultural workers


source('update-immunity-function.R')

AgentGen <- function (N, E0 = 1, IA0 = 0, IP0 = 0, IM0 = 0,
                      initial_recovered = 0, initial_V1 = 0, initial_V2 = 0,
                      ffv_last_five_months, #TBD (eventually) make fraction vs.
                                            #number consistent across parameters
                                            #TBD (eventually): rename to be more
                                            #general?
                      age_probabilities = c(0.04, 0.26, 0.26, 0.21, 0.15, 0.07,
                                            0.01, 0),
                      fraction_boosted_ever = 0,
                      fraction_boosted_last_five_months = 0,
                      protection_functions,
                      kConstants
                      ) {
    #TBD (eventually): Either add back in the ability to use these or remove
    #them from the parameter list.
    if(max(IA0, IP0, initial_V1) > 0) {
        stop('Attempt to use buggy functionality in AgentGen.')
    }

    net_protection = get('net_protection', protection_functions)
    infection_protection = get('infection_protection', protection_functions)
    boosting_interval = get('boosting_interval', kConstants)
    second_shot_interval = get('second_shot_interval', kConstants)
    SEVERE_MULTIPLIER = get('SEVERE_MULTIPLIER', kConstants)
    R_question_period = get('R_question_period', kConstants)
    time_since_first_V2 = get('time_since_first_V2', kConstants)

    Age_Categories = c("10-19", "20-29", "30-39", "40-49", "50-59", "60-69",
                       "70-79", "80+")

    ##duration_E calculations
    #Retrieving and calculating constants
    mu = get('mu', kConstants)
    sd = get('sd', kConstants)

    duration_IA_mean = get('duration_IA_mean', kConstants)
    duration_IA_shape = get('duration_IA_shape', kConstants)
    duration_IA_scale = duration_IA_mean / duration_IA_shape

    duration_IP_mean = get('duration_IP_mean', kConstants)
    duration_IP_shape = get('duration_IP_shape', kConstants)
    duration_IP_scale = duration_IP_mean / duration_IP_shape

    duration_IM_mean = get('duration_IM_mean', kConstants)
    duration_IM_shape = get('duration_IM_shape', kConstants)
    duration_IM_scale = duration_IM_mean / duration_IM_shape

    duration_IS_mean = get('duration_IS_mean', kConstants)
    duration_IS_shape = get('duration_IS_shape', kConstants)
    duration_IS_scale = duration_IS_mean / duration_IS_shape

    duration_IC_mean = get('duration_IC_mean', kConstants)
    duration_IC_shape = get('duration_IC_shape', kConstants)
    duration_IC_scale = duration_IC_mean / duration_IC_shape
    #END RETRIEVING AND CALCULATING CONSTANTS

    sdlog = sqrt(log(sd**2 / mu**2 + 1))
    mulog = log(mu) - (sdlog ** 2) / 2


    duration_IP = sgamma(N, shape=duration_IP_shape, scale=duration_IP_scale)
    duration_E = pmax(slnorm(N, mulog, sdlog) - duration_IP, 0) 
    #Moghadas et al., 2020 & need for a non-negative duration
    #Create a population of susceptibles 
    agents <- data.frame(ID = 1:N,          #Unique ID number,
                                            #Allows contact tracing
                                            #(in a future version)

                         #All Initially susceptible
                         #Will be overwritten later in this
                         #function for some employees.
                         infection_status = 'NI', #not infected
                         immune_status = 'FS', #fully susceptible
                         vax_status = 'NV', #not vaccinated
                         infector_ID = 0,    
                         infector_state = '', #state of source of infection
                                              #(if any)
                         #The following are conceptualized as time *of*
                         #transition into the state in question, not time
                         #*since* transition. Thus, the value for a state that
                         #has not (yet) been reached (and hence, a transition
                         #that has not happened (yet) is infinity, not 0.
                         #This may in the future be changed to NA.
                         time_E = Inf,  
                         time_IA = Inf,  
                         time_IP = Inf,  
                         time_IM = Inf,
                         time_IS = Inf,
                         time_IC = Inf,
                         time_R = Inf,
                         time_V1 = Inf,
                         time_V2 = Inf,
                         time_B = Inf,
                         boosting_on_time = NA, # will be replaced later for
                                                # everyone
                         time_isolated = Inf,  #setup for time in Isolation
                                               #unlike some of the states listed
                                               #above, this is not a mutually
                                               #exclusive state, e.g., someone
                                               #can be both IM and isolated.
                         time_tested = -Inf,
                         #Unlike most times, we want the "last time at which
                         #this person was tested," for someone who has never
                         #been tested, to be LESS than the last time at which
                         #someone was tested who has ever been tested; hence the
                         #use of -Inf instead of Inf.
                         time_last_immunity_event = NA,
                         #using NAs here will trip a debug guard if these are
                         #not updated as they should be when individuals change
                         #immune_status
                         previous_immunity = NA,
                         previous_infection_immunity = NA,
                         isolated = FALSE, #initially
                         Age_Cat = sample(Age_Categories, N, replace = TRUE,
                                          prob = age_probabilities),
                         duration_E = duration_E,
                         duration_IA = sgamma(N, shape=duration_IA_shape, scale=duration_IA_scale),
                         duration_IP = duration_IP,
                         duration_IM = sgamma(N, shape=duration_IM_shape, scale=duration_IM_scale),
                         duration_IS = sgamma(N, shape=duration_IS_shape, scale=duration_IS_scale),
                         duration_IC = sgamma(N, shape=duration_IC_shape, scale=duration_IC_scale),
                         stringsAsFactors = FALSE
                         #"stringsAsFactors = FALSE" is to allow transition into
                         #states that are not present at simulation start
    )

    #pre-calculating all indices for clarity and ease of debugging
    #This block is based on the need to insure that any(index_E & index_IM) ==
    #FALSE . Its complexity is based on the R misfeature that if a == 0, then
    #1:a is not, as one might expect, numeric(), but c(1,0) .
    index_E_or_IM = sample(N, E0 + IM0)
    #There has to be a cleaner way to do this.
    if(E0 > 0 & IM0 > 0) {
        index_E = index_E_or_IM[1:E0]
        index_IM = index_E_or_IM[(E0+1):(E0+IM0)]
    } else if(E0 > 0) {
        index_E = index_E_or_IM
        index_IM = numeric()
    } else if(IM0 > 0) {
        index_E = numeric()
        index_IM = index_E_or_IM
    } else {
        index_E = numeric()
        index_IM = numeric()
    }

    index_R = 1:N %in% sample(N, initial_recovered)
    index_V2 = 1:N %in% sample(N, initial_V2)
    initial_V2_last_five_months = round(ffv_last_five_months * N)
    if(initial_V2_last_five_months == 0) {
        index_V2_last_five_months = rep(FALSE, N)
    } else if(sum(index_V2) == initial_V2_last_five_months) {
        index_V2_last_five_months = index_V2
    } else { #implies sum(index_V2) > 1, so the following sample() call is safe
        index_V2_last_five_months = 1:N %in% sample((1:N)[index_V2],
                                                    initial_V2_last_five_months)
    }
    index_V2_older = index_V2 & !index_V2_last_five_months

    n_V2_older = sum(index_V2_older)
    n_boosted_ever = round(n_V2_older * fraction_boosted_ever)
    n_boosted_last_five_months = round(n_V2_older * fraction_boosted_last_five_months)
    n_boosted_older = n_boosted_ever - n_boosted_last_five_months

    if(n_boosted_ever == 0) {
        index_B = rep(FALSE, N) #= index_B_older = index_B_last_five_months
    } else if(n_boosted_ever == n_V2_older) {
        index_B = index_V2_older
    } else {#implies sum(index_V2_older) > 1, so this sample() call is safe
        index_B = 1:N %in% sample((1:N)[index_V2_older], n_boosted_ever)
    }

    if(n_boosted_last_five_months == 0) {
        index_B_last_five_months = rep(FALSE, N)
    } else if(n_boosted_last_five_months == n_boosted_ever) {
        index_B_last_five_months = index_B
    } else {#implies sum(index_B) > 1, so this sample() call is safe
        index_B_last_five_months = 1:N %in% sample((1:N)[index_B],
                                                   n_boosted_last_five_months)
    }

    index_B_older = index_B & !index_B_last_five_months

    agents$boosting_on_time = ifelse(index_B,
        TRUE,
        ifelse(index_V2_older,
            FALSE,
            sbern(N, fraction_boosted_ever)
        )
    )


    #Note: these can be allowed to not all be N, without eliminating the common
    #random variables benefits, as long as the number of calls remains constant
    #across all interventions.

    #since neither of these is affected by the other factors, it can be left
    #here for now
    agents$infection_status[index_E]= "E"
    agents$time_E[index_E]= -sunif(E0, 0, agents$duration_E[index_E])

    agents$infection_status[index_IM]= "IM"
    agents$time_IM[index_IM]= -sunif(IM0, 0, agents$duration_IM[index_IM])
    # These next two lines shouldn't matter, but seem harmless, and like a good
    # way to avoid the possibility of weird errors.
    agents$time_IP[index_IM] = (agents$time_IM[index_IM] -
                                agents$duration_IP[index_IM])
    agents$time_E[index_IM] = (agents$time_IP[index_IM] -
                                agents$duration_E[index_IM])

    #Using a bunch of uniform distributions. We may wish to change this
    #eventually, but it's fine for now.
    times_V2_last_five_months = times_V2_older = times_V2_B_older = times_V2_B_last_five_months = rep(Inf, N)
    times_V2_last_five_months[index_V2_last_five_months] = sunif(
        initial_V2_last_five_months,
        -boosting_interval,
        0
    )
    times_V2_older[index_V2_older] = sunif(
        initial_V2 - initial_V2_last_five_months,
        -(time_since_first_V2),#fully vax starts in mid-december 2020
        -boosting_interval 
    )
    times_V2_B_older[index_B_older] = sunif(n_boosted_older, -(time_since_first_V2), -(2*boosting_interval+1)) 
    times_V2_B_last_five_months[index_B_last_five_months] = sunif(n_boosted_last_five_months,
                                                     -2*boosting_interval, -boosting_interval)

    times_V2 = ifelse(index_B_older,
        times_V2_B_older,
        ifelse(index_B_last_five_months,
            times_V2_B_last_five_months,
            ifelse(index_V2_last_five_months,
                times_V2_last_five_months,
                ifelse(index_V2_older,
                    times_V2_older,
                    Inf
                )
            )
        )
    )

    

    times_V1 = times_V2 - second_shot_interval
    #again, not perfect, but doesn't actually matter
    #(currently, and probably ever)

    times_B = ifelse(index_B, times_V2 + boosting_interval, Inf)
    
    times_R = rep(Inf, N)
    times_R[index_R]= -sunif(initial_recovered, 0, R_question_period)


    index_V1 = index_V2 #redundant now, but facilitates a slightly more logical coding style

    #NB:
    #   R_first also includes those who have recovered but never been vaccinated,
    #   V1_R_V2 (in principle) also includes those who have recovered and have only been vaccinated once (although these don't actually exist at simulation start in this model),
    #   V2_R_B also includes those who have recovered and have completed primary vaccination,
    #Note also that selective updating (subscripting by the second parameter) occurs within the update_immunity function, for less risk of subscripting error when copy-and-pasting. This is probably somewhat slower, but I have not tested that.
    recovered_complete_protection_time = kConstants$recovered_complete_protection_time
    #Thought: It might be worth creating a generate_update_immunity function so that recovered_complete_protection_time, net_protection, and infection_protection only have to be fed in once
    R_first = times_R < times_V1
    agents = update_immunity(agents, R_first, 'R', times_R, recovered_complete_protection_time, net_protection, infection_protection)
    agents = update_immunity(agents, index_V1, 'V1', times_V1, recovered_complete_protection_time, net_protection, infection_protection)
    V1_R_V2 = times_V1 <= times_R & times_R < times_V2
    agents = update_immunity(agents, V1_R_V2, 'R', times_R, recovered_complete_protection_time, net_protection, infection_protection)
    agents = update_immunity(agents, index_V2, 'V2', times_V2, recovered_complete_protection_time, net_protection, infection_protection)
    V2_R_B = times_V2 <= times_R & times_R < times_B
    agents = update_immunity(agents, V2_R_B, 'R', times_R, recovered_complete_protection_time, net_protection, infection_protection)
    agents = update_immunity(agents, index_B, 'B', times_B, recovered_complete_protection_time, net_protection, infection_protection)
    B_R = times_B < times_R & times_R < Inf
    agents = update_immunity(agents, B_R, 'R', times_R, recovered_complete_protection_time, net_protection, infection_protection)

    #Import text file of disease progression probabilities
    Probability_Matrix <- read.csv('Probability_Matrix.csv')
  
    #Merge parameters for disease progression (Psymptomatic, Psevere, and
    #Pcritical) from table into data frame
    agents=merge.data.frame(agents, Probability_Matrix, by="Age_Cat",
                            all.x=TRUE)

    #pre-calculating, as with times
    agents$p_severe = SEVERE_MULTIPLIER * agents$p_severe
                                                                
    agents = agents[sample(nrow(agents)),] # this is a reshuffled database,
                                           # to randomize order of agents
                                           # reminder: Age_Cat is first column
    return(agents) 
}

# constants.R is part of Food INdustry CoViD Control Tool
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

kConstants = list( #using Google-style notation for now
    isolation_duration = 14, #5, #days
    mu = 5.2, #for log-linear duration_{E+IP}
    sd = .1, #||
    duration_IP_mean = 1.058 * 2.174,
    duration_IP_shape = 1.058,
    duration_IA_mean = 1 * 5,
    duration_IA_shape = 5,
    duration_IM_mean = 0.5 * 16,
    duration_IM_shape = 16,
    duration_IS_mean = 0.4114 * 34.0278,
    duration_IS_shape = 34.0278,
    duration_IC_mean = 0.4114 * 34.0278,
    duration_IC_shape = 34.0278,
    #trying altering p_trans_Ix variables separately
    p_trans_IP = .0575,
    p_trans_IA = .0575 * .11,
    p_trans_IM = .0575 * .44,
    boosting_interval = 152,
    complete_immunity_duration_R = 45,
    second_shot_interval = 21, #days
    max_V1_protection = 0.36,
    V2_ramp_time = 14, #days
    V2_magnitude = 0.9115739,
    V2_decay_rate = 0.08904459 / 7, #1/days; / 7 converts from 1/weeks
    B_ramp_time_1 = 7, # days
    B_ramp_time_2 = 7, # days, between B_ramp_time_1 and end of ramp
    B_mid_ramp_protection = .62,
    B_magnitude_1 = 0.471669758,
    B_magnitude_2 = 0.326600870,
    B_decay_rate_1 = 0.083161719 / 7,
    B_decay_rate_2 = 0.008970573 / 7,
    SEVERE_MULTIPLIER = 1.2,
    R_question_period = 365, #days
    time_since_first_V2 = 365 + 61,
    R_nsp_a = 1.70512,
    R_nsp_b = -0.05211/30.5,
    H_RV12_nsp_a = 3.04736420,
    H_RV12_nsp_b = 0.04724741/30.5,
    H_RB_nsp_a = 4.0685452,
    H_RB_nsp_b = -0.1756493/30.5,
    R_ip_a = 1.2100,
    R_ip_b = -0.1937/30.5,
    H_RV12_ip_a = 1.176188,
    H_RV12_ip_b = -0.125678/30.5,
    H_RB_ip_a = 1.7006945,
    H_RB_ip_b = -0.3068089/30.5,
    hybrid_ramp_time = 30.5,
    recovered_complete_protection_time = 61,
    fraction_ssp_symptomatic = 0.5
)


#Note: Doesn't check for existence of all constants, only those which are used
#in consistency checking (relevant if constructing a set of constants "from
#scratch"
check_consistency = function(constants, altered_single_parameter = NULL, altered_parameters = NULL) {
    p_trans_IP = get('p_trans_IP', constants)
    p_trans_IA = get('p_trans_IA', constants)
    p_trans_IM = get('p_trans_IM', constants)
    if(p_trans_IP + p_trans_IA + p_trans_IM == 0) { #Inversely proportional to contact rates
        return(list(consistent = FALSE, #There isn't really a good solution here
                                        #In general, we are trying to return the
                                        #"least changed" valid set, but that's
                                        #ill-defined
                    fixed = FALSE,
                    fixed_constants = NULL))
    }
    boosting_interval = get('boosting_interval', constants)
    time_since_first_V2 = get('time_since_first_V2', constants)
    if(2 * boosting_interval + 1 > time_since_first_V2) {
        if(!is.null(altered_single_parameter)) { #will eventually want to add
                                                 #altered_parameters logic  as well
            if(altered_single_parameter == 'boosting_interval') {
                constants[['boosting_interval']] = (time_since_first_V2 - 1) / 2
                return(list(consistent = FALSE,
                            fixed = TRUE,
                            fixed_constants = constants))
            } else if(altered_single_parameter == 'time_since_first_V2') {
                constants[['time_since_first_V2']] = 2 * boosting_interval + 1
                return(list(consistent = FALSE, #ultimately, we ought to allow clear identification of what is wrong
                            fixed = TRUE,
                            fixed_constants = constants))
            } else {
                return(list(consistent = FALSE, #ultimately, we ought to allow clear identification of what is wrong
                            fixed = FALSE,
                            fixed_constants = NULL))
            }
        } else {
            return(list(consistent = FALSE,
                        fixed = FALSE,
                        fixed_constants = NULL))
        }
    }
    return(list(consistent = TRUE,
                fixed = TRUE,
                fixed_constants = constants))
}


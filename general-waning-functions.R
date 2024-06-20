# general-waning-functions.R is part of Food INdustry CoViD Control Tool
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


make_V1_protection_general = function (kConstants) {
    second_shot_interval = get('second_shot_interval', kConstants)
    max_V1_protection = get('max_V1_protection', kConstants)
    function(t, prev) {
    #prev unused; included for consistent interface
        ifelse(t < second_shot_interval,
            max_V1_protection * t / second_shot_interval,
            max_V1_protection
        )
    }
}

make_V2_protection_general = function(kConstants) {
    V2_ramp_time = get('V2_ramp_time', kConstants)
    V2_magnitude = get('V2_magnitude', kConstants)
    V2_decay_rate = get('V2_decay_rate', kConstants)
    V2_decay = function(t) {
        pmin(V2_magnitude * exp(-V2_decay_rate * t), 1)
        #pmin(..., 1) is a kludge to handle 0.5, 1.5 sensitivity testing
        #otherwise, maximum protection would be > 1
    }
    V2_max = V2_decay(V2_ramp_time)
    function(t, prev) {
        ifelse(t < V2_ramp_time,
            (t / V2_ramp_time) * V2_max +
                ((V2_ramp_time - t) / V2_ramp_time) * prev,
            V2_decay(t)
        )
    }
}

make_B_decay_general = function(kConstants) {
    B_magnitude_1 = get('B_magnitude_1', kConstants)
    B_magnitude_2 = get('B_magnitude_2', kConstants)
    B_decay_rate_1 = get('B_decay_rate_1', kConstants)
    B_decay_rate_2 = get('B_decay_rate_2', kConstants)
    function(t) {
        (B_magnitude_1 * exp(-B_decay_rate_1 * t) +
            B_magnitude_2 * exp(-B_decay_rate_2 * t))
    }
}

make_B_protection_general = function(B_decay, kConstants) {
    #B_decay is passed in as an already-constructed parameter because it is used
    #both here and in make_R_protection_general
    #NB: Changes here may make for some weird results
    B_ramp_time_1 = get('B_ramp_time_1', kConstants)
    B_ramp_time_2 = get('B_ramp_time_2', kConstants)
    B_mid_ramp_protection = get('B_mid_ramp_protection', kConstants)
    total_B_ramp_time = B_ramp_time_1 + B_ramp_time_2
    B_max = B_decay(total_B_ramp_time)
    function(t, prev) {
        ifelse(t < B_ramp_time_1, 
            t / B_ramp_time_1 * B_mid_ramp_protection +
                (B_ramp_time_1 - t) / B_ramp_time_1 * prev,
            ifelse(t < total_B_ramp_time,
                (total_B_ramp_time - t) / B_ramp_time_2 * B_mid_ramp_protection +
                    (t - B_ramp_time_1) / B_ramp_time_2 * B_max, 
                B_decay(t)
            )
        )
    }
}

make_R_protection_general = function(B_decay, kConstants) {
    complete_immunity_duration_R = get('complete_immunity_duration_R',
                                       kConstants)
    B_ramp_time_1 = get('B_ramp_time_1', kConstants)
    B_ramp_time_2 = get('B_ramp_time_2', kConstants)
    total_B_ramp_time = B_ramp_time_1 + B_ramp_time_2
    
    R_decay = function(t) {
        pmin(1, B_decay(t - complete_immunity_duration_R + total_B_ramp_time))
    }

    function(t, prev) {
    #prev unused; included for consistent interface
        ifelse(t < complete_immunity_duration_R,
            1,
            R_decay(t)
        )
    }
}

#R_protection = B_protection
#default_B_protection = default_R_protection

make_protection_functions = function(V1_protection, V2_protection, B_protection,
                                     R_protection,
                                     kConstants
                                     ) {
    #ignoring R_protection and hardcoding the rest for now
    #leaving V1_protection, V2_protection, and B_protection alone for now
    with(kConstants,
        {
        make_logistic_protection = function(a, b) {
            f = function(t) exp(a + b * t) / (1 + exp(a + b * t))
            f_max = f(hybrid_ramp_time) #whatever
            
            function(t, t_R, prev) {
                ifelse(t_R < recovered_complete_protection_time,
                    1,
                    ifelse(t < hybrid_ramp_time,
                        ((hybrid_ramp_time - t) * prev + t * f_max) / hybrid_ramp_time,
                        f(t)
                    )
                )
            }
        }

        R_nsp <<- make_logistic_protection(R_nsp_a, R_nsp_b)
        H_V1_nsp <<- make_logistic_protection(H_RV12_nsp_a, H_RV12_nsp_b)
        H_V2_nsp <<- make_logistic_protection(H_RV12_nsp_a, H_RV12_nsp_b)
        H_B_nsp <<- make_logistic_protection(H_RB_nsp_a, H_RB_nsp_b)
    
        R_ip <<- make_logistic_protection(R_ip_a, R_ip_b)
        H_V1_ip <<- make_logistic_protection(H_RV12_ip_a, H_RV12_ip_b)
        H_V2_ip <<- make_logistic_protection(H_RV12_ip_a, H_RV12_ip_b)
        H_B_ip <<- make_logistic_protection(H_RB_ip_a, H_RB_ip_b)
        }
    )

    net_severity_protection = function(agents, start_time) {
        ais = agents$immune_status
        tryCatch(
            {
                t = (start_time - agents$time_last_immunity_event)
                t_R = (start_time - agents$time_R)
            },
            error = function(e) {
                print(e)
                browser()
            }
        )
        #t = (start_time - agents$time_last_immunity_event)
        t = pmax(t, 0) #TBD (eventually): find a way to not need this kludge
        prev = agents$previous_immunity
        protection = ifelse(ais == 'FS',
            0,
            ifelse(ais %in% c('V1', 'V2', 'B'),
                NA, #kludge: this are handled in a different fashion
                ifelse(ais == 'R',
                    R_nsp(t, t_R, prev),
                    ifelse(ais == 'H_V1',
                        H_V1_nsp(t, t_R, prev),
                        ifelse(ais == 'H_V2', 
                            H_V2_nsp(t, t_R, prev),
                            ifelse(ais == 'H_B',
                                H_B_nsp(t, t_R, prev),
                                NaN
                            )
                        )
                    )
                )
            )
        )
        if(any(is.nan(protection))) {
            mask = is.nan(protection)
            cat('Problematic agents:\n')
            print(agents[mask,])
            stop('NaNs in net_severity_protection') #debugging
        }
        protection
    }

    RH_infection_protection = function(agents, start_time) {
        ais = agents$immune_status
        tryCatch(
            {
                t = (start_time - agents$time_last_immunity_event)
                t_R = (start_time - agents$time_R)
            },
            error = function(e) {
                print(e)
                browser()
            }
        )
        t = pmax(t, 0) #A kludge designed to avoid "negative ramp" issues when
                       #start_time is the start of an interval, and
                       #time_last_immunity_event occcurs within that interval.
                       #This was required in an earlier version of the code, due
                       #to how vaccination times were handled.
                       #TBD: Check if anything goes wrong if it's removed now.
        prev = agents$previous_infection_immunity
        ifelse(ais %in% c('FS', 'V1', 'V2', 'B'),
            NA, #handled differently
            ifelse(ais == 'R',
                R_ip(t, t_R, prev),
                ifelse(ais == 'H_V1',
                    H_V1_ip(t, t_R, prev),
                    ifelse(ais == 'H_V2', 
                        H_V2_ip(t, t_R, prev),
                        ifelse(ais == 'H_B',
                            H_B_ip(t, t_R, prev),
                            NA
                        )
                    )
                )
            )
        )   
    }

    RH_symptomseverity_protection = function(agents, start_time) {
        nsp = net_severity_protection(agents, start_time)
        ip = RH_infection_protection(agents, start_time)
        #browser()
        ifelse(nsp == 1, 1, 1 - (1 - nsp) / (1 - ip))
    }

    RH_symptom_protection = function(agents, start_time) {
        ssp = RH_symptomseverity_protection(agents, start_time)
        fraction_ssp_symptomatic = kConstants$fraction_ssp_symptomatic
        ifelse(ssp == 1, 1, 1 - (1 - ssp)^fraction_ssp_symptomatic)
    }

    RH_severity_protection = function(agents, start_time) {
        ssp = RH_symptomseverity_protection(agents, start_time)
        fraction_ssp_symptomatic = kConstants$fraction_ssp_symptomatic
        ifelse(ssp == 1, 1, 1 - (1 - ssp)^(1 -fraction_ssp_symptomatic))
    }

    net_symptomatic_protection = function(agents, start_time) {
        ais = agents$immune_status
        t = (start_time - agents$time_last_immunity_event)
        t = pmax(t, 0) #TBD (eventually): find a way to not need this kludge
        prev = agents$previous_immunity
        protection = ifelse(ais == 'FS',
            0,
            ifelse(ais == 'V1',
                V1_protection(t, prev),
                ifelse(ais == 'V2',
                    V2_protection(t, prev),
                    ifelse(ais == 'B',
                        B_protection(t, prev),
                        ifelse(ais %in% c('R', 'H_V1', 'H_V2', 'H_B'),
                            NA, #handled a different way
                            NaN #just not valid
                        )
                    )
                )
            )
        )
        if(any(is.nan(protection))) {
            mask = is.na(protection)
            cat('Problematic agents:\n')
            print(agents[mask,])
            stop('NAs in net_symptomatic_protection') #debugging
        }
        protection
    }

    infection_protection = function(agents, start_time) {
        ais = agents$immune_status
        ifelse(ais %in% c('FS', 'V1', 'V2', 'B'),
            1 - sqrt(1 - net_symptomatic_protection(agents, start_time)),
            RH_infection_protection(agents, start_time)
        )   
    }

    symptom_protection = function(agents, start_time) {
        nsp = net_symptomatic_protection(agents, start_time)
        ip = infection_protection(agents, start_time)
        ifelse(agents$immune_status %in% c('FS', 'V1', 'V2', 'B'),
            ifelse(nsp == 1, 1, 1 - (1 - nsp) / (1 - ip)),
            RH_symptom_protection(agents, start_time)
        )
    }

    severity_protection = function(agents, start_time) {
        ifelse(agents$immune_status %in% c('FS', 'V1', 'V2', 'B'),
            symptom_protection(agents, start_time), 
            RH_severity_protection(agents, start_time)
        )
    }

    net_protection = function(agents, start_time) { #for the purpose of setting previous_immunity
                                                    #in turn, for the purpose of constructing ramps
                                                    #not a valid measure of relative protection between the first and second types
        ifelse(agents$immune_status %in% c('FS', 'V1', 'V2', 'B'),
            net_symptomatic_protection(agents, start_time),
            net_severity_protection(agents, start_time)
        )
    }

    list(net_protection = net_protection,
         infection_protection = infection_protection,
         symptom_protection = symptom_protection,
         severity_protection = severity_protection
    )
}

make_protection_functions_general = function(kConstants) {
    B_decay = make_B_decay_general(kConstants)
    make_protection_functions(
        make_V1_protection_general(kConstants),
        make_V2_protection_general(kConstants),
        make_B_protection_general(B_decay, kConstants),
        NULL,
        kConstants
    )
}


two_level_protection = function(level, duration) {
    function(t, prev) {
    #prev unused; included for consistent interface
        ifelse(t < duration,
            level,
            0
        )
    }
}

#For use when comparing results with v1.1.3
"make_one_one_three = function() {
    V1_net_symptoms = 1 - .37
    V2_susceptibility = 1 - .65
    V2_net_symptoms = 1 - .88
    V2_symptoms = V2_net_symptoms / V2_susceptibility
    V2_exp = log(V2_symptoms) / log(V2_net_symptoms)
    V1_symptoms = V1_net_symptoms ** V2_exp #haven't seen them separated, so doing the best I can
    V1_susceptibility = V1_net_symptoms / V1_symptoms

    protection_functions = make_protection_functions(
        two_level_protection(1 - V1_net_symptoms, Inf),
        two_level_protection(1 - V2_net_symptoms, Inf),
        two_level_protection(1, Inf), #Since the only way to become boosted
                                      #without booster shots is vax + R
        two_level_protection(1, Inf)
    )

    net_symptomatic_protection = protection_functions$net_symptomatic_protection

    infection_protection = function(agents, start_time) {
        ais = agents$immune_status
        t = (start_time - agents$time_last_immunity_event)
        t = pmax(t, 0) #TBD (eventually): find a way to not need this kludge
        prev = agents$previous_immunity
        protection = ifelse(ais == 'FS',
            0,
            ifelse(ais == 'V1',
                1 - V1_susceptibility,
                ifelse(ais == 'V2',
                    1 - V2_susceptibility,
                    ifelse(ais == 'B',
                        1,
                        ifelse(ais == 'R',
                            1,
                            NA
                        )
                    )
                )
            )
        )
        if(any(is.na(protection))) {
            mask = is.na(protection)
            cat('Problematic agents:\n')
            print(agents[mask,])
            stop('NAs in net_symptomatic_protection') #debugging
        }
        protection
    }

    symptom_protection = function(agents, start_time) {
        nsp = net_symptomatic_protection(agents, start_time)
        ip = infection_protection(agents, start_time)
        ifelse(nsp == 1, 1, 1 - (1 - nsp) / (1 - ip))
    }

    protection_functions$infection_protection = infection_protection
    protection_functions$symptom_protection = symptom_protection

    protection_functions
}

one_one_three_protection_functions = make_one_one_three()"

make_protection_functions_2020 = function(V1_protection, V2_protection, B_protection,
                                         R_protection) {
    net_symptomatic_protection = function(agents, start_time) {
        ais = agents$immune_status
        t = (start_time - agents$time_last_immunity_event)
        t = pmax(t, 0) #TBD (eventually): find a way to not need this kludge
        prev = agents$previous_immunity
        protection = ifelse(ais == 'FS',
            0,
            ifelse(ais == 'V1',
                V1_protection(t, prev),
                ifelse(ais == 'V2',
                    V2_protection(t, prev),
                    ifelse(ais == 'B',
                        B_protection(t, prev),
                        ifelse(ais %in% c('R', 'H_V1', 'H_V2', 'H_B'),
                            R_protection(t, prev),
                            NaN
                        )
                    )
                )
            )
        )
        if(any(is.nan(protection))) {
            mask = is.na(protection)
            cat('Problematic agents:\n')
            print(agents[mask,])
            stop('NAs in net_symptomatic_protection') #debugging
        }
        protection
    }

    infection_protection = function(agents, start_time) {
        1 - sqrt(1 - net_symptomatic_protection(agents, start_time))
    }

    symptom_protection = function(agents, start_time) {
        ifelse(nsp == 1, 1, 1 - (1 - nsp) / (1 - ip))
    }

    severity_protection = function(agents, start_time) {
        0
        #ifelse(agents$immune_status %in% c('FS', 'V1', 'V2', 'B'),
        #    symptom_protection(agents, start_time), #TBD 2023-06-27: make this less kludgey
        #    RH_severity_protection(agents, start_time)
        #)
    }

    net_protection = function(agents, start_time) { #for the purpose of setting previous_immunity
                                                    #in turn, for the purpose of constructing ramps
                                                    #not a valid measure of relative protection between the first and second types
        #ifelse(agents$immune_status %in% c('FS', 'V1', 'V2', 'B'),
            net_symptomatic_protection(agents, start_time)#,
        #    net_severity_protection(agents, start_time)
        #)
    }

    list(net_protection = net_protection,
         infection_protection = infection_protection,
         symptom_protection = symptom_protection,
         severity_protection = severity_protection,
         net_symptomatic_protection = net_symptomatic_protection
    )
}


make_2020 = function() {
    V1_susceptibility = 1 - .8
    V2_susceptibility = 1 - .9
    V1_symptoms = (1 - .88) / (1 - .8) #purely default; if I've seen an estimate of this, I do not recall it
    V2_symptoms = (1 - .94) / (1 - .9)
    V1_net_symptoms = 1-.88
    V2_net_symptoms = 1-.94

    protection_functions = make_protection_functions_2020(
        two_level_protection(1 - V1_net_symptoms, Inf),
        two_level_protection(1 - V2_net_symptoms, Inf),
        two_level_protection(1, Inf), #Since the only way to become boosted
                                      #without booster shots is vax + R
        two_level_protection(1, Inf)
    )

    net_symptomatic_protection = protection_functions$net_symptomatic_protection

   infection_protection = function(agents, start_time) {
        ais = agents$immune_status
        t = (start_time - agents$time_last_immunity_event)
        t = pmax(t, 0) #TBD (eventually): find a way to not need this kludge
        prev = agents$previous_immunity
        protection = ifelse(ais == 'FS',
            0,
            ifelse(ais == 'V1',
                1 - V1_susceptibility,
                ifelse(ais == 'V2',
                    1 - V2_susceptibility,
                    ifelse(ais == 'B',
                        1,
                        ifelse(ais == 'R',
                            1,
                            NA
                        )
                    )
                )
            )
        )
        if(any(is.na(protection))) {
            mask = is.na(protection)
            cat('Problematic agents:\n')
            print(agents[mask,])
            stop('NAs in net_symptomatic_protection') #debugging
        }
        protection
    }

    symptom_protection = function(agents, start_time) {
        nsp = net_symptomatic_protection(agents, start_time)
        ip = infection_protection(agents, start_time)
        ifelse(nsp == 1, 1, 1 - (1 - nsp) / (1 - ip))
    }

    protection_functions$infection_protection = infection_protection
    protection_functions$symptom_protection = symptom_protection

    protection_functions
}

protection_2020 = make_2020()


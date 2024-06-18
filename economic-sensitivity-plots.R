source('constants.R')
library(abind)

dir.create('figures-2024-04-16/')

N = 103

eConstants = list(
    threshold = 0.85,
    exponent = 0.437,

    thermometer_cost_each = 20,    # $20 per thermometer
    KN95_cost = 1,                 # $1 per mask per day
    face_shield_cost = 3,          # $3 per face shield. Changing every 30 days. ($0.1/day)
    ts_time = 3,                   # 3 seconds for each screening
    ts_limit = 5,                  #screening should be completed under 5 minute

    vt_kit = 10,                   # $10 per test
    vt_time = 1/4,                 # 15 minutes waiting assumed

    vaccination_time = 0.75,        # in hours, per person

    air_cleaner = 1000,            # 1 air cleaner per 1000 sqft
    life = 3 * 365                 # 3year life of air_cleaner
)

#NB: This will only work if the data has new_internal_infections, and will only
#be *right* if the destination of transmission has been changed from E to R, as
#in branch sensitivity-r0
r_eff = function(df, mask, limited_i, output_per_shift, hourly_wage, eConstants) {
    df; mask; limited_i; output_per_shift; hourly_wage; eConstants
    apply(df[,'new_internal_infections',], 2, sum)
}

symptomatic_infections = function(df, mask, limited_i,
    output_per_shift, hourly_wage, eConstants
) { 
    df; mask; limited_i; output_per_shift; hourly_wage; eConstants
    apply(df[,'new_symptomatic_infections',], 2, sum)
}

shiftwise_unavailable = function(df, mask, limited_i,
    output_per_shift, hourly_wage, eConstants
) {
    df; mask; limited_i; output_per_shift; hourly_wage; eConstants
    apply(df[,'qn_absent',], 2, sum)
}

shiftwise_scheduled = function(data) {
    data[,'qn_scheduled',]
}

shiftwise_unavailable_fraction = function(df) {
    #print(summary(df))
    #print(dimnames(df))
    df[,'qn_absent',] / df[,'qn_scheduled',]
}

#shiftwise_short = function(data) {
#    shiftwise_unavailable_fraction(data) > .15
#}

shiftwise_production_loss = function(df, mask,
    output_per_shift, eConstants
) {
    df; mask; output_per_shift; eConstants
    #print(summary(df))
    #print(dimnames(df))
    #print('---')
    #print(mask)
    #df = df[mask,,]
    runs = dim(mask)[2]
    df = abind(sapply(1:runs, function(run) df[mask[,run],, run][1:128,], simplify = FALSE), along = 3)

    fraction_available = 1 - shiftwise_unavailable_fraction(df)
    #print(mean(fraction_available))
    #print(eConstants$threshold)
    #browser()
    adjusted_fraction_available = pmin(fraction_available / eConstants$threshold, 1)
    fractional_production = adjusted_fraction_available^eConstants$exponent
    fractional_loss = 1 - fractional_production

    #apply(
    fractional_loss * output_per_shift
    #, 2, sum)
}


limited_runs_index = c(1,2,4,9,13)


temperature_screening_cost = function(data,
    mask,
    hourly_wage,
    eConstants
) {
    data; mask; hourly_wage; eConstants
    scheduled = shiftwise_scheduled(data)
    #we can use NA for limited_i and output_per_shift because they're not actually used there
    #available = scheduled - shiftwise_unavailable(data, mask, NA,
    #                                              NA, hourly_wage,
    #                                              eConstants)
    available = scheduled - data[,'qn_absent',]
    screeners = ceiling(scheduled) / (eConstants$ts_limit * 60 / eConstants$ts_time)
    ts_time <- available * eConstants$ts_time / screeners / 3600   # Actual daily screening time in hours
        #the above should be renamed, for clarity
    compensation <- ts_time * screeners * hourly_wage * 2 # have to pay the screeners, and the people being screened
    screener_training_cost = ceiling(N/100) * hourly_wage #max(screeners) * hourly_wage # 1hour training cost for screeners
    thermometer_cost <- max(screeners) * eConstants$thermometer_cost_each

    initial_cost = screener_training_cost + thermometer_cost
    ongoing_cost = compensation + (eConstants$KN95_cost + eConstants$face_shield_cost/30) * screeners

    ongoing_cost[1] = ongoing_cost[1] + initial_cost

    #browser()

    ifelse(is.na(ongoing_cost), 0, ongoing_cost) #needs modification if we ever end up plotting over time
}

virus_testing_cost = function(data,
    mask,
    hourly_wage,
    eConstants
) {
    data; hourly_wage; eConstants
    #array(0, c(dim(data)[1], dim(data)[3]))
    #vt_prod <- output_per_week / 5 / 8 * vt_time #production value during 15 min ts_time (5days/wk, 8hr/day)

    # Average wage compensation + kit cost over simulation
    vt <- data[,'tests',] # number of tests
    vt_cost <- vt * (eConstants$vt_time * hourly_wage + eConstants$vt_kit) # total cost

    vt_cost
}

vaccination_cost = function(data, mask, hourly_wage, eConstants) {
    data; hourly_wage; eConstants
    #array(0, c(dim(data)[1], dim(data)[3]))
    ############## Vaccination ############
    # 0.75 hour paid sick leave per vaccination  
    # no production loss
    data[,'doses',] * hourly_wage * eConstants$vaccination_time
}

R0_reduction_cost = function(data, kludge_index,
    mask, hourly_wage, eConstants
) {
    data; kludge_index; hourly_wage; eConstants
    "bi_available <- shiftwise_scheduled(data) - shiftwise_unavailable(data,
                                                                      mask,
                                                                      limited_i,
                                                                      output_per_shift,
                                                                      hourly_wage,
                                                                      eConstants)"
    bi_available = data[,'qn_scheduled',] - data[,'qn_absent',]
    bi_available = ifelse(is.na(bi_available), 0, bi_available)
#    array(0, c(dim(data)[1], dim(data)[3]))
    if(kludge_index == 8) {
        bi_cost = eConstants$KN95_cost * bi_available
    } else if(kludge_index == 9 || (kludge_index == 10 && farm_or_facility == 'farm')) {
        bi_cost = ((eConstants$KN95_cost + eConstants$face_shield_cost/30) * bi_available)
    } else if(kludge_index == 10) {
        #TBD: introduce size, if it ever matters
        bi_cost = ((eConstants$KN95_cost + eConstants$face_shield_cost/30) * bi_available) + size/1000 * eConstants$air_cleaner / eConstants$life 
    } else {
        stop(kludge_index)
    }
    #print(dim(bi_cost))
    #print(dim(bi_cost))

    bi_cost
}

intervention_expenses_function = function(data, limited_i,
    output_per_shift, mask, hourly_wage, eConstants
) {
    data; limited_i; output_per_shift; hourly_wage; eConstants
    i = limited_runs_index[limited_i]
    #print('in')
    if(i == 1) {
        array(0, c(dim(data)[1], dim(data)[3]))
    } else if(i == 2) {
        #print('TEMP')
        temperature_screening_cost(data, mask, hourly_wage, eConstants)
        #print('PMET')
    } else if(i %in% 3:5) {
        virus_testing_cost(data, mask, hourly_wage, eConstants)
    } else if(i %in% c(6:7, 11:13)) {
        vaccination_cost(data, mask, hourly_wage, eConstants)
    } else {
        R0_reduction_cost(data, i, mask, hourly_wage, eConstants)
        #print(dim(x))
        #x
    }
}

g = function(data, mask, limited_i,
    output_per_shift, hourly_wage, eConstants
) {
    data; mask; limited_i; output_per_shift; hourly_wage; eConstants
    #cat('\nlimited_i:', limited_i, '\noutput_per_shift:', output_per_shift, '\nmask:', mask, '\nhourly_wage:', hourly_wage, '\n\n')
    fd = shiftwise_production_loss(data, mask, output_per_shift, eConstants)
    r = intervention_expenses_function(data, limited_i, output_per_shift, mask, hourly_wage, eConstants)
    #print(dim(r))
    #if(limited_i != 1) browser()
    #browser()
    total = apply(fd, 2, sum) + apply(r, 2, sum)
    "cat(
        '\n\n\nlimited_i:', limited_i, '\ni:', limited_runs_index[limited_i],
        # '\nfd:', fd,
        '\napply(fd, 2, sum):', apply(fd, 2, sum),
        '\nmean(apply(fd, 2, sum)):', mean(apply(fd, 2, sum)),
        # '\nr:', r,
        '\napply(r, 2, sum):', apply(r, 2, sum), '\nmean(apply(r, 2, sum)):',
        mean(apply(r, 2, sum)), '\ntotal:', total, '\n\n\n'
    )"
    total
}

c4 = c('black', 'blue3', 'lightblue1', 'red2', 'gray80', 'darkgreen', 'yellow2')

colors = c(
    c4[1],
    c4[2],
    c4[3], c4[3], c4[3],
    c4[4], c4[4],
    c4[5], c4[5], c4[5],
    c4[6], c4[6],
    c4[7]
)[limited_runs_index]

colors[3] = 'red2'
colors[4] = 'darkgreen'

row.names<-c(
    "Baseline",
    "Temperature Screening, 38.0°C",
    "Virus Test, p = 0.05 / Working Day",
    "Virus Test, p = 0.3 / Working Day",
    "Virus Test, p = 1.0 / Working Day",
    "Vaccination, p = 0.02 / Day",
    "Vaccination, p = 0.04 / Day",
    "Phys. Dist./Biosafety: -20% R₀",
    "Phys. Dist./Biosafety: -40% R₀",
    "Phys. Dist./Biosafety: -80% R₀",
    'Boosting, p = 0.02 / day',
    'Boosting, p = 0.04 / day',
    'Vax + Boosting, p = 0.02/day'
)[limited_runs_index]


rds_filename = function(
    prepended_key, index_i, index_j,
    max_j,
    folder_name,
    unique_ids,
    is_baselines,
    community_transmissions,
    work_R0s,
    dormitory_R0s,
    E0,
    initial_recovereds,
    initial_V2s,
    n_sims
) { #j being the index to all the parameters that are plurals
    prepended_key; index_i; index_j; max_j; folder_name; unique_ids; is_baselines; community_transmissions; work_R0s; dormitory_R0s; E0; initial_recovereds; initial_V2s; n_sims
    work_id_multiplier = ifelse(index_i == 4,
        'x(1-0.4)',
        ''
    )
    testing_string = ifelse(index_i == 2,
        ',T.test-38',
        ifelse(index_i == 3,
            ',v.test-0.3-rational',
            ''
        )
    )
    vax_string = ifelse(index_i == 5,
        ',vax-rate0.02',
        ''
    )
    dormitory_R0_string = ifelse(dormitory_R0s[index_j] == 0,
        '',
        paste0(',dormitory_R0-', dormitory_R0s[index_j])
    )
    baseline_string = ifelse(is_baselines[index_j],
        'baseline',
        ''
    )
    initial_V2_string = ifelse(initial_V2s[index_j] == 0,
        '',
        paste0(',initial_V2-', initial_V2s[index_j])
    )
    initial_recovered_string = ifelse(initial_recovereds[index_j] == 0,
        '',
        paste0(',initial_recovered-', initial_recovereds[index_j])
    )

    paste0(
        folder_name, '/', unique_ids[index_j],
        prepended_key,
        baseline_string,
        '_community-', community_transmissions[index_j],
        ',work_R0-', work_R0s[index_j],
        work_id_multiplier,
        dormitory_R0_string,
        ',E0-', E0,
        vax_string,
        testing_string,
        initial_recovered_string,
        initial_V2_string,
        ',n_sims-', n_sims,
        'index_i-', index_i,
        '_full-output.rds'
    )
}

make_batch = function(d, key, index_j,
    max_j,
    folder_name,
    unique_ids,
    is_baselines,
    community_transmissions,
    work_R0s,
    dormitory_R0s,
    E0,
    initial_recovereds,
    initial_V2s,
    n_sims,
    summary_names,
    summary_fns, #CONTINUE HERE
    masks,
    output_per_shifts,
    hourly_wages,
    eConstants
) {
    d; key; index_j; max_j; folder_name; unique_ids; is_baselines; community_transmissions; work_R0s; dormitory_R0s; E0; initial_recovereds; initial_V2s; n_sims; summary_names; summary_fns; masks; output_per_shifts; hourly_wages; eConstants
    prepended_key = ifelse(key == '',
        '',
        paste0('-', key)
    )
    for(i in 1:5) {
        #print('Confirm')
        #cat(prepended_key, ':', i, ':', dormitory_R0s, '\n')
        filename = rds_filename(prepended_key, i, index_j, max_j, folder_name,
                                unique_ids, is_baselines,
                                community_transmissions, work_R0s,
                                dormitory_R0s, E0, initial_recovereds,
                                initial_V2s, n_sims)
        #print('Deconfirm')
        data_ = readRDS(filename)
        key_ = paste0(i, key)
        d[[key_]] = list()
        #print(summary_names)
        #print(summary_fns)
        for(k in 1:length(summary_names)) {
            #print('chnug')
            #print(masks)
            #print(masks[index_j])
            #print('gunch')
            d[[key_]][[summary_names[k]]] = summary_fns[[k]](
                data_, masks[[index_j]], i, output_per_shifts[[index_j]],
                hourly_wages[[index_j]], eConstants
            )
            #cat('outcomes:', d[[key_]][[summary_names[k]]], '\n')
        }
    }
    d
}


get_real_multiplier = function(sensitivity_variable,
                               theoretical_multiplier,
                               kConstants) {
    sensitivity_variable; theoretical_multiplier; kConstants
    kConstants_ = kConstants
    kConstants_[[sensitivity_variable]] = theoretical_multiplier * kConstants_[[sensitivity_variable]]
    ccl = check_consistency(kConstants_, altered_single_parameter = sensitivity_variable)
    kConstants_fixed = get('fixed_constants', ccl)
    if(!get('consistent', ccl) && !get('fixed', ccl)) {
        stop('Unfixable constants')
    }
    get(sensitivity_variable, kConstants_fixed) / get(sensitivity_variable, kConstants)
}

make_dd = function(
    max_j,
    sensitivity_multipliers,
    kConstants,
    folder_name,
    unique_ids,
    is_baselines,
    community_transmissions,
    work_R0s,
    dormitory_R0s,
    E0,
    initial_recovereds,
    initial_V2s,
    n_sims,
    summary_names,
    summary_fns,
    masks,
    output_per_shifts, hourly_wages, eConstants
) {
    max_j; sensitivity_multipliers; kConstants; folder_name; unique_ids; is_baselines; community_transmissions; work_R0s; dormitory_R0s; E0; initial_recovereds; initial_V2s; n_sims; summary_names; summary_fns; masks; output_per_shifts; hourly_wages; eConstants
    dd = list()
    for(index_j in 1:max_j) {
        #cat('\n\n', index_j, ':', dormitory_R0s, '\n\n')
        dd[[index_j]] = make_batch(list(), '', index_j, max_j, folder_name,
                                   unique_ids, is_baselines,
                                   community_transmissions, work_R0s,
                                   dormitory_R0s, E0, initial_recovereds,
                                   initial_V2s, n_sims, summary_names, summary_fns,
                                   masks, output_per_shifts, hourly_wages,
                                   eConstants)
    }

    for(sensitivity_variable in names(kConstants)) {
        real_multipliers = sapply(sensitivity_multipliers, function(m) get_real_multiplier(sensitivity_variable, m, kConstants))
        for(sensitivity_multiplier in real_multipliers) {
            if(sensitivity_multiplier != 1) {
                key = paste0(sensitivity_variable, '-', sensitivity_multiplier)
                for(index_j in 1:max_j) {
                    dd[[index_j]] = make_batch(dd[[index_j]], key, index_j,
                                               max_j, folder_name, unique_ids,
                                               is_baselines,
                                               community_transmissions,
                                               work_R0s, dormitory_R0s, E0,
                                               initial_recovereds, initial_V2s,
                                               n_sims, summary_names,
                                               summary_fns, masks,
                                               output_per_shifts, hourly_wages,
                                               eConstants)
                }
            }
        }
    }
    dd
}


make_paneled_plot = function(filename, outcome_name, ylab, dd, kConstants,
                             sensitivity_multipliers, max_j,
                             csv_filename, unique_ids) {#, selection_mode) {  ######
    #print(outcome_name)
    filename; outcome_name; ylab; dd; kConstants; sensitivity_multipliers; max_j#; selection_mode
    png(filename, height = 200*5, width = 200*7)
    #layout(matrix(c(1:29, 34, 30:34), ncol = 7))
    layout(matrix(c(1:49), ncol = 7))
    #figuring out bounds:
    greatest_positive_difference = 0
    greatest_negative_difference = 0
    variables_to_exclude = list()
    values_df = data.frame(parameter_set = character(), sensitivity_variable = character(), multiplier = numeric(), intervention = character(), value = numeric()) ######
    for(sensitivity_variable in names(kConstants)) {
        real_multipliers = sapply(
            sensitivity_multipliers,
            function(m) get_real_multiplier(sensitivity_variable, m, kConstants)
        )

        for(i in 1:5) {
            intervention = row.names[i] ######
            null_value = dd[[1]][[paste0(i)]][[outcome_name]]
            keys = sapply(
                real_multipliers,
                function(m) {
                    ifelse(m == 1,
                        paste0(i),
                        paste0(i, sensitivity_variable, '-', m)
                    )
                }
            )
            for(j in 1:max_j) {
                values = sapply(keys, function(key) dd[[j]][[key]][[outcome_name]])
                for(multiplier in real_multipliers) { ######
                    if(multiplier == 1) { ######
                        key = paste0(i) ######
                    } else { ######
                        key = paste0(i, sensitivity_variable, '-', multiplier) ######
                    } ######
                    value = dd[[j]][[key]][[outcome_name]] ######
                    values_df = rbind(values_df, data.frame(parameter_set = unique_ids[j], sensitivity_variable = sensitivity_variable, multiplier = multiplier, intervention = intervention, value = value)) ######
                } ######
                #print('on')
                #print(values)
                #print('off')
                #browser()
                this_greatest_positive_difference = max(0, log(values[3] / values[2]), log(values[1] / values[2]), na.rm = TRUE)
                this_greatest_negative_difference = min(0, log(values[3] / values[2]), log(values[1] / values[2]), na.rm = TRUE)
                if(this_greatest_positive_difference >= greatest_positive_difference) {
                    #print('Positive')
                    #print(values)
                    #print(this_greatest_positive_difference)
                    #print('End positive')
                    if(this_greatest_positive_difference == Inf) {
                        cat('\nPositive Infinite:\n', sensitivity_variable, '\n', i, '\n', j, '\n', values, '\n\n')
                        variables_to_exclude = c(variables_to_exclude, sensitivity_variable)
                    } else {
                        greatest_positive_difference = this_greatest_positive_difference
                    }
                }
                if(this_greatest_negative_difference <= greatest_negative_difference) {
                    #print('Negative')
                    #print(values)
                    #print(this_greatest_negative_difference)
                    #print('End negative')
                    if(this_greatest_negative_difference == Inf) {
                        cat('\nNegative Infinite:\n', sensitivity_variable, '\n', i, '\n', j, '\n', values, '\n\n')
                        variables_to_exclude = c(variables_to_exclude, sensitivity_variable)
                    } else {
                        greatest_negative_difference = this_greatest_negative_difference
                    }
                }
            }
        }
    }
    write.csv(values_df, csv_filename) ######

#print(greatest_negative_difference)
#print(greatest_positive_difference)

    greatest_differences = c()
    greatest_difference_indices_matrix = c()

    #actually doing it
    #for(sensitivity_variable in names(kConstants)) {
    print(variables_to_exclude)
    for(sensitivity_index in 1:length(kConstants)) {
        sensitivity_variable = names(kConstants)[sensitivity_index]
        #cat(sensitivity_index, ':', sensitivity_variable, '\n')
        real_multipliers = sapply(
            sensitivity_multipliers,
            function(m) get_real_multiplier(sensitivity_variable, m, kConstants)
        )
        
        greatest_difference_all_5 = 0
        greatest_difference_indices = rep(0, 5)

        for(i in 1:5) {
            greatest_difference = 0
            gd_j = NULL
        
            keys = sapply(
                real_multipliers,
                function(m) {
                    ifelse(m == 1,
                        paste0(i),
                        paste0(i, sensitivity_variable, '-', m)
                    )
                }
            )

            for(j in 1:max_j) {
                values = sapply(keys, function(key) dd[[j]][[key]][[outcome_name]])
                this_greatest_difference = max(sapply(
                    values,
                    function(value) {
                        max(0,
                            sapply(
                                values,
                                function(value_) {
                                    abs(log(value) - log(value_))
                                }
                            ),
                            na.rm = TRUE
                        )
                    }
                ))
                #print(this_greatest_difference)
                if(this_greatest_difference >= greatest_difference) {
                    greatest_difference = this_greatest_difference
                    gd_j = j
                    #print(gd_j)
                }
            }
            if(greatest_difference > greatest_difference_all_5) {
                greatest_difference_all_5 = greatest_difference
            }
            greatest_difference_indices[i] = gd_j


            values = sapply(keys, function(key) dd[[gd_j]][[key]][[outcome_name]])
            null_value = dd[[gd_j]][[paste0(i)]][[outcome_name]]
            if(greatest_negative_difference == greatest_positive_difference ||
               greatest_negative_difference == -Inf ||
               greatest_positive_difference == Inf) {
                browser()
                stop('FAILURE.')
            #} else if(sensitivity_variable %in% variables_to_exclude) {
            #    if(i == 1) {
            #        cat('SKIPPED:', sensitivity_index, ':', sensitivity_variable, '\n')
            #        plot.new()
            #    } else {
            #        cat('\tand again\n')
            #    }
            } else {
                #cat('PLOTTED:', sensitivity_index, ':', sensitivity_variable, '\n')
                #print('on')
                #print(greatest_negative_difference)
                #print('mid')
                #print(greatest_positive_difference)
                #print(values)
                #print(log(values))
                #print('off')
                if(i == 1) {
                    if(null_value == 0) {
                        log_differences = sign(values) * 10 * (greatest_positive_difference - greatest_negative_difference)
                    } else {
                        log_differences = log(values) - log(null_value)
                    }
                    plot(
                        real_multipliers,
                        #log(values) - log(null_value),
                        log_differences,
                        ylim = c(greatest_negative_difference, greatest_positive_difference),
                        xlim = c(0.5, 1.5),
                        xlab = paste0(sensitivity_variable, ' (multiplier)'),
                        ylab = ylab,
                        type = 'b',
                        lwd = 4
                    )
                } else {
                    points(
                        real_multipliers,
                        log(values) - log(null_value),
                        type = 'b',
                        col = colors[i],
                        lwd = 4
                    )
                }
            }
        }
        greatest_differences = c(greatest_differences, greatest_difference_all_5)
        greatest_difference_indices_matrix = rbind(greatest_difference_indices_matrix, greatest_difference_indices)
    }
    plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("bottom", row.names, lwd = 4, col = colors)
    dev.off()
    print('PLOT COMPLETE')
    list(gd = greatest_differences, gdim = greatest_difference_indices_matrix)
}


panelwise_interesting_sensitivity_fn = function(
    folder_name,
    unique_ids,
    is_baselines,
    community_transmissions,
    work_R0s,
    dormitory_R0s,
    E0,
    initial_recovereds,
    initial_V2s,
    n_sims,
    dd = NULL,
    summary_names,
    summary_fns,
    masks,
    output_per_shifts,
    hourly_wages,
    eConstants
) {
    folder_name; unique_ids; is_baselines; community_transmissions; work_R0s; dormitory_R0s; E0; initial_recovereds; initial_V2s; n_sims; dd = NULL; summary_names; summary_fns; masks; output_per_shifts; hourly_wages; eConstants
    max_j = length(unique_ids)
    sensitivity_multipliers = c(0.5, 1, 1.5)

    if(is.null(dd)) {
        dd = make_dd(max_j, sensitivity_multipliers, kConstants,
                     folder_name, unique_ids, is_baselines,
                     community_transmissions, work_R0s, dormitory_R0s, E0,
                     initial_recovereds, initial_V2s, n_sims,
                     summary_names, summary_fns, masks,
                     output_per_shifts, hourly_wages, eConstants)
    }
    #dd <<- dd
    l_si = make_paneled_plot('figures-2024-04-16/shared-summary-sensitivity-plots-si.png',
                             'symptomatic_infections',
                             'Symptomatic infections (multiplier)', dd,
                             kConstants, sensitivity_multipliers, max_j,
                             'figures-2024-04-16/shared-sensitivity-symptomatic-infections.csv', ######
                             unique_ids) ######
    l_su = make_paneled_plot('figures-2024-04-16/shared-summary-sensitivity-plots-su.png',
                             'shifts_unavailable', 
                             'Shifts unavailable (multiplier)', dd,
                             kConstants, sensitivity_multipliers, max_j,
                             'figures-2024-04-16/shared-sensitivity-shifts-unavailable.csv', ######
                             unique_ids) ######
#    l_tc = make_one_parameter_paneled_plots('v16-summary-sensitivity-plots-tc.png',
#                             'total_cost',
#                             'Total cost (multiplier)', dd,
#                             kConstants, sensitivity_multipliers, max_j,
#                             'v16-sensitivity-total-cost.csv', ######
#                             unique_ids) ######
    l_tc = make_paneled_plot('figures-2024-04-16/summary-sensitivity-plots-tc.png',
                             'total_cost',
                             'Total cost (multiplier)', dd,
                             kConstants, sensitivity_multipliers, max_j,
                            'figures-2024-04-16/sensitivity-total-cost.csv', ######
                             unique_ids) ######


    #list(gd_tc = l_tc$gd, gdim_tc = l_tc$gdim, dd = dd)
    #list(gd_si = l_si$gd, gd_su = l_su$gd, gd_tc = l_tc$gd, gdim_si = l_si$gdim, gdim_su = l_su$gdim, gdim_tc = l_tc$gdim, dd = dd)
    list(gd_si = l_si$gd, gd_su = l_su$gd, gd_tc = l_tc$gd, gdim_si = l_si$gdim, gdim_su = l_su$gdim, gdim_tc = l_tc$gdim, dd = dd)
}

panelwise_r_eff_sensitivity_fn = function(
    folder_name,
    unique_ids,
    is_baselines,
    community_transmissions,
    work_R0s,
    dormitory_R0s,
    E0,
    initial_recovereds,
    initial_V2s,
    n_sims,
    dd = NULL,
    summary_names,
    summary_fns,
    masks,
    output_per_shifts,
    hourly_wages,
    eConstants
) {
    folder_name; unique_ids; is_baselines; community_transmissions; work_R0s; dormitory_R0s; E0; initial_recovereds; initial_V2s; n_sims; dd = NULL; summary_names; summary_fns; masks; output_per_shifts; hourly_wages; eConstants
    max_j = length(unique_ids)
    sensitivity_multipliers = c(0.5, 1, 1.5)

    if(is.null(dd)) {
        dd = make_dd(max_j, sensitivity_multipliers, kConstants,
                     folder_name, unique_ids, is_baselines,
                     community_transmissions, work_R0s, dormitory_R0s, E0,
                     initial_recovereds, initial_V2s, n_sims,
                     summary_names, summary_fns, masks,
                     output_per_shifts, hourly_wages, eConstants)
    }
    #dd <<- dd
    l_r_eff = make_paneled_plot('figures-2024-04-16/shared-summary-sensitivity-plots-r_eff--90-days.png',
                             'r_eff',
                             'Effective reproduction number (internal)', dd,
                             kConstants, sensitivity_multipliers, max_j,
                             'figures-2024-04-16/shared-sensitivity-r_eff--90-days.csv', ######
                             unique_ids) ######

    #list(gd_tc = l_tc$gd, gdim_tc = l_tc$gdim, dd = dd)
    #list(gd_si = l_si$gd, gd_su = l_su$gd, gd_tc = l_tc$gd, gdim_si = l_si$gdim, gdim_su = l_su$gdim, gdim_tc = l_tc$gdim, dd = dd)
    #list(gd_si = l_si$gd, gd_su = l_su$gd, gdim_si = l_si$gdim, gdim_su = l_su$gdim, dd = dd)
    list(gd = l_r_eff$gd, gdim = l_r_eff$gdim, dd = dd)
}


get_real_economic_multiplier = function(sensitivity_variable,
                               theoretical_multiplier,
                               eConstants) {
    sensitivity_variable; theoretical_multiplier; kConstants
    eConstants_ = eConstants
    eConstants_[[sensitivity_variable]] = theoretical_multiplier * eConstants_[[sensitivity_variable]]
    if(sensitivity_variable == 'threshold' && eConstants_[['threshold']] > 1) {
        eConstants_[['threshold']] = 1
    }
    #browser()
    list(get(sensitivity_variable, eConstants_) / get(sensitivity_variable, eConstants), eConstants_)
    #ccl = check_economic_consistency(eConstants_, altered_single_parameter = sensitivity_variable)
    #eConstants_fixed = get('fixed_constants', ccl)
    #if(!get('consistent', ccl) && !get('fixed', ccl)) {
    #    stop('Unfixable constants')
    #}
    #get(sensitivity_variable, eConstants_fixed) / get(sensitivity_variable, eConstants)
}


make_paneled_economic_plot = function(filename, outcome_name, ylab, dd, eConstants,
                             sensitivity_multipliers, max_j,
                             csv_filename, unique_ids) {#, selection_mode) {  ######
    filename; outcome_name; ylab; dd; eConstants; sensitivity_multipliers; max_j#; selection_mode
    png(filename, height = 200*5, width = 200*7)
    layout(matrix(c(1:29, 34, 30:34), ncol = 7))

    #figuring out bounds:
    greatest_positive_difference = 0
    greatest_negative_difference = 0
    variables_to_exclude = list()
    values_df = data.frame(parameter_set = character(), sensitivity_variable = character(), multiplier = numeric(), intervention = character(), value = numeric()) ######
    for(sensitivity_variable in names(eConstants)) {
        real_multipliers = sapply(
            sensitivity_multipliers,
            function(m) get_real_economic_multiplier(sensitivity_variable, m, eConstants)[[1]]
        )

        for(i in 1:5) {
            intervention = row.names[i] ######
            null_value = dd[[1]][[paste0(i)]][[outcome_name]]
            keys = sapply(
                real_multipliers,
                function(m) {
                    ifelse(m == 1,
                        paste0(i),
                        paste0(i, sensitivity_variable, '-', m)
                    )
                }
            )
            for(j in 1:max_j) {
                values = sapply(keys, function(key) dd[[j]][[key]][[outcome_name]])
                for(multiplier in real_multipliers) { ######
                    if(multiplier == 1) { ######
                        key = paste0(i) ######
                    } else { ######
                        key = paste0(i, sensitivity_variable, '-', multiplier) ######
                    } ######
                    value = dd[[j]][[key]][[outcome_name]] ######
                    values_df = rbind(values_df, data.frame(parameter_set = unique_ids[j], sensitivity_variable = sensitivity_variable, multiplier = multiplier, intervention = intervention, value = value)) ######
                } ######
                this_greatest_positive_difference = max(0, log(values[3] / values[2]), log(values[1] / values[2]), na.rm = TRUE)
                this_greatest_negative_difference = min(0, log(values[3] / values[2]), log(values[1] / values[2]), na.rm = TRUE)
                if(this_greatest_positive_difference >= greatest_positive_difference) {
                    if(this_greatest_positive_difference == Inf) {
                        cat('\nPositive Infinite:\n', sensitivity_variable, '\n', i, '\n', j, '\n', values, '\n\n')
                        variables_to_exclude = c(variables_to_exclude, sensitivity_variable)
                    } else {
                        greatest_positive_difference = this_greatest_positive_difference
                    }
                }
                if(this_greatest_negative_difference <= greatest_negative_difference) {
                    if(this_greatest_negative_difference == Inf) {
                        cat('\nNegative Infinite:\n', sensitivity_variable, '\n', i, '\n', j, '\n', values, '\n\n')
                        variables_to_exclude = c(variables_to_exclude, sensitivity_variable)
                    } else {
                        greatest_negative_difference = this_greatest_negative_difference
                    }
                }
            }
        }
    }
    write.csv(values_df, csv_filename) ######

    greatest_differences = c()
    greatest_difference_indices_matrix = c()

    #actually doing it
    print(variables_to_exclude)
    for(sensitivity_index in 1:length(eConstants)) {
        sensitivity_variable = names(eConstants)[sensitivity_index]
        real_multipliers = sapply(
            sensitivity_multipliers,
            function(m) get_real_economic_multiplier(sensitivity_variable, m, eConstants)[[1]]
        )
        
        greatest_difference_all_5 = 0
        greatest_difference_indices = rep(0, 5)

        for(i in 1:5) {
            greatest_difference = 0
            gd_j = NULL
        
            keys = sapply(
                real_multipliers,
                function(m) {
                    ifelse(m == 1,
                        paste0(i),
                        paste0(i, sensitivity_variable, '-', m)
                    )
                }
            )

            for(j in 1:max_j) {
                values = sapply(keys, function(key) dd[[j]][[key]][[outcome_name]])
                #cat(i, ':', j, '\n')
                #browser()
                this_greatest_difference = max(sapply(
                    values,
                    function(value) {
                        max(0,
                            sapply(
                                values,
                                function(value_) {
                                    abs(log(value) - log(value_))
       sensitivity-2022-11-22                         }
                            ),
                            na.rm = TRUE
                        )
                    }
                ))
                if(this_greatest_difference >= greatest_difference) {
                    greatest_difference = this_greatest_difference
                    gd_j = j
                }
            }
            if(greatest_difference > greatest_difference_all_5) {
                greatest_difference_all_5 = greatest_difference
            }
            greatest_difference_indices[i] = gd_j


            values = sapply(keys, function(key) dd[[gd_j]][[key]][[outcome_name]])
            null_value = dd[[gd_j]][[paste0(i)]][[outcome_name]]
            if(greatest_negative_difference == greatest_positive_difference ||
               greatest_negative_difference == -Inf ||
               greatest_positive_difference == Inf) {
                stop('FAILURE.')
            } else {
                if(i == 1) {
                    if(null_value == 0) {
                        log_differences = sign(values) * 10 * (greatest_positive_difference - greatest_negative_difference)
                    } else {
                        log_differences = log(values) - log(null_value)
                    }
                    plot(
                        real_multipliers,
                        log_differences,
                        ylim = c(greatest_negative_difference, greatest_positive_difference),
                        xlim = c(0.5, 1.5),
                        xlab = paste0(sensitivity_variable, ' (multiplier)'),
                        ylab = ylab,
                        type = 'b',
                        lwd = 4
                    )
                } else {
                    points(
                        real_multipliers,
                        log(values) - log(null_value),
                        type = 'b',
                        col = colors[i],
                        lwd = 4
                    )
                }
            }
        }
        greatest_differences = c(greatest_differences, greatest_difference_all_5)
        greatest_difference_indices_matrix = rbind(greatest_difference_indices_matrix, greatest_difference_indices)
    }
    plot(NULL, xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend("bottom", row.names, lwd = 4, col = colors)
    dev.off()
    print('PLOT COMPLETE')
    list(gd = greatest_differences, gdim = greatest_difference_indices_matrix)
}


make_economic_batch = function(d, keys, index_j,
    max_j,
    folder_name,
    unique_ids,
    is_baselines,
    community_transmissions,
    work_R0s,
    dormitory_R0s,
    E0,
    initial_recovereds,
    initial_V2s,
    n_sims,
    summary_names,
    summary_fns, 
    masks,
    output_per_shifts,
    hourly_wages,
    eConstants
) {
    d; keys; index_j; max_j; folder_name; unique_ids; is_baselines; community_transmissions; work_R0s; dormitory_R0s; E0; initial_recovereds; initial_V2s; n_sims; summary_names; summary_fns; masks; output_per_shifts; hourly_wages; eConstants
    #prepended_key = ifelse(key == '',
    #    '',
    #    paste0('-', key)
    #)
    for(i in 1:5) {
        filename = rds_filename('', i, index_j, max_j, folder_name,
                                unique_ids, is_baselines,
                                community_transmissions, work_R0s,
                                dormitory_R0s, E0, initial_recovereds,
                                initial_V2s, n_sims)
        data_ = readRDS(filename)

        for(key in keys) {
            if(length(key) == 1 && key == '') {
                key_ = paste0(i)
                eConstants_ = eConstants
            } else {
                key_ = paste0(i, key[[1]], '-', key[[2]])
                eConstants_ = key[[3]]
                #browser()
            }
            
            d[[key_]] = list()
            #browser()
            for(k in 1:length(summary_names)) {
                d[[key_]][[summary_names[k]]] = summary_fns[[k]](
                    data_, masks[[index_j]], i, output_per_shifts[[index_j]],
                    hourly_wages[[index_j]], eConstants_#key[[3]]
                )
            }
        }
    }
    d
}

make_economic_dd = function(
    max_j,
    sensitivity_multipliers,
    kConstants,
    folder_name,
    unique_ids,
    is_baselines,
    community_transmissions,
    work_R0s,
    dormitory_R0s,
    E0,
    initial_recovereds,
    initial_V2s,
    n_sims,
    summary_names,
    summary_fns,
    masks,
    output_per_shifts, hourly_wages, eConstants
) {
    max_j; sensitivity_multipliers; kConstants; folder_name; unique_ids; is_baselines; community_transmissions; work_R0s; dormitory_R0s; E0; initial_recovereds; initial_V2s; n_sims; summary_names; summary_fns; masks; output_per_shifts; hourly_wages; eConstants
    dd = list()

    keys = list('')
    for(sensitivity_variable in names(eConstants)) {
        real_multipliers = lapply(sensitivity_multipliers, function(m) get_real_economic_multiplier(sensitivity_variable, m, eConstants))
        for(sensitivity_multiplier in real_multipliers) {
            if(sensitivity_multiplier[[1]] != 1) {
                keys[[length(keys) + 1]] = list(sensitivity_variable, sensitivity_multiplier[[1]], sensitivity_multiplier[[2]])
            }
        }
    }   

    for(index_j in 1:max_j) {
        #cat('\n\n', index_j, ':', dormitory_R0s, '\n\n')
        dd[[index_j]] = make_economic_batch(list(), keys, index_j, max_j, folder_name,
                                   unique_ids, is_baselines,
                                   community_transmissions, work_R0s,
                                   dormitory_R0s, E0, initial_recovereds,
                                   initial_V2s, n_sims, summary_names, summary_fns,
                                   masks, output_per_shifts, hourly_wages,
                                   eConstants)
    }

"    for(sensitivity_variable in names(kConstants)) {
        real_multipliers = sapply(sensitivity_multipliers, function(m) get_real_multiplier(sensitivity_variable, m, kConstants))
        for(sensitivity_multiplier in real_multipliers) {
            if(sensitivity_multiplier != 1) {
                key = paste0(sensitivity_variable, '-', sensitivity_multiplier)
                for(index_j in 1:max_j) {
                    dd[[index_j]] = make_batch(dd[[index_j]], key, index_j,
                                               max_j, folder_name, unique_ids,
                                               is_baselines,
                                               community_transmissions,
                                               work_R0s, dormitory_R0s, E0,
                                               initial_recovereds, initial_V2s,
                                               n_sims, summary_names,
                                               summary_fns, masks,
                                               output_per_shifts, hourly_wages,
                                               eConstants)
                }
            }
        }
    }"
    dd
}


pi_economic_sensitivity_fn = function(
    folder_name,
    unique_ids,
    is_baselines,
    community_transmissions,
    work_R0s,
    dormitory_R0s,
    E0,
    initial_recovereds,
    initial_V2s,
    n_sims,
    dd = NULL,
    summary_names,
    summary_fns,
    masks,
    output_per_shifts,
    hourly_wages,
    eConstants
) {
    folder_name; unique_ids; is_baselines; community_transmissions; work_R0s; dormitory_R0s; E0; initial_recovereds; initial_V2s; n_sims; dd = NULL; summary_names; summary_fns; masks; output_per_shifts; hourly_wages; eConstants
    max_j = length(unique_ids)
    sensitivity_multipliers = c(0.5, 1, 1.5)

    if(is.null(dd)) {
        dd = make_economic_dd(max_j, sensitivity_multipliers, kConstants,
                     folder_name, unique_ids, is_baselines,
                     community_transmissions, work_R0s, dormitory_R0s, E0,
                     initial_recovereds, initial_V2s, n_sims,
                     summary_names, summary_fns, masks,
                     output_per_shifts, hourly_wages, eConstants)
    }
    #dd <<- dd
    #return here
    l_tc = make_paneled_economic_plot('figures-2024-04-16/summary-sensitivity-plots-economic-tc.png',
                             'total_cost',
                             'Total cost (multiplier)', dd,
                             eConstants, sensitivity_multipliers, max_j,
                             'figures-2024-04-16/sensitivity-economic-only-parameters-total-cost.csv', ######
                             unique_ids) ######

    list(gd_tc = l_tc$gd, gdim_tc = l_tc$gdim, dd = dd)
}


#dd = readRDS('saved_dd.RDS')

#facility_production_mask = c(sapply(0:13, function(y) 21 * y + c(sapply(0:4, function(x) 3*x + 1:2))))[1:130]
#farm_production_mask = c(sapply(0:13, function(y) 21 * y + c(sapply(0:4, function(x) 3*x + 1))))[1:65]

workday = c('work', 'work', 'work')
day_off = c('home', 'home', 'sleep')
#week = c(rep(workday, 5), rep(day_off, 2))
days = 90

work_shifts = function(start_day) {
    #print(start_day)
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
    #production_shifts = work_shifts = (schedule == 'work')
    (schedule == 'work')
}

production_shifts = function(start_day) {
    ws = work_shifts(start_day)
    l = length(ws)
    (ws & ((1:l) %% 3 != 0))
}

production_shifts_mask_fn = function(start_days) {
    sapply(1:100, function(x) production_shifts(start_days[x]))
}

#start_days = readRDS(paste0(fragments[1], '/start_days--', fragments[2]))
start_days = readRDS('H_R_V2-check--sensitivity/start_days--farmlike-facility_community-0,work_R0-6,dormitory_R0-2,E0-1,initial_recovered-71,initial_V2-73,n_sims-100index_i-1_full-output.rds')

masks = list(production_shifts_mask_fn(start_days)) #note that this is now a list of 270 x 100 matrices, not 
#masks = list(facility_production_mask)#rep(list(farm_production_mask, facility_production_mask), 8)

mean_fn = function(fn) {
    function(x, mask, limited_i, output_per_shift, hourly_wage, eConstants) {
        mean(fn(x, mask, limited_i, output_per_shift, hourly_wage, eConstants))
    }
}

output_per_shifts = rep(784346.67 / 10, 4) 
hourly_wages = rep(16.57, 4)
#output_per_shifts = rep(1680000 / 10, 16)
#hourly_wages = rep(13.89, 16)
#eConstants

l = panelwise_interesting_sensitivity_fn(
    'H_R_V2-check--sensitivity',
    c('farmlike-facility'#,
      # 'facility-no-vax',
      # 'facility-no-recovered',
      # 'facility-start-of-epidemic'
      ),
    c(FALSE#, FALSE, FALSE, FALSE
      ),
    c(0#,0.002,0.002,0.002
      ),
    c(6#, 6, 6, 6
      ),
    c(2#,0,0,0
      ),
    1,
    c(71#, 71, 0, 0
      ),
    c(73#, 0, 73, 0
      ),
    100,
    dd = NULL,
    c('symptomatic_infections', 'shifts_unavailable', 'total_cost'),
    c(mean_fn(symptomatic_infections),
      mean_fn(shiftwise_unavailable),
      mean_fn(g)),
    masks,
    output_per_shifts, hourly_wages, eConstants
)
#v = pmax(l$gd_si, l$gd_su, l$gd_tc)
#print(sort(v))
#cutoff = sort(v)[15]
#print(names(kConstants)[v >= cutoff])

dd = l$dd
dir.create('figures-2024-04-16/')
saveRDS(dd, 'figures-2024-04-16/saved_dd.RDS')


discordant = which(
    names(dd[[1]])[6:485] != c(
        sapply(
            names(kConstants),
            function(s) sapply(
                c(0.5, 1.5),
                function(m) sapply(
                    1:5, function(i) paste0(i,s,'-',m)
                )
            )
        )
    )
)
discordant
# [1] 166 167 168 169 170 321 322 323 324 325
names(dd[[1]])[discordant + 5]
# [1] "1boosting_interval-1.39802631578947"   
# [2] "2boosting_interval-1.39802631578947"   
# [3] "3boosting_interval-1.39802631578947"   
# [4] "4boosting_interval-1.39802631578947"   
# [5] "5boosting_interval-1.39802631578947"   
# [6] "1time_since_first_V2-0.715962441314554"
# [7] "2time_since_first_V2-0.715962441314554"
# [8] "3time_since_first_V2-0.715962441314554"
# [9] "4time_since_first_V2-0.715962441314554"
#[10] "5time_since_first_V2-0.715962441314554"

generate_name = function(i, s, m) {
    if(s == 'boosting_interval' && m == 1.5) {
        m = 1.39802631578947
    } else if(s == 'time_since_first_V2' && m == 0.5) {
        m = 0.715962441314554
    }
    paste0(i,s,'-',m)
}

difference = function(dd, i, ss, outcome) {
    d = exp(abs(log(dd[[1]][[ss]][[outcome]] / dd[[1]][[i]][[outcome]])))
}

max_difference_ = function(dd, i, ss) {
    max(
        sapply(
            c(
                'symptomatic_infections',
                'shifts_unavailable',
                'total_cost'
            ),
            function(outcome) {
                difference(dd, i, ss, outcome)
            }
        ),
        na.rm = TRUE
    )
}

max_difference = function(dd, s) {
    max(
        sapply(
            c(0.5, 1.5),
            function(m) sapply(
                1:5, function(i) max_difference_(dd, i, generate_name(i,s,m))
            )
        )
    )
}

#differences = sapply(names(kConstants), function(s) max_difference(dd, s))
#o = order(differences)
#cat(sapply(1:length(o), function(i) paste0(names(kConstants)[o][i], ':\t', differences[o][i])), sep = '\n')
#...
#duration_IM_mean:       2.02835820895522
#duration_IM_shape:      2.05597579425114
#mu:     2.06159169550173
#R_question_period:      2.06818181818182
#V2_decay_rate:  2.09237875288684
#duration_IA_mean:       2.16666666666667
#H_RB_nsp_a:     2.18840579710145
#duration_IP_mean:       2.22479462285288
#SEVERE_MULTIPLIER:      2.38212094653812

#but wait. This doesn't take account of 0.5 vs. 1.5 

#stop('Not using r_eff right now.')

generate_name = function(i, s, m) {
    if(s == 'boosting_interval' && m == 1.5) {
        m = 1.39802631578947
    } else if(s == 'time_since_first_V2' && m == 0.5) {
        m = 0.715962441314554
    }
    paste0(i,s,'-',m)
}

#difference = function(dd, i, ss, outcome) {
#    d = exp(abs(log(dd[[1]][[ss]][[outcome]] / dd[[1]][[i]][[outcome]])))
#}

max_difference_ = function(dd, s, i, outcome) {
    a = dd[[1]][[generate_name(i,s,0.5)]][[outcome]]
    b = dd[[1]][[i]][[outcome]]
    c = dd[[1]][[generate_name(i,s,1.5)]][[outcome]]
    max(exp(abs(log(c(a/b, b/c, a/c)))),
        na.rm = TRUE
    )
}

max_difference = function(dd, s) {
    max(
        sapply(
            c(
                'symptomatic_infections',
                'shifts_unavailable',
                'total_cost'
            ),
            function(outcome) sapply(1:5, function(i) max_difference_(dd, s, i, outcome))
        ),
        na.rm = TRUE
    )
}

differences = sapply(names(kConstants), function(s) max_difference(dd, s))
o = order(differences)
cat(sapply(1:length(o), function(i) paste0(names(kConstants)[o][i], ':\t', differences[o][i])), sep = '\n')
#...
#V2_magnitude:   1.80345710627401
#p_trans_IM:     1.94074074074074
#R_nsp_a:        1.97115384615385
#-------------------------------------------------
#mu:     2.06159169550173
#R_question_period:      2.06818181818182
#duration_IM_shape:      2.17391304347826
#H_RB_nsp_a:     2.18840579710145
#V2_decay_rate:  2.18937644341801
#duration_IP_mean:       2.22479462285288
#duration_IA_mean:       2.2595399188092
#fraction_ssp_symptomatic:       2.47796143250689
#duration_IS_mean:       2.50054171180932
#SEVERE_MULTIPLIER:      2.55302366345311
#duration_IM_mean:       2.63333333333333

#values of sensitivities that were > 2 before bobrovitz, but aren't now

#for comparison with older results
"dd_old = readRDS('saved_dd_17-shared.RDS')
old_differences = sapply(names(kConstants), function(s) max_difference(dd_old, s))
oo = order(old_differences)
cat(sapply(1:length(oo), function(i) paste0(names(kConstants)[oo][i], ':\t', old_differences[oo][i])), sep = '\n')
#...
#V2_decay_rate:  1.16885553470919                   *
#sd:     1.39328063241107
#boosting_interval:      1.43953813790413
#p_trans_IA:     1.45474137931034
#duration_IM_shape:      1.46575457951524           *
#mu:     1.64516129032258                           *
#B_decay_rate_1: 1.66149870801034
#duration_IP_shape:      1.69432986828613
#duration_IA_mean:       1.78862431769638           *
#p_trans_IP:     1.80821917808219
#complete_immunity_duration_R:   1.90909090909091

#surge in V2_decay_rate is weird, kinda concerning

cat(sapply(1:5, function(i) paste0(dd[[1]][[i]][['symptomatic_infections']], '\t', dd[[1]][[generate_name(i, 'V2_decay_rate', 0.5)]][['symptomatic_infections']], '\t', dd[[1]][[generate_name(i, 'V2_decay_rate', 1.5)]][['symptomatic_infections']])), sep='\n')
#16.83   15.53   17.04
#15.39   14.34   15.59
# 1.82    1.01    1.92                                          #1.90099
#11.99    8.57   12.36
#14.78   12.69   14.97

cat(sapply(1:5, function(i) paste0(dd_old[[1]][[i]][['symptomatic_infections']], '\t', dd_old[[1]][[generate_name(i, 'V2_decay_rate', 0.5)]][['symptomatic_infections']], '\t', dd_old[[1]][[generate_name(i, 'V2_decay_rate', 1.5)]][['symptomatic_infections']])), sep='\n')
#35.83   35.15   35.82
#34.79   34.02   34.82
#5.58    5.33    6.23                                           #1.168856 virus testing at p = 0.3
#30.28   29.34   30.89
#33.52   32.77   33.66

cat(sapply(1:5, function(i) paste0(dd[[1]][[i]][['shifts_unavailable']], '\t', dd[[1]][[generate_name(i, 'V2_decay_rate', 0.5)]][['shifts_unavailable']], '\t', dd[[1]][[generate_name(i, 'V2_decay_rate', 1.5)]][['shifts_unavailable']])), sep='\n')
#11.47              6.96333333333333        11.52               #1.65438
#27.4066666666667   22.3633333333333        27.9533333333333    
#25.2733333333333   17.39                   26.2833333333333
# 9.06               4.33                    9.48               #2.189376 in 40% R0 reduction
# 9.93               4.99333333333333        9.78
#So _reducing_ V2_decay_rate can have a big impact on unavailability, especially in the absence of testing. This might be chance? What does the equivalent look like for dd_old? (It also could be a systemic effect, due to r_eff being kinda marginal to begin with . . . hard to say.
cat(sapply(1:5, function(i) paste0(dd_old[[1]][[i]][['shifts_unavailable']], '\t', dd_old[[1]][[generate_name(i, 'V2_decay_rate', 0.5)]][['shifts_unavailable']], '\t', dd_old[[1]][[generate_name(i, 'V2_decay_rate', 1.5)]][['shifts_unavailable']])), sep='\n')
#29.6833333333333   28.5366666666667        30.2
#70.65              67.6966666666667        70.7833333333333
#59.48              57.42                   64.9966666666667    #1.131952
#25.94              25.0966666666667        25.8366666666667    #1.033603
#28.8466666666667   27.3433333333333        28.97"

#Let's take a moment to make sure we haven't somehow greatly increased the number of V2 . . .

#okay, so nothing surprising here . . .

l_r_eff = panelwise_r_eff_sensitivity_fn(
    'H_R_V2-check--sensitivity-r0s--90-days',
    c('farmlike-facility'#,
      # 'facility-no-vax',
      # 'facility-no-recovered',
      # 'facility-start-of-epidemic'
      ),
    c(FALSE#, FALSE, FALSE, FALSE
      ),
    c(0#,0.002,0.002,0.002
      ),
    c(6#, 6, 6, 6
      ),
    c(2#,0,0,0
      ),
    1,
    c(71#, 71, 0, 0
      ),
    c(73#, 0, 73, 0
      ),
    100,
    dd = NULL,
    c('r_eff'),
    c(mean_fn(r_eff)),
    masks,
    output_per_shifts, hourly_wages, eConstants
)
#v = pmax(l$gd_si, l$gd_su, l$gd_tc)
#print(sort(v))
#cutoff = sort(v)[15]
#print(names(kConstants)[v >= cutoff])

dd_l = l_r_eff$dd
saveRDS(dd_l, 'figures-2024-04-16/saved_dd-r_eff--90-days.RDS')

"cat(sapply(1:5, function(i) paste0(dd_l[[1]][[i]][['r_eff']], '\t', dd_l[[1]][[generate_name(i, 'V2_decay_rate', 0.5)]][['r_eff']], '\t', dd_l[[1]][[generate_name(i, 'V2_decay_rate', 1.5)]][['r_eff']])), sep='\n')
#2.75    1.88    2.59
#2.43    1.73    2.3
#0.94    0.81    0.95
#2.13    1.45    2.01   #1.468966
#2.71    1.87    2.55

#and to compare ...
dd_l_old = readRDS('saved_dd_17-shared-r_eff.RDS')
cat(sapply(1:5, function(i) paste0(dd_l_old[[1]][[i]][['r_eff']], '\t', dd_l_old[[1]][[generate_name(i, 'V2_decay_rate', 0.5)]][['r_eff']], '\t', dd_l_old[[1]][[generate_name(i, 'V2_decay_rate', 1.5)]][['r_eff']])), sep='\n')"
#4.48    4.03    4.55
#4.18    3.73    4.25
#1.35    1.2     1.37
#3.35    3.03    3.37
#4.47    4.03    4.54

#Okay, broadly speaking, this seems reasonable. Tentatively. We're just seeing a more marginal system be more likely to tip over the edge.#Here's a thought: What is the general trend comparing maxima for the same variable?
#Pause a moment and _notice_: Even at baseline, R_eff is < 1 for p = 0.3 testing! This is an important change that likely explains some of what we are seeing.

max_difference_r_eff = function(dd, s) {
    max(sapply(1:5, function(i) max_difference_(dd, s, i, 'r_eff')))
}

r_eff_differences = sapply(names(kConstants), function(s) max_difference_r_eff(dd_l, s))

combined_differences = pmax(differences, r_eff_differences)
o_combined = order(combined_differences)

"r_eff_old_differences = sapply(names(kConstants), function(s) max_difference_r_eff(dd_l_old, s))
old_combined_differences = pmax(old_differences, r_eff_old_differences)
oo_combined = order(old_combined_differences, combined_differences)
cat(sapply(1:length(oo_combined), function(i) paste0(names(kConstants)[oo_combined][i], ':\t', old_combined_differences[oo_combined][i], '\t', combined_differences[oo_combined][i])), sep = '\n')"
#somewhat confusing results

print('And combined:')

#now, practical question: What are our top sensitivity variables?
cat(sapply(1:length(o_combined), function(i) paste0(names(kConstants)[o_combined][i], ':\t', combined_differences[o_combined][i])), sep = '\n')
#below is old; now we have:
#fraction_ssp_symptomatic:       2.21001757469244
#p_trans_IM:     2.47552447552448
#duration_IS_mean:       2.90022675736961
#duration_IM_mean:       3.09848484848485
#SEVERE_MULTIPLIER:      3.36772486772487


#mu:     2.06159169550173
#R_question_period:      2.06818181818182
#duration_IM_shape:      2.17391304347826
#H_RB_nsp_a:     2.18840579710145
#V2_decay_rate:  2.18937644341801
#duration_IP_mean:       2.22479462285288
#duration_IA_mean:       2.2595399188092
#fraction_ssp_symptomatic:       2.47796143250689
#duration_IS_mean:       2.50054171180932
#SEVERE_MULTIPLIER:      2.55302366345311
#duration_IM_mean:       2.63333333333333

#Reordering for presentation

#duration_IM_mean:       2.63333333333333
#duration_IS_mean:       2.50054171180932
#duration_IA_mean:       2.2595399188092
#duration_IP_mean:       2.22479462285288
#mu:                     2.06159169550173        #this is essentially a proxy for duration_E_mean, probably
                                                #is this a case of longer -> higher probability of detection?
                                                #no, because we don't detect then. So why?
                                                #specifically x-1.5 gives a reduced fraction of any symptomatic under p = 0.3 testing. why?
#duration_IM_shape:      2.17391304347826

#SEVERE_MULTIPLIER:      2.55302366345311

#fraction_ssp_symptomatic:       2.47796143250689
#V2_decay_rate:  2.18937644341801
#H_RB_nsp_a:     2.18840579710145
#R_question_period:      2.06818181818182


#stop('Good enough for the moment.')

#Bits to do by hand, probably
#this bit probably actually linear

#custom_linear_panel('r_eff', 'Effective Reproduction Number (internal)', l_r_eff$dd, kConstants, c(0.5, 1, 1.5), 'farmlike-facility', 'B_magnitude_1', c(0,10), 1)

custom_linear_panel = function(outcome_name, ylab, dd, kConstants,
                             sensitivity_multipliers, unique_id,
                             sensitivity_variable, ylim, gd_j = 1, log_scale = FALSE) {#, selection_mode) {  ######
    #print(outcome_name)
    outcome_name; ylab; dd; kConstants; sensitivity_multipliers; unique_id

    #actually doing it
    real_multipliers = sapply(
        sensitivity_multipliers,
        function(m) get_real_multiplier(sensitivity_variable, m, kConstants)
    )
    for(i in 1:5) {        
        keys = sapply(
            real_multipliers,
            function(m) {
                ifelse(m == 1,
                    paste0(i),
                    paste0(i, sensitivity_variable, '-', m)
                )
            }
        )
        values = sapply(keys, function(key) dd[[gd_j]][[key]][[outcome_name]])
        if(log_scale) {
            values = log(values + 1)
        }
        #null_value = dd[[gd_j]][[paste0(i)]][[outcome_name]]
        #browser()
        if(i == 1) {
            plot(
                real_multipliers,
                values,
                ylim = ylim,
                xlim = c(0.5, 1.5),
                xlab = paste0(sensitivity_variable, ' (multiplier)'),
                ylab = ylab,
                type = 'b',
                lwd = 4
            )
        } else {
            points(
                real_multipliers,
                values,
                type = 'b',
                col = colors[i],
                lwd = 4
            )
        }
    }
}

png('figures-2024-04-16/master-summary--90-days.png', height = 200*5, width = 200*4)
#layout(matrix(1:45, ncol = 5, byrow = TRUE))
#layout(matrix(1:44, ncol = 4, byrow = TRUE))
layout(matrix(1:20, ncol = 4, byrow = TRUE))
parameters = c('SEVERE_MULTIPLIER', 'duration_IM_mean', 'duration_IS_mean', 'p_trans_IM',  'fraction_ssp_symptomatic')
#parameters = c('duration_IM_mean', 'duration_IS_mean', 'duration_IA_mean', 'duration_IP_mean', 'mu', 'duration_IM_shape', 'SEVERE_MULTIPLIER', 'fraction_ssp_symptomatic', 'V2_decay_rate', 'H_RB_nsp_a', 'R_question_period')
#parameters = c('SEVERE_MULTIPLIER', 'duration_IM_mean', 'p_trans_IM', 'B_magnitude_2', 'duration_IS_mean', 'isolation_duration', 'R_question_period', 'duration_IP_mean', 'B_magnitude_1')
#parameters = c('SEVERE_MULTIPLIER', 'duration_IM_mean', 'duration_IS_mean', 'duration_IP_mean', 'p_trans_IM', 'B_magnitude_2', 'B_magnitude_1', 'isolation_duration', 'R_question_period')

for(i in 1:5) {
    parameter = parameters[i]
    custom_linear_panel('symptomatic_infections', 'Symptomatic Infections', l$dd, kConstants, c(0.5, 1, 1.5), 'farmlike-facility', parameter, c(0,30), 1)
    custom_linear_panel('shifts_unavailable', 'Worker-Shifts Unavailable', l$dd, kConstants, c(0.5, 1, 1.5), 'farmlike-facility', parameter, c(0,60), 1)
    custom_linear_panel('r_eff', 'Effective Reproduction Number (internal)', l_r_eff$dd, kConstants, c(0.5, 1, 1.5), 'farmlike-facility', parameter, c(0,4), 1)
    custom_linear_panel('total_cost', 'Total Cost ($)', l$dd, kConstants, c(0.5, 1, 1.5), 'farmlike-facility', parameter, c(0,40000), 1) #was 51000
    #custom_linear_panel('total_cost', 'log(Total Cost) (log($))', l$dd, kConstants, c(0.5, 1, 1.5), 'farmlike-facility', parameter, c(0,11), 1, TRUE) #was 51000
}
dev.off()

stop('Return to normal running?')

l2 = pi_economic_sensitivity_fn(
    'sensitivity-2022-11-22',
    c('farmlike-facility',
      'farmlike-facility-start-of-epidemic',
      'farmlike-facility-no-vax',
      'farmlike-facility-no-recovered'),
    c(FALSE, FALSE, FALSE, FALSE),
    c(0, 0,
      0,
      0),
    c(6,
      6, 6, 6),
    c(2, 2, 2, 2),
    1,
    c(71,
      0, 
      71,0),
    c(73, 0, 0,
      73,),
    100,
    dd = NULL,
    c('total_cost'),
    c(mean_fn(g)),
    masks,
    output_per_shifts, hourly_wages, eConstants
)

dd2 = l2$dd
saveRDS(dd2, 'saved_dd2_14.RDS')

# [1] 0.02943518 0.03385564 0.05364339 0.05606389 0.07001799 0.07627487
# [7] 0.07748810 0.08334730 0.10040950 0.11874043 0.12788328 0.18079600
#[13] 0.19454219 0.20439436 0.26785635 0.30954087 0.31035215 0.31138418
#[19] 0.33153395 0.34520592 0.35681882 0.38269800 0.41976294 0.46789730
#[25] 0.49411372 0.59414166 0.63206588 0.65465846 0.77133011 0.94025349
#[31] 0.94748234 1.20980732 1.22516004
# [1] "isolation_duration"           "mu"
# [3] "duration_IP_mean"             "duration_IP_shape"
# [5] "duration_IA_mean"             "duration_IM_mean"
# [7] "duration_IS_mean"             "p_trans_IP"
# [9] "p_trans_IA"                   "p_trans_IM"
#[11] "boosting_interval"            "complete_immunity_duration_R"
#[13] "B_magnitude_1"                "B_magnitude_2"
#[15] "B_decay_rate_1"               "B_decay_rate_2"
#[17] "SEVERE_MULTIPLIER"            "R_question_period"


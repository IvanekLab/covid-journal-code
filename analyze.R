# analyze.R is part of Food INdustry CoViD Control Tool
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

#limited_runs_index = c(1,2,4,9,13)
output_per_shift = output_per_week / (5 * (1 + (supervisors > 1 && tolower(farm_or_facility) == 'facility')))

library('vioplot')

######## analyze model predictions
analyze_fn = function() {  #this may, in the future, be revised to provide
                           #better encapsulation; for now, it simply serves
                           #to provide more meaningful debugging data

ANALYZE = TRUE
source('double-wrapped.R', local = TRUE)
list_ = double_wrapped_fn()
row.names = list_[[1]]
colors = list_[[2]]
ltys = list_[[3]]
full_output_filenames = list_[[4]]

if(farm_or_facility == 'farm') {
    production_shift_size = N
    cleaning_shift_size = 0
    workday = c('work', 'home', 'home')
    day_off = c('home', 'home', 'home')
} else {

    psX_only_size = 1 + workers_per_crew * crews_per_supervisor + n_shift_floaters
    if(supervisors > 1) {
        on_ps_1 = c(1/3, rep(1, psX_only_size), rep(0, psX_only_size),
                    rep(0, n_cleaners), rep(1/3, n_all_floaters))
        on_ps_2 = c(1/3, rep(0, psX_only_size), rep(1, psX_only_size),
                    rep(0, n_cleaners), rep(1/3, n_all_floaters))
        on_cs = c(1/3, rep(0, 2 * psX_only_size), rep(1, n_cleaners),
                rep(1/3, n_all_floaters))
        workday = c('work', 'work', 'work')
    } else {
        on_ps_1 = c(1/2, rep(1, psX_only_size), rep(0, n_cleaners),
                    rep(1/2, n_all_floaters))
        on_ps_2 = rep(0, 1 + psX_only_size + n_cleaners + n_all_floaters)
        on_cs = c(1/2, rep(0, psX_only_size), rep(1, n_cleaners),
                rep(1/2, n_all_floaters))
        workday = c('work', 'home', 'work')
    } 
    day_off = c('home', 'home', 'sleep')
    week = c(rep(workday, 5), rep(day_off, 2))
    
    production_shift_size = sum(on_ps_1)
    cleaning_shift_size =  sum(on_cs)
}

work_shifts = function(start_day) {
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
    (schedule == 'work')
}

production_shifts = function(start_day) {
    ws = work_shifts(start_day)
    l = length(ws)
    (ws & ((1:l) %% 3 != 0))
}

cleaning_shifts = function(start_day) {
    ws = work_shifts(start_day)
    l = length(ws)
    (ws & ((1:l) %% 3 == 0))
}


#summary plots
combine = function(data, outcome_fn, summary_fn, mask, default_value = NA, default_detector = is.na) { 
    #mask has dimensions of times (e.g., 270) x runs (e.g., 100 or 1000)
    dimnames(data) = list(rep(NA, dim(data)[1]),
                          colnames(data),
                          rep(NA, dim(data)[3])
    )
    #above is a bit kludgey, but it works -- at some point in the future,
    #we may explicitly save data with dimnames
    outcomes = outcome_fn(data)
    summarize = function(i) { #time index
        if(sum(mask[i,]) == 0) {
            # Meaning: for this time i, there is no run j for which
            # mask[i,j] == TRUE
            # Realistically, this shouldn't happen solely by chance, due to
            # random start days, even at only 100 runs, let alone 1000. But it
            # can happen deterministically, e.g., for cleaning shifts when
            # calculating production plots, because we always start at the same
            # time of day, just not necessarily the same day of the week.
            default_value
        } else {
            summary_fn(outcomes[i,][mask[i,]])
        }
    }
    if(!identical(mask, NA)) {
        summarized = sapply(1:dim(outcomes)[1], summarize)
        sanity_check = sapply(1:dim(outcomes)[1], function(i) sum(mask[i,]) != 0)
        if(any(default_detector(summarized[sanity_check]))) {
            # This means that we are getting a default value from summarized,
            # for one or more times that *do* have runs for which that time is
            # included.
            browser() # For manual debugging; this may at some point be replaced with an error message
        }
    } else {
        summarized = apply(outcomes, 1, summary_fn)
    }

    summarized
}

#outcome_fn's
infected = function(data) {
    data[,'IA',] + data[,'IP',] + data[,'IM',] + data[,'IS',] + data[,'IC',]
}

symptomatic = function(data) {
    data[,'IM',] + data[,'IS',] + data[,'IC',]
}

day_add_all = function(ys) {
    len = length(ys)
    if(len %% 3) {
        stop('Incorrect ys length:', len)
    }
    index = (1:(len/3)) * 3
    ys[index - 2] + ys[index - 1] + ys[index]
}

day_last = function(step_index) {
    len = length(step_index)
    if(len %% 3) {
        stop('Incorrect step_index length:', len)
    }
    index = (1:(len/3)) * 3
    step_index[index]
}

day_average_all = function(ys) {
    len = length(ys)
    if(len %% 3) {
        stop('Incorrect ys length:', len)
    }
    index = (1:(len/3)) * 3
    (ys[index - 2] + ys[index - 1] + ys[index]) / 3
}

production_add_two = function(ys) {
    len = length(ys)
    if(len %% 2) {
        stop('Incorrect ys length:', len)
    }
    index = (1:(len/2)) * 2
    ys[index - 1] + ys[index]
}

production_average_two = function(ys) {
    len = length(ys)
    if(len %% 2) {
        stop('Incorrect ys length:', len)
    }
    index = (1:(len/2)) * 2
    (ys[index - 1] + ys[index]) / 2
}


#The following several functions may be combined at some point in the future.
oneplot = function(
                   filename,
                   outcome_fn,
                   primary_summary_fn,
                   ylim,
                   ylab,
                   main_title = NULL,
                   mask_fn = function(d) NA,
                   step_combiner = function(x) x,
                   ys_combiner = function(x) x
                   ) {
    png(paste(subdirectory, unique_id, '_', filename, '_', VERSION, '.png',
              sep = ''),
        height = 1000, width = 1000)

    #bit of a kludge, but should ensure sane limits
    ys = list()
    step_indices = list()
    for (i in 1:length(full_output_filenames)) {
        full_output = readRDS(full_output_filenames[i])
        fragments = unlist(strsplit(full_output_filenames[i], '/'))
        start_days = readRDS(paste0(fragments[1], '/start_days--', fragments[2]))
        mask = mask_fn(start_days)

        ys[[i]] = combine(full_output, outcome_fn, primary_summary_fn, mask)
        step_indices[[i]] = step_index
        len = length(step_indices[[i]])
        if(!identical(mask, NA)) {
            include = sapply(1:dim(mask)[1], function(i) sum(mask[i,]) != 0)
            ys[[i]] = ys[[i]][include]
            step_indices[[i]] = step_indices[[i]][include]
        }
        step_indices[[i]] = step_combiner(step_indices[[i]])
        ys[[i]] = ys_combiner(ys[[i]])
    }
    for(i in 1:length(full_output_filenames)) {
        if(i == 1) {
            par(mar = c(5,5,4,2))
            plot(step_indices[[i]], ys[[i]], type = 'l', col = colors[i],
                 ylim = c(min(ylim[1], min(sapply(ys, function(x) min(x)))),
                          max(ylim[2], max(sapply(ys, function(x) max(x))))
                        ),
                 #xlim = c(0, days),
                 lwd = 4,
                 xlab = "Day", ylab = ylab, cex.axis = 1.5, cex.lab = 1.5,
                 lty = ltys[i])
            title(main=main_title, cex.main = 3)
            existing_ticks = axTicks(1)

            #If the maximum x-value being plotted is far enough from the maximum
            #automatically-generated tick label, add a tick and corresponding
            #label at the maximum x-value. The value of .036 was found
            #experimentally, and may need to be adjusted for a different
            #resolution, different cex.axis, different taste, etc. Or just drop
            #this part if you'd rather just allow the data to extend beyond the
            #final tick (as it does on the y-axis).
            if((days - existing_ticks[length(existing_ticks)]) / days > .036) { 
                axis(1, at=days, cex.axis = 1.5)
            }
        } else {
            points(step_indices[[i]], ys[[i]], col = colors[i], lwd = 4, type = 'l', lty = ltys[i])
        }
    }
    legend("topright",inset = .06, row.names, lwd = 4,
           col = colors, lty = ltys, y.intersp = 1, cex = 1.5)
    dev.off()
    return(max(sapply(ys,max))) #to get maxes for the plots with forced same
                                #axes (in special purpose internal versions;
                                #kept here for consistency)
}

shiftwise_unavailable = function(data) {
    data[,'qn_absent',] 
}

shiftwise_scheduled = function(data) {
    data[,'qn_scheduled',]
}

shiftwise_unavailable_fraction = function(data) {
    data[,'qn_absent',] / data[,'qn_scheduled',]
}

shiftwise_short = function(data) {
    shiftwise_unavailable_fraction(data) > .15
}

shiftwise_production_loss = function(data) {
    fraction_available = 1 - shiftwise_unavailable_fraction(data)
    adjusted_fraction_available = pmin(fraction_available / 0.85, 1)
    fractional_production = adjusted_fraction_available^0.437
    fractional_loss = 1 - fractional_production

    fractional_loss * output_per_shift
}

new_infections = function(data) {
    data[,'new_infections',]
}

new_symptomatic_infections = function(data) {
    data[,'new_symptomatic_infections',]
}

new_unavailables = function(data) {
    data[,'new_unavailables',]
}

temperature_screening_cost = function(data) {
    thermometer_cost_each <- 20 # $20 per thermometer 
    KN95_cost <- 1 # $1 per mask per day 
    face_shield_cost <- 3 # $3 per face shield. Changing every 30 days. ($0.1/day) 
    ts_time <- 3 # 3 seconds for each screening
    ts_limit <- 5 #screening should be completed under 5 minute

    scheduled = shiftwise_scheduled(data)
    available = scheduled - shiftwise_unavailable(data)
    screeners = ceiling(scheduled) / (ts_limit * 60 / ts_time)
    ts_time <- available * ts_time / screeners / 3600   # Actual daily screening time in hours
    compensation <- ts_time * screeners * hourly_wage * 2 # have to pay the screeners, and the people being screened
    screener_training_cost = ceiling(N/100) * hourly_wage #max(screeners) * hourly_wage # 1hour training cost for screeners
    thermometer_cost <- max(screeners) * thermometer_cost_each

    initial_cost = screener_training_cost + thermometer_cost
    ongoing_cost = compensation + (KN95_cost + face_shield_cost/30) * screeners

    ongoing_cost[1] = ongoing_cost[1] + initial_cost

    ifelse(is.na(ongoing_cost), 0, ongoing_cost) #needs modification if we ever end up plotting over time
}

virus_testing_cost = function(data) {
    vt_kit <- 10 # $10 per test
    vt_time <- 1/4 # 15 minutes waiting assumed

    # Average wage compensation + kit cost over simulation
    vt <- data[,'tests',] # number of tests
    vt_cost <- vt * (vt_time * hourly_wage + vt_kit) # total cost

    vt_cost
}

vaccination_cost = function(data) {
    ############## Vaccination ############
    # 0.75 hour paid sick leave per vaccination  
    # no production loss
    data[,'doses',] * hourly_wage * 0.75
}

R0_reduction_cost = function(data, kludge_index) {
    face_shield <- 3 # $3 per face shield. Changing every month (30 days)
    KN95 <- 1 # $1 per N95 per shift
    air_cleaner <- 1000 # 1 air cleaner per 1000 sqft
    life <- 3 * 365 # 3year life of air_cleaner
    bi_available <- shiftwise_scheduled(data) - shiftwise_unavailable(data)
    bi_avilable = ifelse(is.na(bi_available), 0, bi_available)
    if(kludge_index == 8) {
        bi_cost = KN95 * bi_available
    } else if(kludge_index == 9 || (kludge_index == 10 && farm_or_facility == 'farm')) {
        bi_cost = ((KN95 + face_shield/30) * bi_available)
    } else if(kludge_index == 10) {
        bi_cost = ((KN95 + face_shield/30) * bi_available) + size/1000 * air_cleaner / life 
    } else {
        stop(kludge_index)
    }

    bi_cost
}

generate_intervention_expenses_function = function() {
    i = 0

    function(data) {
        i <<- i + 1
        i_ = ceiling(i / double_wrap_num_sims)
        if(i_ == 1) {
            array(0, c(dim(data)[1], dim(data)[3]))
        } else if(i_ == 2) {
            temperature_screening_cost(data)
        } else if(i_ %in% 3:5) {
            virus_testing_cost(data)
        } else if(i_ %in% c(6:7, 11:13)) {
            vaccination_cost(data)
        } else {
            R0_reduction_cost(data, i_)
        }
    }
}

source('vioplot-multi-h.R') #a little kludgey, but better balances the needs of the plots for different interventions
end_boxplot = function(
                       filename,
                       outcome_fn,
                       xlab,
                       average = FALSE,
                       xlim = NULL,
                       percent = FALSE,
                       main_title = NULL,
                       mask_fn = NULL,
                       function_ = boxplot,
                       ys_combiner = function(x) x,
                       pairwise_differences = FALSE, #implicitly vs. i = 1 for now; can elaborate later
                       run_mask = TRUE,
                       percent_differences = FALSE,
                       areaEqual = FALSE,
                       h = 12, #NULL
                       outlier_kludge = FALSE,
                       plus.minus.100 = FALSE
                       ) {
    png(paste(subdirectory, unique_id, '_', filename, '_', VERSION, '.png', sep = ''), height = 1000, width = 1000)
    if(sum(run_mask) == 0) { #can't do anything meaningful with no runs!
                             #but we still want to make the nature of the problem clear
        plot.new()
        if(is.null(main_title)) {
            title('NO SUCH RUNS')
        } else {
            title(paste0(main_title, ': NO SUCH RUNS'))
        }
        dev.off()
        return()
    }
    print(filename)
    means = NULL
    all_outcomes = NULL
    min_h = Inf
    max_h = 0
    hs = NULL
    for (i in 1:length(full_output_filenames)) {
        intervention_start = Sys.time()
        full_output = readRDS(full_output_filenames[i])
        full_output = full_output[,,run_mask, drop = FALSE] #added drop = FALSE in case only _1_ run is selected
                                                            #if full_output[,,run_mask] selects multiple runs _or_ zero runs, then the dimension is preserved
                                                            #but if full_output[,,run_mask] selects exactly 1 run, then without drop = FALSE, its dimension will be reduced
        
        fragments = unlist(strsplit(full_output_filenames[i], '/'))
        start_days = readRDS(paste0(fragments[1], '/start_days--', fragments[2]))
        start_days = start_days[run_mask]
        if((i == 1) && is.null(mask_fn)) {
            trivial_mask = matrix(TRUE, nrow = dim(full_output)[1], ncol = dim(full_output)[3])
            mask_fn = function(x) trivial_mask
        }
        mask = mask_fn(start_days)
        dimnames(full_output) = list(rep(NA, dim(full_output)[1]), colnames(full_output), rep(NA, dim(full_output)[3])) #kludge
        obtain_value = function(j) { #run index; note difference from approach in combine()
                                     #TBD: create a better function name (and
                                     #better intermediate vector names
            tryCatch({v = full_output[mask[,j], , j, drop = FALSE]}, # prevent reduction in dimensions, so the same outcome_fn can be used,
                error = function(e) {
                    print(e)
                    browser()
                }
            )
            vv = outcome_fn(v)
            vvv = ys_combiner(vv)
            len <<- length(vvv)
            sum(vvv)
        }
        
        final = sapply(1:(dim(full_output)[3]), obtain_value)
        
        if(average) {
            final = final / len
        }
        
        if(pairwise_differences) {
            if(i == 1) {
                this_h = NULL
                final_1 = final
            } else {
                if(percent_differences) {
                    final_this = (final - final_1) / final_1
                } else {
                    final_this = final - final_1
                }
                all_outcomes = rbind(all_outcomes, data.frame(intervention = row.names[i], outcome = final_this))
                if(all(final_this == final_this[1])) { #i.e., all equal (or none present)
                    this_h = 0
                } else {
                    tryCatch(
                        {
                            if(outlier_kludge) {
                                min_final_this_index = which(final_this == min(final_this))[1]
                                max_final_this_index = which(final_this == max(final_this))[1]
                                trimmed_final_this = final_this[-c(min_final_this_index, max_final_this_index)]
                                this_h = sm.density(trimmed_final_this, display='none')$h
                            } else {
                                this_h = sm.density(final_this, display='none')$h
                            }
                        },
                        error = function(e) {
                            print(e)
                            browser()
                        }
                    )
                }
                hs = c(hs, this_h)
                means[i - 1] = mean(final_this) 
            }
        } else {
            all_outcomes = rbind(all_outcomes, data.frame(intervention = row.names[i], outcome = final))
            if(all(final == 0)) {
                this_h = NULL
                hs = c(hs, 0)
            } else {
                this_h = sm.density(final, display='none')$h
                hs = c(hs, this_h)
            }
            means[i] = mean(final) 
        }
        if(!is.null(this_h)) {
            cat('\t', this_h, '\t', row.names[i], '\n')
            max_h = max(max_h, this_h)
            min_h = min(min_h, this_h)
        }
    }
    cat('\t\t', min_h, 'MIN', '\n\t\t', max_h, 'MAX\n')

    all_outcomes$intervention = factor(all_outcomes$intervention, levels = unique(all_outcomes$intervention), ordered = TRUE)

    if(pairwise_differences) {
        col = colors[-1]
        full_output_filenames = full_output_filenames[-1]
    } else {
        if(identical(function_, ecdfs)) {
            col = colors
        } else {
            col = c('white', colors[-1])
        }
    }

    par(mar = c(5,23,4,2))
    if(identical(function_, vioplot)) {
        par(cex.lab = 1.5)
    }
    if(percent) {
        if(identical(function_, vioplot)) {
            function_(outcome ~ intervention, data = all_outcomes, horizontal = TRUE, las = 1, xlab = xlab, ylim = xlim, col = col, cex.axis = 1.5, cex.names=1.5, cex.lab=1.5, ylab = '', na.action = na.pass, yaxt='n', areaEqual = areaEqual, h = hs / 5)
        } else {
            function_(outcome ~ intervention, data = all_outcomes, horizontal = TRUE, las = 1, xlab = xlab, ylim = xlim, col = col, cex.axis = 1.5, cex.names=1.5, cex.lab=1.5, ylab = '', na.action = na.pass, xaxt='n')
        }
        if(plus.minus.100) {
            axis(1, at = c(-1, -0.5, 0, 0.5, 1), paste0(lab=c('-100%', '-50%', '0%', '+50%', '+100%')), las=TRUE, cex.axis = 1.5, cex.lab=1.5)
        } else {
            axis(1, at=pretty(c(all_outcomes$outcome,xlim)), paste0(lab=pretty(c(all_outcomes$outcome,xlim)) * 100, ' %'), las=TRUE, cex.axis = 1.5, cex.lab=1.5)
        }
    } else {
        function_(outcome ~ intervention, data = all_outcomes, horizontal = TRUE, las = 1, xlab = xlab, ylim = xlim, col = col, cex.axis = 1.5, cex.names=1.5, cex.lab=1.5, ylab = '', na.action = na.pass, areaEqual = areaEqual, h = hs / 5)
    }
    title(main=main_title, cex.main = 3, outer = TRUE, line = -3)
    points(means, 1:length(full_output_filenames), cex =2, pch = 8)
    abline(v = 0)
    dev.off()
}

#parameters to match pattern in end_boxplot
ecdfs = function(f, data, horizontal, las, xlab, ylim, col, cex.axis, cex.names, cex.lab, ylab, na.action) {
    all_outcomes = data
    interventions = all_outcomes$intervention
    levels_ = levels(interventions)
    outcomes = all_outcomes$outcome
    n = length(levels_)
    xlim = c(min(outcomes), max(outcomes))
    for(i in 1:n) {
        this_intervention = levels_[i]
        these_outcomes = outcomes[interventions == this_intervention]
        plot(ecdf(these_outcomes), xlab = xlab, xlim = xlim, col = col[i], cex.axis = 1.5, cex.names=1.5, cex.lab=1.5, ylab = '', na.action = na.pass, main = '')
        if(i != n) {
            par(new = TRUE)
        }
    }
}

scatter_plot = function(filename,
                       outcome_fn_x,
                       xlab,
                       mask_x = NA,
                       outcome_fn_y,
                       ylab,
                       mask_y = NA,
                       main_title = NULL
                       ) {
    png(paste(subdirectory, unique_id, '_', filename, '_', VERSION, '.png', sep = ''), height = 1000, width = 1300)

    step_index_x = step_index_y = step_index

    if(!is.na(mask_x)[1]) {
        step_index_x = step_index[mask_x]
    }
    if(!is.na(mask_y)[1]) {
        step_index_y = step_index[mask_y]
    }

    ###
        
    means_x = numeric(length(full_output_filenames))        
    means_y = numeric(length(full_output_filenames))
    for (i in 1:length(full_output_filenames)) {
        full_output_x = full_output_y = readRDS(full_output_filenames[i])
        if(!is.na(mask_x)[1]) {
            full_output_x = full_output_x[mask_x,,]
        }
        if(!is.na(mask_y)[1]) {
            full_output_y = full_output_y[mask_y,,]
        }

        dimnames(full_output_x) = list(rep(NA, dim(full_output_x)[1]), colnames(full_output_x), rep(NA, dim(full_output_x)[3])) #kludge
        dimnames(full_output_y) = list(rep(NA, dim(full_output_y)[1]), colnames(full_output_y), rep(NA, dim(full_output_y)[3])) #kludge
        outcomes_x = outcome_fn_x(full_output_x)
        outcomes_x = apply(outcomes_x, 2, cumsum)
        outcomes_y = outcome_fn_y(full_output_y)
        outcomes_y = apply(outcomes_y, 2, cumsum)

        final_x = as.vector(outcomes_x[dim(full_output_x)[1],])
        final_y = as.vector(outcomes_y[dim(full_output_y)[1],])
        
        means_x[i] = mean(final_x, na.rm = TRUE)
        means_y[i] = mean(final_y, na.rm = TRUE)
    }

    par(mar = c(5,5,4,32), xpd=TRUE)
    plot(means_x, means_y, xlab = xlab, ylab = ylab, col = colors, cex.axis = 2, #cex.names=1.5,
         cex.lab=2, pch = ltys, lwd = 12)
    title(main=main_title, cex.main = 3)
    legend("topright", row.names, lty = 0, lwd = 12,
           col = colors, pch = ltys, y.intersp = 1, cex = 2, inset = c(-0.585,0))
    dev.off()

}



first_x_boxplot = function(
                           filename,
                           outcome_fn,
                           xlab,
                           xlim = NULL,
                           mask = NA,
                           function_ = boxplot
                           ) {
    if(!is.null(filename)) {
        png(paste(subdirectory, unique_id, '_', filename, '_', VERSION, '.png', sep = ''), height = 1000, width = 1000)
    }

    if(!is.na(mask)[1]) {
        step_index = step_index[mask]
    }

    means = numeric(length(full_output_filenames))
    for (i in 1:length(full_output_filenames)) {
        full_output = readRDS(full_output_filenames[i])

        if(!is.na(mask)[1]) {
            full_output = full_output[mask,,]
        }

        dimnames(full_output) = list(rep(NA, dim(full_output)[1]), colnames(full_output), rep(NA, dim(full_output)[3])) #kludge
        outcomes = outcome_fn(full_output)
        
        first = apply(outcomes, 2, function(v) ifelse(length(which(v)) > 0, step_index[which(v)[1]], NA))
        
        means[i] = mean(first, na.rm = TRUE)

        if(i == 1) {
            all_outcomes = data.frame(intervention = row.names[i], outcome = first)
        } else {
            all_outcomes = rbind(all_outcomes, data.frame(intervention = row.names[i], outcome = first))
        }
    }

    all_outcomes$intervention = factor(all_outcomes$intervention, levels = unique(all_outcomes$intervention), ordered = TRUE)

    par(mar = c(5,23,4,2))
    function_(outcome ~ intervention, data = all_outcomes, horizontal = TRUE, las = 1, xlab = xlab, ylim = xlim, col = c('white', colors[-1]), cex.axis = 1.5, cex.names=1.5, cex.lab=1.5, ylab = '', na.action = na.pass)
    debug_all_outcomes <<- all_outcomes
    debug_means <<- means
    points(means, 1:length(full_output_filenames), cex =2, pch = 8)
    if(!is.null(filename)) {
        dev.off()
    }
}

end_barplot = function(
                       filename, 
                       outcome_fn, 
                       xlab, 
                       summary_fn,
                       xlim = NULL, 
                       percent = FALSE,
                       main_title = NULL,
                       mask_fn = NULL,
                       ys_combiner = max
                       ) {
    if(!is.null(filename)) {
        png(paste(subdirectory, unique_id, '_', filename, '_', VERSION, '.png', sep = ''), height = 1000, width = 1000)
    }
    

    all_outcomes = numeric(length(full_output_filenames))
    names(all_outcomes) = row.names
    for (i in 1:length(full_output_filenames)) {
        full_output = readRDS(full_output_filenames[i])
        fragments = unlist(strsplit(full_output_filenames[i], '/'))
        start_days = readRDS(paste0(fragments[1], '/start_days--', fragments[2]))

        if((i == 1) && is.null(mask_fn)) {
            trivial_mask = matrix(TRUE, nrow = dim(full_output)[1], ncol = dim(full_output)[3])
            mask_fn = function(x) trivial_mask
        }
        mask = mask_fn(start_days)
        dimnames(full_output) = list(rep(NA, dim(full_output)[1]), colnames(full_output), rep(NA, dim(full_output)[3])) #kludge
        obtain_value = function(j) { #run index; note difference from approach in combine()
                                     #TBD: create a better function name (and
                                     #better intermediate vector names
            v = full_output[mask[,j], , j, drop = FALSE] # prevent reduction in dimensions, so the same outcome_fn can be used
            vv = outcome_fn(v)
            vvv = ys_combiner(vv)
            len <<- length(vvv)
            sum(vvv)
        }
        outcomes = sapply(1:(dim(full_output)[3]), obtain_value)

        fraction = summary_fn(outcomes)
        
        all_outcomes[i] = fraction
    }
    par(mar = c(5,23,4,2))
    if(percent) {
        par(xaxt="n")
    }
    barplot(all_outcomes, horiz = TRUE, las = 1, xlab = xlab, xlim = xlim, col = colors, cex.axis = 1.5, cex.names=1.5, cex.lab=1.5)#, main_title = main_title)
    title(main=main_title, cex.main = 3, outer = TRUE, line = -3)
    if(percent) {
        par(xaxt='s')
        axis(1, at=pretty(c(all_outcomes,xlim)), paste0(lab=pretty(c(all_outcomes,xlim)) * 100, ' %'), las=TRUE, cex.axis = 1.5, cex.lab=1.5)
    }
    if(!is.null(filename)) {
        dev.off()
    }
}


production_shifts_mask_fn = function(start_days) {
    sapply(start_days, production_shifts)
}

cleaning_shifts_mask_fn = function(start_days) {
    sapply(start_days, cleaning_shifts)
}
work_shifts_mask_fn = function(start_days) {
    sapply(start_days, work_shifts)
}

end_barplot(filename = 'Symptomatic-Fraction-Non-Zero', outcome_fn = symptomatic, xlab = 'Fraction of runs where symptomatic infections > 0', summary_fn = mean, xlim = c(0, 1), percent = TRUE, main_title = '(C) Fraction of Runs > 0', mask_fn = NULL, ys_combiner = function(x) sum(x) > 0)
end_barplot(filename = 'Unavailable-production-Fraction-Non-Zero', outcome_fn = shiftwise_unavailable, xlab = 'Fraction of runs where worker-shifts missed > 0', summary_fn = mean, xlim = c(0, 1), percent = TRUE, main_title = '(C) Fraction of Runs > 0', mask_fn = production_shifts_mask_fn, ys_combiner = function(x) sum(x) > 0)

oneplot('Symptomatic', symptomatic, mean, c(0,0), paste('People symptomatically infected', sep = ''), step_combiner = day_average_all, ys_combiner = day_average_all, main_title = '(B) Mean Prevalence at Each Time Point')

oneplot('Symptomatic-incidence', new_symptomatic_infections, mean, c(0,0), paste('Incidence of symptomatic infection', sep = ''), step_combiner = day_average_all, ys_combiner = day_average_all, main_title = '(A) Mean Incidence at Each Time Point')


if(farm_or_facility == 'facility') {
    if(supervisors > 1) {
        production_step_combiner = production_average_two
        production_ys_combiner = production_add_two
        all_step_combiner = day_average_all
        all_ys_combiner = day_add_all
        n_shifts = 2
    } else {
        production_step_combiner = function(x) x
        production_ys_combiner = function(x) x
        all_step_combiner = production_average_two
        all_ys_combiner = production_add_two
        n_shifts = 1
    }
} else {
    production_step_combiner = function(x) x
    production_ys_combiner = function(x) x
    n_shifts = 1
}
df1 = readRDS(full_output_filenames[1])
zero_symptomatic_mask = apply(df1[,'new_symptomatic_infections',], 2, sum) == 0
baseline_fragments = unlist(strsplit(full_output_filenames[1], '/'))
sd1 = readRDS(paste0(baseline_fragments[1], '/start_days--', baseline_fragments[2]))
production_shifts_mask = production_shifts_mask_fn(sd1)
obtain_baseline_production_shifts_absences = function(j) sum(df1[production_shifts_mask[,j],'qn_absent', j]) 
zero_unavailable_mask = sapply(1:double_wrap_num_sims, obtain_baseline_production_shifts_absences) == 0
zero_unavailable_mask


end_boxplot('non-zero-pairwise-percent-differences-Total-Symptomatic-Infections-violin--cut-and-trimmed', new_symptomatic_infections, xlab = paste('Total symptomatic infections: pairwise fractional difference'), average = FALSE, main_title = '(F) P. F. Change, Non-Zero Baseline Runs', function_ = vioplot, pairwise_differences = TRUE, run_mask = !zero_symptomatic_mask, percent_differences = TRUE, percent = TRUE, xlim = c(-1, 1), outlier_kludge = TRUE, , plus.minus.100 = TRUE)

print('ping')
end_boxplot('non-zero-pairwise-percent-differences-Total-Unavailable-production-violin', shiftwise_unavailable, xlab = paste('Total Worker-shifts Missed: pairwise fractional difference'), average = FALSE, main_title = '(F) P. F. Change, Non-Zero Baseline Runs', mask_fn = production_shifts_mask_fn, percent = TRUE, function_ = vioplot, pairwise_differences = TRUE, run_mask = !zero_unavailable_mask, percent_differences = TRUE)
print('pong')

end_boxplot('non-zero-pairwise-differences-Total-Symptomatic-Infections-violin', new_symptomatic_infections, xlab = paste('Total symptomatic infections: pairwise difference'), average = FALSE, main_title = '(E) P. Differences, Non-Zero Baseline Runs', function_ = vioplot, pairwise_differences = TRUE, run_mask = !zero_symptomatic_mask)
print('pung')


end_boxplot('zero-pairwise-differences-Total-Unavailable-production-violin', shiftwise_unavailable, xlab = paste('Total worker-shifts missed: pairwise difference'), average = FALSE, main_title = '(D) P. Differences, Zero Baseline Runs', mask_fn = production_shifts_mask_fn, function_ = vioplot, pairwise_differences = TRUE, run_mask = zero_unavailable_mask)
end_boxplot('non-zero-pairwise-differences-Total-Unavailable-production-violin', shiftwise_unavailable, xlab = paste('Total Worker-shifts missed: pairwise difference'), average = FALSE, main_title = '(E) P. Differences, Non-Zero Baseline Runs', mask_fn = production_shifts_mask_fn, function_ = vioplot, pairwise_differences = TRUE, run_mask = !zero_unavailable_mask)
print('pang')

oneplot('Unavailable-production', shiftwise_unavailable, mean, c(0,0), paste('Mean number at each time point', sep = ''), mask_fn = production_shifts_mask_fn, step_combiner = production_step_combiner, ys_combiner = production_ys_combiner, main_title = '(A) People Unavailable to Work Their Scheduled Shift')
print('pong')


main_title = ''

sys_time_start = Sys.time()
end_boxplot('Total-Unavailable-production-violin', shiftwise_unavailable, xlab = paste('Total worker-shifts missed'), average = FALSE, main_title = '(B) Cumulative Worker-Shifts Missed', mask_fn = production_shifts_mask_fn, function_ = vioplot)

end_boxplot('Total-Symptomatic-Infections', new_symptomatic_infections, xlab = paste('Total symptomatic infections'), average = FALSE, main_title = '')
end_boxplot('Total-Symptomatic-Infections-violin', new_symptomatic_infections, xlab = paste('Total symptomatic infections'), average = FALSE, main_title = '(D) Cumulative Incidence, Distribution', function_ = vioplot)
#end_boxplot('Total-Symptomatic-Infections-prevalence-violin', symptomatic, xlab = paste('Total Worker-Days Symptomatically Infected'), average = FALSE, ys_combiner = day_average_all, main_title = '(F) Cumulative Prevalence, distribution', function_ = vioplot)
#end_boxplot('pairwise-differences-Total-Symptomatic-Infections-violin', new_symptomatic_infections, xlab = paste('Total Symptomatic Infections (among', N, 'total workers)'), average = FALSE, main_title = '(C) P. Differences, all runs', function_ = vioplot, pairwise_differences = TRUE)
#browser()
#end_boxplot('Fraction-Short-production', shiftwise_short, xlab = 'Percentage of Production Shifts Short (> 15% of workers absent)', average = TRUE, xlim = c(0,1), percent = TRUE, main_title = main_title, mask = production_shifts)
end_boxplot('v4b-Fraction-Short-production-violin', shiftwise_short, xlab = 'Percentage of production shifts with a shortage (> 15% of workers absent)', average = TRUE, xlim = c(0,0.2), percent = TRUE, main_title = '(D) Fraction of Production Shifts Short', mask_fn = production_shifts_mask_fn, function_ = vioplot)

"oneplot('v4b-Production-Loss', shiftwise_production_loss, mean, c(0,0), 'Production Loss (Dollars ($) per production shift)', mask = production_shifts)
end_boxplot('v4b-Total-Production-Loss', shiftwise_production_loss, xlab = 'Total Production Loss in Dollars ($)', mask = production_shifts)"
end_boxplot('Total-Production-Loss-violin', shiftwise_production_loss, xlab = 'Total production loss in dollars ($)', mask_fn = production_shifts_mask_fn, ys_combiner = production_ys_combiner, function_ = vioplot, main_title = '(B) Total Production Loss')
end_boxplot('Total-Intervention-Expenses-violin', generate_intervention_expenses_function(), xlab = 'Total intervention expenses in dollars ($)', function_ = vioplot, main_title = '(A) Total Intervention Expenses')



intervention_expenses_function = generate_intervention_expenses_function()
#below is massively kludged, to deal with production loss fn not handling
#cleaning shifts, and cost needing to handle all shifts
#follow by the need to have two dimensions for handling in end_boxplot
g = function(data) {
    #TBD: Make this more elegant
    ad_hoc_production_mask = rep(c(TRUE, TRUE, FALSE), days)
    fd = shiftwise_production_loss(data[ad_hoc_production_mask,,, drop = FALSE])
    fd = ifelse(is.na(fd), 0, fd)
    r = intervention_expenses_function(data)
    r[ad_hoc_production_mask] = r[ad_hoc_production_mask] + fd
    r
}
end_boxplot('Total-Cost-violin', g, xlab = 'Total cost (intervention expenses + production losses) in dollars ($)', function_ = vioplot, main_title = '(C) Total Cost')


ANALYZE = FALSE

} #analyze_fn

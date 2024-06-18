# update-immunity-function.R is part of Food INdustry CoViD Control Tool
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


#going to start this out the comparatively easy way: assume all agents get the same event, and times are preset
update_immunity = function(all_agents, mask, event, times, recovered_complete_protection_time, net_protection, infection_protection) {#, debug_string = NULL) {
    if(length(event) > 1) {
        stop('Multiple events in a single update_immunity call: ', paste(event, collapse = ', '))
    }
    if(!any(mask)) { #there's nothing to do
        return(all_agents)
    }
    agents = all_agents[mask,] #we have to keep the full data frame in order to return it
    times = times[mask] #these, on the other hand, we don't have to worry about
    if(any(is.na(times) | times == Inf | times == -Inf)) {
        stop('Invalid event times: ', paste(times, collapse = ', '), ' for event: ', event)
    }

    atlie = agents$time_last_immunity_event
    ais = agents$immune_status
    pi_ = net_protection(agents, times)
    pii_ = infection_protection(agents, times)
    if(event %in% c('V1', 'V2', 'B')) {
        agents[, paste0('time_', event)] = times
        agents[, 'time_last_immunity_event'] = times
        agents[, 'previous_immunity'] = pi_
        agents[, 'previous_infection_immunity'] = pii_
        agents$immune_status = ifelse(ais %in% c('FS', 'V1', 'V2', 'B'),
            event,
            ifelse(ais %in% c('R', 'H_V1', 'H_V2', 'H_B'),
                paste0('H_', event),
                NA #flagging for error
            )
        )
        if(any(is.na(agents$immune_status))) {
            stop('Unknown immune status.')
        }
        agents$vax_status = event
    } else if(event == 'R') {
        agents$time_R = agents$time_last_immunity_event = times
        agents$previous_immunity = pi_
        agents$previous_infection_immunity = pii_
        agents$immune_status = ifelse(ais %in% c('FS', 'R'),
            'R',
            ifelse(ais %in% c('V1', 'H_V1'),
                'H_V1',
                ifelse(ais %in% c('V2', 'H_V2'),
                    'H_V2',
                    ifelse(ais %in% c('B', 'H_B'),
                        'H_B',
                        NA
                    )
                )
            )
        )
    } else {
        stop('Invalid event: ', event)
    }
    #if(!is.null(debug_string)) {
    #    agents$update_record = paste(agents$update_record, debug_string)
    #}
    all_agents[mask,] = agents
    all_agents
}

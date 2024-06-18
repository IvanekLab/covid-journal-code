# custom-contacts-gen-general.R is part of Food INdustry CoViD Control
# Tool (FInd CoV Control), version 3.0.
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


facility_contacts_gen = function(workers_per_line = 10,
                                 n_lines = 3,
                                 n_production_shifts = 2, # or 1
                                 n_shift_floaters = 10,
                                 n_cleaners = 10,
                                 n_all_floaters = 10) {

    tt = matrix(1, nrow = workers_per_line, ncol = workers_per_line)
    ttd = diag(1, nrow = workers_per_line)

    crew = tt - ttd
    diff_crew = tt * 0.1 # for different lines

    for(i in 1:n_lines) {
        for(j in 1:n_lines) {
            if(i == j) {
                this_submatrix = crew
            } else {
                this_submatrix = diff_crew
            }
            if(j == 1) {
                these_rows = this_submatrix
            } else {
                these_rows = cbind(these_rows, this_submatrix)
            }
        }
        if(i == 1) {
            team_minus_supervisor = these_rows
        } else {
            team_minus_supervisor = rbind(team_minus_supervisor, these_rows)
        }
    }
    
    team_minus_sup = team_minus_supervisor
    
    #possibly, this should just be a static number. but keeping with the current approach for now . . .
    suprow = matrix((rowSums(crew)[1] + (n_lines - 1) * rowSums(diff_crew)[1]) / (workers_per_line * n_lines - 1), nrow = 1, ncol = workers_per_line * n_lines)
    supcol = t(suprow)
    
    team = rbind(cbind(0, suprow),
                 cbind(supcol, team_minus_sup))
    
    shift_floaters_contact_rate = rowSums(team)[1] / (workers_per_line * n_lines)
    
    shift_floaters = (matrix(1, nrow = n_shift_floaters, ncol = n_shift_floaters) - diag(1, nrow = n_shift_floaters)) * shift_floaters_contact_rate
    
    shift_floater_tall = matrix(shift_floaters_contact_rate, nrow = (workers_per_line * n_lines + 1), ncol = n_shift_floaters)
    shift_floater_wide = t(shift_floater_tall)
    
    production_shift = rbind(cbind(team, shift_floater_tall),
                             cbind(shift_floater_wide, shift_floaters))
        
    #NB (important): This and some other lines can choke on low, but
    #potentially valid, numbers (e.g., here if n_cleaners is 1, NaNs end up
    #propagating (and if it's 0, things can get even weirder)).
    #This can be avoided by adhering to the lower limits of the possible ranges
    #listed in the manuscript.
    cleaning_shift_contact_rate = rowSums(production_shift)[1] / (n_cleaners - 1)
    
    cleaning_shift = matrix(cleaning_shift_contact_rate, nrow = n_cleaners, ncol = n_cleaners) - diag(cleaning_shift_contact_rate, nrow = n_cleaners)
    
    production_shift_size = 1 + workers_per_line * n_lines + n_shift_floaters
    
    all_floaters_contact_rate = rowSums(cleaning_shift)[1] / (n_production_shifts * production_shift_size + n_cleaners + (1 + n_production_shifts) * (n_all_floaters) - (n_all_floaters + 1))
    
    all_floaters = matrix(all_floaters_contact_rate, nrow = n_all_floaters, ncol = n_all_floaters) - diag(all_floaters_contact_rate, nrow = n_all_floaters)
    
    
    almost_all_size = production_shift_size * (n_production_shifts) + n_cleaners + n_all_floaters
    
    
    ps1_1x1 = production_shift
    ps1_1x2 = matrix(0, nrow = production_shift_size, ncol = production_shift_size * (n_production_shifts > 1))
    ps1_1xC = matrix(0, nrow = production_shift_size, ncol = n_cleaners)
    ps1_1xF = matrix(all_floaters_contact_rate, nrow = production_shift_size, ncol = n_all_floaters)
    
    ps1_1 = cbind(ps1_1x1, ps1_1x2, ps1_1xC, ps1_1xF)
    
    ps1_2 = matrix(0, nrow = production_shift_size * (n_production_shifts > 1), ncol = almost_all_size)
    ps1_C = matrix(0, nrow = n_cleaners, almost_all_size)
    
    ps1_Fx1 = matrix(all_floaters_contact_rate, nrow = n_all_floaters, ncol = production_shift_size)
    ps1_Fx2 = matrix(0, nrow = n_all_floaters, ncol = production_shift_size * (n_production_shifts > 1))
    ps1_FxC = matrix(0, nrow = n_all_floaters, ncol = n_cleaners)
    ps1_FxF = all_floaters
    
    ps1_F = cbind(ps1_Fx1, ps1_Fx2, ps1_FxC, ps1_FxF)
    
    production_shift_1_minus_manager = rbind(ps1_1,
                                             ps1_2,
                                             ps1_C,
                                             ps1_F)
    
    ps1_Mx1 = matrix(all_floaters_contact_rate, nrow = 1, ncol = production_shift_size)
    ps1_Mx2 = matrix(0, nrow = 1, ncol = production_shift_size * (n_production_shifts > 1))
    ps1_MxC = matrix(0, nrow = 1, ncol = n_cleaners)
    ps1_MxF = matrix(all_floaters_contact_rate, nrow = 1, ncol = n_all_floaters)
    
    ps1_M = cbind(ps1_Mx1, ps1_Mx2, ps1_MxC, ps1_MxF)
    
    production_shift_1 = rbind(cbind(0, ps1_M),                                 
                               cbind(t(ps1_M), production_shift_1_minus_manager))
    
    if(n_production_shifts == 1) {
        production_shift_2 = matrix(0, nrow = 1 + almost_all_size, ncol =  1 + almost_all_size)
    } else if(n_production_shifts == 2) {
        ps2_1 = ps1_2
        ps2_2 = cbind(ps1_1x2, ps1_1x1, ps1_1xC, ps1_1xF) #first two swapped from above
        ps2_C = ps1_C
    
        ps2_F = cbind(ps1_Fx2, ps1_Fx1, ps1_FxC, ps1_FxF) #first two swapped from above
    
        production_shift_2_minus_manager = rbind(ps2_1,
                                                 ps2_2,
                                                 ps2_C,
                                                 ps2_F)
    
        ps2_Mx1 = ps1_Mx2
        ps2_Mx2 = ps1_Mx1
        ps2_MxC = ps1_MxC
        ps2_MxF = ps1_MxF
    
        ps2_M = cbind(ps2_Mx1, ps2_Mx2, ps2_MxC, ps2_MxF)
    
        production_shift_2 = rbind(cbind(0, ps2_M),                                 
                                   cbind(t(ps2_M), production_shift_2_minus_manager))
    } else {
        stop('n_teams must be 1 or 2')
    }
    
    C_1 = matrix(0, nrow = production_shift_size, ncol = almost_all_size)
    C_2 = ps1_2
    
    C_Cx1 = matrix(0, nrow = n_cleaners, ncol = production_shift_size)
    C_Cx2 = matrix(0, nrow = n_cleaners, ncol = production_shift_size * (n_production_shifts > 1))
    C_CxC = cleaning_shift
    C_CxF = matrix(all_floaters_contact_rate, nrow = n_cleaners, ncol = n_all_floaters)
    
    C_C = cbind(C_Cx1, C_Cx2, C_CxC, C_CxF)
    
    C_Fx1 = matrix(0, nrow = n_all_floaters, ncol = production_shift_size)
    C_Fx2 = ps1_Fx2
    C_FxC = t(C_CxF)
    C_FxF = ps1_FxF
    
    C_F = cbind(C_Fx1, C_Fx2, C_FxC, C_FxF)
    
    cleaning_shift_without_manager = rbind(C_1,
                                           C_2,
                                           C_C,
                                           C_F)
    
    C_Mx1 = matrix(0,nrow = 1, ncol = production_shift_size)
    C_Mx2 = matrix(0,nrow = 1, ncol = production_shift_size * (n_production_shifts > 1))
    C_MxC = matrix(all_floaters_contact_rate, nrow = 1, ncol = n_cleaners)
    C_MxF = matrix(all_floaters_contact_rate, nrow = 1, ncol = n_all_floaters)
    
    C_M = cbind(C_Mx1, C_Mx2, C_MxC, C_MxF)
    
    cleaning_shift_full = rbind(cbind(0, C_M),
                                cbind(t(C_M), cleaning_shift_without_manager))
    
    shift_sum = production_shift_1 + production_shift_2 + cleaning_shift_full
    png('sum.png', height=1000, width=1000) #probably disable this later
    image(shift_sum, main = 'Sum across all shifts')
    dev.off()
    list(production_shift_1 = production_shift_1,
         production_shift_2 = production_shift_2,
         cleaning_shift_full = cleaning_shift_full,
         shift_sum = shift_sum)
}

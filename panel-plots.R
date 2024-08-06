# panel-plots.R is part of Food INdustry CoViD Control Tool
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

library(png)
library(ggplot2)
library(grid)
library(cowplot)

dir.create('figures-after-additional-fixes/')
filenames = c(
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Symptomatic-incidence_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Symptomatic_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Symptomatic-Fraction-Non-Zero_2.2.0.png',

    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Total-Symptomatic-Infections-violin_2.2.0.png',



    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_non-zero-pairwise-differences-Total-Symptomatic-Infections-violin_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_non-zero-pairwise-percent-differences-Total-Symptomatic-Infections-violin--cut-and-trimmed_2.2.0.png'
)

png('figures-2024-04-29/both-figure-2-rearranged-stealing-issue-resolved-ABM-1000x--cut-and-trimmed.png', width = 2000, height = 3000)
l = NULL
for(filename in filenames) {
    this_image = readPNG(filename)
    this_plot = ggplot() + annotation_custom(
        grid::rasterGrob(
            this_image,
            width = unit(1, 'npc'),
            height = unit(1, 'npc')
        )
    )
    l = c(l, list(this_plot))
}
print(plot_grid(plotlist = l, nrow = 3, byrow = FALSE))
dev.off()


#stop('Got the first plot!')

filenames = c(
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Unavailable-production_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Total-Unavailable-production-violin_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Unavailable-production-Fraction-Non-Zero_2.2.0.png',

    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_zero-pairwise-differences-Total-Unavailable-production-violin_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_non-zero-pairwise-differences-Total-Unavailable-production-violin_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_non-zero-pairwise-percent-differences-Total-Unavailable-production-violin_2.2.0.png'
)



png('figures-2024-04-29/both-figure-3-rearranged-stealing-issue-resolved-ABM-1000x.png', width = 2000, height = 3000)
l = NULL
for(filename in filenames) {
    if(is.na(filename)) {
        this_plot = ggplot()
    } else {
        this_image = readPNG(filename)
        this_plot = ggplot() + annotation_custom(
            grid::rasterGrob(
                this_image,
                width = unit(1, 'npc'),
                height = unit(1, 'npc')
            )
        )
        l = c(l, list(this_plot))
    }
}
print(plot_grid(plotlist = l, nrow = 3, byrow = FALSE))
dev.off()

filenames = c(
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Total-Intervention-Expenses-violin_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Total-Production-Loss-violin_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_Total-Cost-violin_2.2.0.png',
    '2024-flexible-housing/stealing-issue-resolved-ABM-1000xfacility-dormitory_TRUE-community_TRUE-vaccinated_TRUE-recovered_TRUE_v4b-Fraction-Short-production-violin_2.2.0.png'
)

png('figures-2024-04-29/both-figure-4-rearranged-stealing-issue-resolved-ABM-1000x.png', width = 2000, height = 2000)
l = NULL
for(filename in filenames) {
    this_image = readPNG(filename)
    this_plot = ggplot() + annotation_custom(
        grid::rasterGrob(
            this_image,
            width = unit(1, 'npc'),
            height = unit(1, 'npc')
        )
    )
    l = c(l, list(this_plot))
}
print(plot_grid(plotlist = l, nrow = 2, byrow = FALSE))#2))
dev.off()


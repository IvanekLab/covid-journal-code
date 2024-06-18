# heatmaps.R is part of Food INdustry CoViD Control Tool
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

library(viridis)
library(fields)
civ = cividis(1000)
civ[1] = '#000000'
source('ContactsGen.R')
farm = ContactsGen(rep(3,3), rep(10,9), example_rates, 1)
#so farm structure is 1,(2,(3,4:13),(14,15:24),(25))
#leaving farm structure alone for now
png('farm-heatmap.png', height = 1000, width = 1000)
#previously had parameter col = civ[1:513] (probably for scaling?) only for this one plot
image.plot(farm, col = civ, xaxt = 'n', yaxt = 'n', legend.lab = 'Relative Contact Rate', legend.cex = 4, legend.width = 3, axis.args = list(lwd.ticks = 0, at = c(0, max(farm)), labels = rep('',2)), zlim = c(0,max(farm)))
dev.off()
source('custom-contacts-gen-general.R')
shift_sum = facility_contacts_gen()$shift_sum
v = c(1:22, 33:42, 23:32,94:103,43:63,74:83,64:73,84:93)
ssra = shift_sum[v,v]
png('facility-heatmap.png', height = 1000, width = 1000)
image.plot(ssra, col = civ, xaxt = 'n', yaxt = 'n', legend.lab = 'Relative Contact Rate', legend.cex = 4, legend.width = 3, axis.args = list(lwd.ticks = 0, at = c(0, max(shift_sum)), labels = rep('',2)), zlim = c(0,max(shift_sum)))
dev.off()

#> max(farm/sum(farm))/max(shift_sum/sum(shift_sum))
#[1] 0.5131786

png('facility-heatmap-original-order.png', height = 1000, width = 1000)
image.plot(shift_sum, col = civ, xaxt = 'n', yaxt = 'n', legend.lab = 'Relative Contact Rate', legend.cex = 4, legend.width = 3, axis.args = list(lwd.ticks = 0, at = c(0, max(shift_sum)), labels = rep('',2)), zlim = c(0,max(shift_sum)))
dev.off()

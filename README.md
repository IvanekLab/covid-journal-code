Food INdustry CoViD Control Tool (FInd CoV Control) version 3.0, with the
exception of the file vioplot-multi-h.R (see below), is Copyright (C) 2020-2024
Cornell University. This program is free software. You can redistribute it
and/or modify it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

vioplot-multi-h.R is adapted from vioplot.R, which is Copyright (C) 2004
Daniel Adler. You can redistribute it and/or modify it under the terms of
3-Clause BSD License (included both at the top of that file and as the file
LICENSE-BSD), with the following parameters:
YEAR: 2004
COPYRIGHT HOLDER: Daniel Adler
ORGANIZATION: University of Goettingen
 
This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.

If this collection of files is unmodified, this license will be in the file
LICENSE-GPL.

This code, apart from vioplot-multi-h.R, was written by Christopher Henry,
Renata Ivanek, Michelle Wemette, and Ece Bulut. vioplot-multi-h.R is based on
vioplot.R, written by Daniel Adler; subsequent modifications were written by
Christopher Henry. Additional work on writing and reviewing documentation,
catching erroneous output, making suggestions, etc. was done by Sarah I. Murphy,
Sebastian Llanos-Soto, Claire Zoellner, and Diane Wetherington. Additional
thanks is due to iFoodDecisionSciences (iFoodDS), who developed a web interface
for this model, and who have helped with background information, debugging,
suggestions and requests for output features, etc.

This model was developed for:
R version 4.3.2 (2023-10-31) -- "Eye Holes"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

The following package versions were used in development and testing:
doParallel 1.0.17
foreach 1.5.2
Matrix 1.6.1.1
vioplot 0.4.0
abind 1.4.5
viridis 0.6.5
fields 15.2
png 0.1.8
ggplot2 3.4.4
grid 4.3.2
cowplot 1.1.1

Some of these may be safe to omit, with light code editing.

Efforts have been made to keep comments and other documentation up to date,
but it is possible that some comments may be out of date, in accordance with
the principle that documentation is *always* out of date.

Suggested citation: “C. Henry, D. Tang, E. Bulut, S.I. Murphy, C. Zoellner, A.
Adalja, D. Wetherington, M. Wiedmann, S. Alcaine, R. Ivanek. Infection control
strategies in essential industries: using COVID-19 in the food industry to model
economic and public health trade-offs. (Under review)”

Usage:
To generate Figures 1.B.iii, 2, 3, 4, 5, and S1, from branch 'main':
source('iFoodDS-wrapper.R')
source('panel-plots.R')
source('rpart-scratch-work.R')
source('heatmaps.R')

Additional specific usage details to be added.

For more general use, from branch 'main', modify iFoodDS-wrapper.R to remove
calls to the full_run(...), and to add calls of your own.

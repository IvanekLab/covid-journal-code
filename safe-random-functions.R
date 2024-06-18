# safe-random-functions.R is part of Food INdustry CoViD Control Tool
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

#safe versions of runif and rbinom (restricted to a Bernoulli distribution
#since that's all I actually need) to allow a guarantee that the same number of
#calls to the underlying PRNG (hence rand() calls) will be made for a given N,
#regardless of the other parameters used. This is necessary because otherwise, a
#call to runif(N, 0, 0), for example, will not actually make any calls to the
#underlying PRNG, while a call to runif(N, 0, 1) will. This also allows for the
#possibility of eventually adding some manner of latin hypercube functionality.


random_functions_and_counter = function() {
    set_seed = NA
    runif_0_1_calls_counter = NA

    safe_set_seed = function(seed) {
        set.seed(seed)
        set_seed <<- seed
        runif_0_1_calls_counter <<- 0
    }

    sunif = function(N, min, max) {
        X = runif(N, 0, 1)
        runif_0_1_calls_counter <<- runif_0_1_calls_counter + N
        X * (max - min) + min
    }

    sbern = function(N, prob) {
        X = runif(N, 0, 1)
        runif_0_1_calls_counter <<- runif_0_1_calls_counter + N
        1 * (X < prob) #strict < ensures that prob = 0 always gives 0,
                       #prob = 1 always gives 1
    }

    sgamma = function(N, shape, scale) {
        X = runif(N, 0, 1)
        runif_0_1_calls_counter <<- runif_0_1_calls_counter + N
        qgamma(X, shape=shape, scale = scale)
    }

    slnorm = function(N, meanlog, sdlog) {
        X = runif(N, 0, 1)
        runif_0_1_calls_counter <<- runif_0_1_calls_counter + N
        qlnorm(X, meanlog, sdlog)
    }

    print_rand_state = function(s) {
        cat('\n\n', s, '\nseed:', set_seed, '\ncalls:', runif_0_1_calls_counter,
            '\ncurrent state, abbreviated:', .Random.seed[1:3], '\n\n')
    }

    list(
         safe_set_seed = safe_set_seed,
         sunif = sunif,
         sbern = sbern,
         sgamma = sgamma,
         slnorm = slnorm,
         print_rand_state = print_rand_state
    )
}

if(!exists('SAFE_RANDOM_FUNCTIONS_INITIALIZED')) {
    rfac_l = random_functions_and_counter()
    safe_set_seed = rfac_l[['safe_set_seed']]
    sunif = rfac_l[['sunif']]
    sbern = rfac_l[['sbern']]
    sgamma = rfac_l[['sgamma']]
    slnorm = rfac_l[['slnorm']]
    sbinom = rfac_l[['sbinom']]
    print_rand_state = rfac_l[['print_rand_state']]
    SAFE_RANDOM_FUNCTIONS_INITIALIZED = TRUE
}

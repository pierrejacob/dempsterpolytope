# dempsterpolytope - Gibbs sampler for Dempster's inference in Categorical distributions 
# Copyright (C) 2019 Pierre E. Jacob
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


# for a given phi_0, and assertion {phi < phi_0}
# upper probability is proportion of time random interval [phi_lower, phi_upper] intersects with [0, phi_0]
# lower probability is proportion of time random interval [phi_lower, phi_upper] is contained in [0, phi_0]
#'@export
linkage_cdf_lowerupper <- function(phi_0, lu_chain){
  # intersect if phi_lower < phi_0
  cdf_upper <- mean(apply(lu_chain, 1, function(v) v[1] < phi_0))
  # contains if phi_upper < phi_0
  cdf_lower <- mean(apply(lu_chain, 1, function(v) v[2] < phi_0))
  return(c(cdf_lower, cdf_upper))
}


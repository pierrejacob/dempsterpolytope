# dempsterpolytope - Gibbs sampler for Dempster's inference in Categorical distributions 
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

## implements the TV upper bounds of https://arxiv.org/abs/1905.09971
## based on meeting times
## from vector of meeting times, lag L and time t, obtain upper bound on TV
#'@export
tv_upper_bound <- function(meetingtimes, L, t){
  return(mean(pmax(0,ceiling((meetingtimes-L-t)/L))))
}

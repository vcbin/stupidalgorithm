#!/bin/bash - 
#===============================================================================
#
#          FILE:  gen_para.sh
# 
#         USAGE:  ./gen_para.sh 
# 
#   DESCRIPTION: generate specific AP experiment parameter file 
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR: YOUR NAME (), 
#       COMPANY: 
#       CREATED: 2011年01月21日 22时22分43秒 CST
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

algo_list="ede eda de-eda dmde"
func_list="rosenbrock foxhole langerman michaelwicz"

for func in $func_list
do
	for algo in $algo_list
	do
		common_file=${func}"_common_para"
		uniq_file=${func}"_"${algo}"_uniq_para"
		dst_file=${algo}"_"${func}"_hb.txt"
		cat $common_file $uniq_file>$dst_file
	done
done

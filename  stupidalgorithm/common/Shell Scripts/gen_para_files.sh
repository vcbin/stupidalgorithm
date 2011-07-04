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

algo_list="bde bbde sde my_sde spde jde"
func_list="F14 F15"

for func in $func_list
do
	for algo in $algo_list
	do
		common_file=${func}"_common_para"
		uniq_file=${algo}"_uniq_para"
		dst_file="parameters_"${algo}"_"${func}"_AP.txt"
		cat $common_file $uniq_file>$dst_file
	done
done


set key below center
set key box on


set title "frontier comparison"
plot 'mode_zdt3_best_pop.out' index 0 t 'mode' w p lt 19 pt 19 ps 1 , \
		'nsga2_zdt3_best_pop.out' index 0 t 'nsga2'w p pt 19 ps 1 # , \
		#'all_pop.out' index 1 t 'connecting line'w l lt 3
set xlabel "y1"
set ylabel "y2"
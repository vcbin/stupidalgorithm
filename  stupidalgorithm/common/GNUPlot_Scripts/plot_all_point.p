cd '.\Output'
set key below center
plot 'all_pop.out' index 0 t 'population' w p lt 19 pt 19 ps 1, \
		'all_pop.out' index 1 t 'archive'w p pt 19 ps 1, \
		'all_pop.out' index 1 t 'connecting line'w l lt 3
set xlabel "y1"
set ylabel "y2"
set title "all point"
pause 10
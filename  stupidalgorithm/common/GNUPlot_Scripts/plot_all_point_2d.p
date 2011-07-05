cd '.\Output'
set key below center
set key box on

if (!exists("gen")) print "gen is not defined"
if (!exists("arc_size")) print "arc_size is not defined"
if (!exists("stag_gen")) print "stag_gen is not defined"
if (!exists("max_gen")) print "max_gen is not defined"
set title "all points,current generation:".gen.",archive number:".arc_size. \
					",stagnation count=".stag_gen
plot 'all_pop.out' index 0 t 'population' w p pt 2 ps 1, \
		'all_pop.out' index 1 t 'archive'w p pt 1 ps 1 #, \
		# 'all_pop.out' index 1 t 'connecting line'w l lt 3
set xlabel "y1"
set ylabel "y2"
if ( !max_gen ) pause 5; else pause -1
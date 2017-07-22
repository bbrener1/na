./header.py $1 $2 mouse
./gold_standard.py $1 $3 mouse
./deviation_matrix.py $1 50
./threshold_analysis.py $1

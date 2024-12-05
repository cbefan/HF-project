# look upon my shitcode, ye mighty, and despair
#gcc hf.c -o hf.out -lm
i=10
rm -f data_${1}${2}.csv
while [ "$i" -le 500 ];
do
	r=`printf "%.2f\n" $(($i))e-2`
	#echo $r
	output=$(./hf.out $1 $2 $r)
	IN=$output
	IFS='	'
	for x in $IN
	do
		e_tot=$x
	done
	echo $r,$e_tot >> data_${1}${2}.csv
	i=$((i+1))
done

python make_plot.py data_${1}${2}.csv

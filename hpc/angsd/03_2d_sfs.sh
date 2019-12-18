module load angsd/0.928
module load parallel

export PATH=${PATH}:/usr/local/sci-irc/sw/angsd-0.928/angsd/misc/

do_2dsfs(){
	pair=$1
	p1=$(echo $pair | cut -f1 -d_)
	p2=$(echo $pair | cut -f2 -d_)

	echo "$p1 $p2"

	if [ ! -f ${p1}.${p2}.ml ];then
		realSFS ${p1}_af.saf.idx ${p2}_af.saf.idx -P 6 > ${p1}.${p2}.ml
	fi
}


pairs=$(
set -- DI FI PI PR MI
for a; do
    shift
    for b; do
        printf "%s_%s " "$a" "$b"
    done
done
)

export -f do_2dsfs

#parallel -j 10 do_2dsfs ::: ${pairs}

do_2dsfs "MI_north"

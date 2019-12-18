amodule load angsd/0.928
module load parallel

do_fst(){
	pair=$1
	p1=$(echo $pair | cut -f1 -d_)
	p2=$(echo $pair | cut -f2 -d_)

	echo "$p1 $p2"
	realSFS fst index ${p1}_af.saf.idx ${p2}_af.saf.idx -sfs ${p1}.${p2}.ml -fstout ${p1}.${p2}
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

export -f do_fst

parallel -j 10 do_fst ::: ${pairs}

do_fst "MI_north"
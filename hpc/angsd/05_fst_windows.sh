module load angsd/0.928
module load parallel

export PATH=${PATH}:/usr/local/sci-irc/sw/angsd-0.928/angsd/misc/


print_fst(){
	pair=$1
	p1=$(echo $pair | cut -f1 -d_)
	p2=$(echo $pair | cut -f2 -d_)

	echo "$p1 $p2"

	realSFS fst stats2 ${p1}.${p2}.fst.idx -win 50000 -step 10000 > ${p1}.${p2}.fst_windows.txt
	realSFS fst stats ${p1}.${p2}.fst.idx > ${p1}.${p2}.fst_stats.txt
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

export -f print_fst

#parallel -j 10 print_fst ::: ${pairs}

print_fst "MI_north"

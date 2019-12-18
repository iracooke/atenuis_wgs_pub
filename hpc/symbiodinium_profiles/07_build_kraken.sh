module load jellyfish/1.1.11
module load kraken

DBNAME=kraken_coral_symbionts


kraken-build --download-taxonomy --db $DBNAME

# Include standard bacteria in the database
#
kraken-build --download-library bacteria --db $DBNAME

find ./genome_krakenfa -name '*.fa' -print0 | \
    xargs -0 -I{} -n1 kraken-build --add-to-library {} --db $DBNAME

kraken-build --build --threads 16 --db $DBNAME


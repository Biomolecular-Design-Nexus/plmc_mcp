PLMC_DIR=./repo
DATASET=notebooks/example

# reformat a3m to a2m (a3m is the mmseqs2 output format)
reformat.pl a3m a2m  $DATASET/DHFR.a3m $DATASET/DHFR.a2m
python notebooks/rm_a2m_query_gaps.py $DATASET/DHFR.a2m $DATASET/alignment.a2m

$PLMC_DIR/plmc/bin/plmc \
    -o $DATASET/plmc/uniref100.model_params \
    -c $DATASET/plmc/uniref100.EC \
    -f $PROTEIN \
    -le 16.2 -lh 0.01 -m 200 -t 0.2 \
    -g $DATASET/alignment.a2m
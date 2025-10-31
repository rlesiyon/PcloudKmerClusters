SEQ_ID=(CP006059.1)
core_cutoffs=(5 10 20 40 )

for core_cutoff in "${core_cutoffs[@]}"
do
    poetry run python src/pcloudkmerclusters/cli.py \
        --seq_id CP006059.1 \
        --save_dir data \
        --core "$core_cutoff"
done


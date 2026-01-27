touch /etc/startup_script_launched

# TODO change workdir
cd $HOME_DIR
source mutatex_env/bin/activate

cd workdir


# try checking override dir
gsutil cp gs://$BATCH_PROCESSING_BUCKET/source_pdbs/$PDB_FILENAME .


# clean up old mutation PDBs and rotabases if current run has filled up memory
cleanUpMutations() {
    printf "Trying to clean up past mutations\n" &>> out.log
    # 1. check if mutations dir exists
    if [ -d mutations ]; then
        find mutations ! -name '*Repair.pdb' -name '*.pdb' -exec rm -f {} \; &>> out.log
        find mutations -name rotabase.txt -exec rm -f {} \; &>> out.log
    fi
}

touch /etc/mutatex_setup_succeeded

retry=0
while [[ ${retry:-} -le 3 ]]; do
    retry=$(( retry + 1 ))

    cleanUpMutations

    mutatex $PDB_FILENAME \
        -p $NUM_CORES \
        -m mutation_list.txt \
        -x $HOME_DIR/foldx/$FOLDX_BINARY_NAME \
        -f suite5 \
        -u \
        -R repair_runfile_template.txt \
        -M mutate_runfile_template.txt  &>> out.log
done

RESULT=${?:-}
if [[ ${RESULT:-} -ne 0 ]]; then
  printf "mutatex ERROR\n" > error.txt
  sudo poweroff; # stop this instance
  exit 1;
fi
touch /etc/mutatex_run_succeeded

# If command has succeeded, copy the right outputs
# to a cloud bucket and terminate the instance
find mutations -type f -name "*.pdb" -exec rm -f {} \;
find mutations -type f -name "*.txt" -exec rm -f {} \;
cd ..
# change this for any run name. you can change later names too
mv workdir ${GENE_NAME}_${UNIPROT_ID}_mutatex_${RUN_NAME}_results
tar -czvf ${GENE_NAME}_${UNIPROT_ID}_mutatex_${RUN_NAME}_results.tar.gz \
     ${GENE_NAME}_${UNIPROT_ID}_mutatex_${RUN_NAME}_results &&> out.log
gsutil cp ${GENE_NAME}_${UNIPROT_ID}_mutatex_${RUN_NAME}_results.tar.gz \
    gs://$BATCH_PROCESSING_BUCKET/mutatex_results/ &&> out.log

# stop the instance after job completion
sudo poweroff

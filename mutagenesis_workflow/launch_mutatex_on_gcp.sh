#!/bin/bash

export GENE_NAME=$1
export UNIPROT_ID=$2
export NUM_CORES=$3
export RUN_NAME=${4:-af2}

DRY_RUN=${5:-false} 
OVERRIDE_PDB_FILENAME=$6

# GCP settings
export GCP_PROJECT=$7
export GCP_REGION=${8:-us-central1}
export GCP_ZONE=${9:-us-central1-c}
export GCP_SERVICE_ACCOUNT=${10}
export GCP_MACHINE_TYPE=${11:-c3d-highcpu}

# Bucket/VM settings
export BATCH_PROCESSING_BUCKET=$12
export HOME_DIR=$13
export FOLDX_BINARY_NAME=$14

if [[ -z $GENE_NAME || -z $UNIPROT_ID || -z $NUM_CORES || -z $GCP_PROJECT || -z $GCP_SERVICE_ACCOUNT || -z $GCP_MACHINE_TYPE || -z $BATCH_PROCESSING_BUCKET || -z $HOME_DIR || -z $FOLDX_BINARY_NAME ]]; then
  printf "Usage:\n./launch_mutatex_on_gcp.sh [GENE_NAME] [UNIPROT_ID] [NUM_CORES] [RUNTYPE (default af2 else pdb)] [DRY_RUN (default false)] [OVERRIDE_PDB_FILENAME (default none)] [GCP_PROJECT] [GCP_REGION (default us-central1)] [GCP_ZONE (default us-central1-c)] [GCP_SERVICE_ACCOUNT] [GCP_MACHINE_TYPE (default c3d-highcpu)] [BATCH_PROCESSING_BUCKET] [HOME_DIR] [FOLDX_BINARY_NAME]\n"
  exit 1;
fi

gcloud config set project $GCP_PROJECT
gcloud config set compute/region $GCP_REGION
gcloud config set compute/zone $GCP_ZONE


INSTANCE_NAME=$(echo ${GENE_NAME}-mutatex-run | awk '{print tolower($0)}')
export PDB_FILENAME=${OVERRIDE_PDB_FILENAME:-"AF-${UNIPROT_ID}-F1-model_v4.pdb"}

printf "GENE_NAME: $GENE_NAME \nUNIPROT_ID: $UNIPROT_ID \nNUM_CORES: $NUM_CORES\n"
printf "RUN_NAME: $RUN_NAME \nINSTANCE_NAME: $INSTANCE_NAME \nPDB: $PDB_FILENAME\n"

envsubst < mutatex_template_run.sh > run.sh

if $DRY_RUN; then
  printf "Dry run: exiting now\n"
  exit 0;
fi

gcloud compute instances create "$INSTANCE_NAME" --project=$GCP_PROJECT \
    --zone=$GCP_ZONE --machine-type=$GCP_MACHINE_TYPE-$NUM_CORES \
    --network-interface=network-tier=PREMIUM,subnet=default --no-restart-on-failure \
    --maintenance-policy=TERMINATE --provisioning-model=SPOT \
    --instance-termination-action=STOP \
    --service-account=$GCP_SERVICE_ACCOUNT \
    --scopes=https://www.googleapis.com/auth/devstorage.read_write,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append \
    --create-disk=auto-delete=yes,boot=yes,device-name="$INSTANCE_NAME",mode=rw,size=100,image=projects/$GCP_PROJECT/global/images/mutatex-run-2025,type=projects/$GCP_PROJECT/zones/$GCP_ZONE/diskTypes/pd-balanced \
    --reservation-affinity=any --metadata-from-file startup-script=run.sh



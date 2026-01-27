## MutateX workflow (local + GCP)

### Setup: GCP Image Prep
1. Get a licensed FoldX binary by following the instructrions on the FoldX website (https://foldxsuite.crg.eu/).
2. Create a new VM in GCP.
3. Copy the FoldX binary to the GCP VM.
4. Clone mutateX to the GCP VM (https://github.com/ELELAB/mutatex) and install the requirements.txt in a virtual env called mutatex_env (pip install).
5. Make a directory called `workdir` and copy all of the foldx_cfg files into the directory.
5. Save a GCP disk image from the VM, and then delete the VM.

### Setup: Storage prep
1. Store any PDB structures needed in a bucket in GCP.
2. Create a path in a GCP bucket for results.

### Running

Run `launch_mutatex_on_gcp.sh` with the target gene and protein name, PDB structure name if relevent, and GCP project and the number of cores desired.



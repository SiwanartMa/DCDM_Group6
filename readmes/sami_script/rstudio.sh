#!/bin/bash -l

#SBATCH --job-name=ops-rstudio
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --signal=USR2
#SBATCH --cpus-per-task=1

module load rstudio/v2023.03.0_386-gcc-13.2.0-r-4.3.0-python-3.11.6

# get unused socket per https://unix.stackexchange.com/a/132524
export PASSWORD=$(openssl rand -base64 15)
readonly IPADDRESS=$(hostname -I | tr ' ' '\n' | grep '10.211.4.')
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:

   Linux and MacOS:
   ssh -NL 8787:${HOSTNAME}:${PORT} ${USER}@hpc.create.kcl.ac.uk

   Windows:
   ssh -m hmac-sha2-512 -NL 8787:${HOSTNAME}:${PORT} ${USER}@hpc.create.kcl.ac.uk

   and point your web browser to http://localhost:8787

2. Login to RStudio Server using the following credentials:

   user: ${USER}
   password: ${PASSWORD}

When done using the RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f ${SLURM_JOB_ID}

END

# Create custom database config
DBCONF=$TMPDIR/database.conf
if [ ! -e $DBCONF ]
then
printf "\nNOTE: creating $DBCONF database config file.\n\n"
echo "directory=$TMPDIR/var-rstudio-server" > $DBCONF
fi

rserver --server-user ${USER} --www-port ${PORT} --server-data-dir $TMPDIR/data-rstudio-server \
--secure-cookie-key-file $TMPDIR/data-rstudio-server/secure-cookie-key \
--database-config-file=$DBCONF --auth-none=0 \
--auth-pam-helper-path=pam-env-helper

printf 'RStudio Server exited' 1>&2

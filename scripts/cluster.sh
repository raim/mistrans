
## just a record of data and scripts for discovery

mkdir programs
cd programs

mkdir data
mkdir data/mammary

## from exon
## rsync -avz /home/raim/data/mammary/processedData r.machne@login.discovery.neu.edu:data/mammary/
## rsync -avz /home/raim/data/mammary/originalData r.machne@login.discovery.neu.edu:data/mammary/

## 
cd programs
git clone https://github.com/lh3/bioawk 
cd bioawk
make


cd programs
git clone https://github.com/psipred/s4pred
cd s4pred
wget http://bioinfadmin.cs.ucl.ac.uk/downloads/s4pred/weights.tar.gz
tar -xvzf weights.tar.gz


### TEST

echo -e ">test\nMQRISSLIHLSLFWAGVMSAIELVPEHQTVPVSIGVPATLRCSMKGEAIGNYYINWYRKTQGNTMTFIYREKDIYGPGFKDNFQGDIDIAKNLAVLKILAPSERDEGSYYCACDT" > test.fas
echo -e ">test\nMQRISSLIHLSLFWAGVMSAIELVPEHQTVPVSIGVPATLRCSMKGEAIGNYYINWYRKTQGNTMTFIYREKDIYGPGFKDNFQGDIDIAKNLAVLKILAPSERDEGSYYCACDTTTT" >> test.fas
echo -e ">test\nMQRISSLIHLSLFWAGVMSAIELVPEHQTVPVSIGVPATLRCSMKGEAIGNYYINWYRKTQGNTMTFIYREKDIYGPGFKDNFQGDIDIAKNLAVLKILAPSERDEGSYYCACDTGGG" >> test.fas
echo -e ">test\nMQRISSLIHLSLFWAGVMSAIELVPEHQTVPVSIGVPATLRCSMKGEAIGNYYINWYRKTQGNTMTFIYREKDIYGPGFKDNFQGDIDIAKNLAVLKILAPSERDEGSYYCACDTAAA" >> test.fas

module load python/3.6.6
srun -p short --nodes=1 --ntasks 1 python run_model.py test.fas > test.ss2

fasta=Homo_sapiens.GRCh38.pep.large.fa
awk "/^>/ {n++} (n>600 && n<701) {print}" $fasta > multi_test.fas
srun -p debug --nodes=10 --ntasks 1 python run_model.py test.fas --threads 10 > test.ss2

## time it; job took 60min on exon , w/o -threads option
##  StartTime=2024-01-26T09:10:48 EndTime=2024-01-26T09:30:48 Deadline=N/A
srun -p debug --nodes=1 --ntasks 1 --cpus-per-task 10 python run_model.py multi_test.fas --threads 10 > multi_test.ss2
scontrol show jobid 40492094
seff 40492094 ## Job Wall-clock time: 00:10:36
## testing many cpus
ncpu=20
srun -p short --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB python run_model.py multi_test.fas --threads $ncpu -z -s -o ./pred/
## 40492807 <- faster allocation with less memory requested
## Job Wall-clock time: 00:10:54
ncpu=2
srun -p short --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB python run_model.py multi_test.fas --threads $ncpu -z -s -o ./pred/
# 40493151
# 2024-01-26T10:29:05
ncpu=10
srun -p short --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=1.5GB python run_model.py multi_test.fas --threads $ncpu -z -s -o ./pred/
## 40493370
## StartTime=2024-01-26T11:02:06
## Job Wall-clock time: 00:21:53

## test largest 600 seq file
##rsync -avz s4pred_177.fas r.machne@login.discovery.neu.edu:programs/s4pred/
ncpu=20
srun -p short --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB python run_model.py s4pred_177.fas --threads $ncpu -z -s -o ./pred/
## 40493602
## StartTime=2024-01-26T11:36:00
## Job Wall-clock time: 00:38:1 !!
ncpu=40
srun -p short --nodes=1 --ntasks 1 --cpus-per-task $ncpu --mem=2GB python run_model.py s4pred_177.fas --threads $ncpu -z -s -o ./pred/
## 40494178
## StartTime=2024-01-26T12:27:35
## Job Wall-clock time: 00:46:46 - strangely no improvement!

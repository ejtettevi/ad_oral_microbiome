# download_sra_data.sh

# Install SRA Toolkit if not already installed
# conda install -c bioconda sra-tools

# Download all SRA runs listed in SRR_Acc_List
# Make sure SRR_Acc_List contains one SRR accession per line

while read run; do
    prefetch "$run"
    fasterq-dump "$run"
done < SRR_Acc_List.txt

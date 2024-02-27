# **** Variables ****
configfile: "config/cluster.yaml"

# **** Imports ****
import glob

# Cluster run template
#nohup snakemake --snakefile expanded_search.smk -j 5 --cluster "sbatch -t {cluster.time}" --cluster-config config/cluster.yaml --latency-wait 120 --use-conda &

# noinspection SmkAvoidTabWhitespace
rule all:
	input:
		# Generate BLAST database
		# expand("{run}/custom_db/{db_prefix}.phr",
		# 	run=config["run"],db_prefix=config["db_prefix"]),
		# PSI-Blast the MSA input using a list of sequence databases
		expand("{run}/psiblast_out_cluster{cluster}/{input_prefix}/db-{db_prefix}_query.blastout",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"], cluster=config["cluster"]),
		# Count TaxID occurrences in PSI-BLAST output
		expand("{run}/krona_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_taxid_counts.tsv",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"], cluster=config["cluster"]),
		# Export FASTA sequences of PSIBLAST Hits containing taxid labels
		expand("{run}/psiblast_out_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_hits.fasta",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"], cluster=config["cluster"]),
		# Generate Krona plot
		expand("{run}/krona_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_krona-plot.html",
			run=config["run"],db_prefix=config["db_prefix"], input_prefix=config["input_prefix"], cluster=config["cluster"]),
		# Re-align the initial MSA with the PSI-BLAST hit sequences
		expand("{run}/clustalo_cluster{cluster}/{input_prefix}/db-{db_prefix}_query.msa.fasta",
			run=config["run"], db_prefix=config["db_prefix"], input_prefix=config["input_prefix"], cluster=config["cluster"]),


rule iterative_search:
	input:
		fasta_in = lambda wildcards: glob.glob("{input_dir}/{input_prefix}.fasta".format(
			input_dir=config["input_dir"], input_prefix=wildcards.input_prefix)),
		# sequence_db = lambda wildcards: glob.glob("{shared_db_path}/{db_prefix}.01.phr".format(
		# 	shared_db_path=config["shared_db_path"], db_prefix=wildcards.db_prefix))
	output:
		psiblast_out = "{run}/psiblast_out_cluster{cluster}/{input_prefix}/db-{db_prefix}_query.blastout"
	params:
		db = config["shared_db_path"],
		taxid = config["taxid"],
		custom_cols = config["blast_custom_cols"],
	conda:
		"envs/blast.yaml"
	message:
		"""
Blasting query :\n {input.fasta_in}
Against database:\n {params.db}/{wildcards.db_prefix} \nGenerating:\n {output.psiblast_out}
Wildcards: {wildcards}
		"""
	shell:
		"""
        psiblast \
        -query {input.fasta_in} \
        -db {params.db}/{wildcards.db_prefix} \
        -outfmt "10 {params.custom_cols}" \
        -num_threads {threads} \
        -taxids {params.taxid} \
        -num_iterations 10 \
        -max_hsps 1 \
        -subject_besthit \
        -gapopen 9 \
        -inclusion_ethresh 1e-15 \
        -evalue 1e-10 \
        -qcov_hsp_perc 70 \
        -out {output.psiblast_out}
        """

# noinspection SmkAvoidTabWhitespace
rule taxid_parse:
	input:
		psiblast_out = "{run}/psiblast_out_cluster{cluster}/{input_prefix}/db-{db_prefix}_query.blastout"
	output:
		taxid_counts = "{run}/krona_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_taxid_counts.tsv",
		hits_fasta= "{run}/psiblast_out_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_hits.fasta"
	params:
		blast_col_names = config["blast_custom_cols"]
	conda:
		"envs/bio.yaml"
	script:
		"py/blastout_taxid_count.py"

# noinspection SmkAvoidTabWhitespace
rule krona:
	input:
		taxid_counts = "{run}/krona_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_taxid_counts.tsv"
	output:
		krona_chart = "{run}/krona_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_krona-plot.html",
	params:
		taxdump_path = config["taxdump_path"]
	conda:
		"envs/krona.yaml"
	message:
		"""
This rule implements Krona for interactive visualization of taxonomic distribution of the sequences.
The first five lines of the shell section are intended to correct Krona's issues with handling of taxonomic
background data.
Input data: {input.taxid_counts}
Output: {output.krona_chart}
Wildcards used in this rule: {wildcards}
		"""
	shell:
		"""		
mkdir $CONDA_PREFIX/bin/scripts || true
mkdir $CONDA_PREFIX/bin/taxonomy || true
ln -s $CONDA_PREFIX/opt/krona/scripts/extractTaxonomy.pl $CONDA_PREFIX/bin/scripts || true 
ln -s $CONDA_PREFIX/opt/krona/scripts/taxonomy.make $CONDA_PREFIX/bin/scripts || true

# Retry logic for ktUpdateTaxonomy.sh
retry_count=0
max_retries=10
while true; do
    if ktUpdateTaxonomy.sh; then
        break  # If the command succeeds, exit the loop
    else
        retry_count=$((retry_count+1))
        if [ $retry_count -ge $max_retries ]; then
            echo "Reached maximum retry attempts. Exiting."
            exit 1  # Exit script if maximum retries reached
        fi
        echo "Error occurred. Retrying ($retry_count of $max_retries)..."
        sleep 10  # Add a delay before retrying (adjust as needed)
    fi
done

ktImportTaxonomy -m 2 -t 1 -tax $CONDA_PREFIX/bin/taxonomy -o {output.krona_chart} {input.taxid_counts}
		"""

# noinspection SmkAvoidTabWhitespace
rule realignment:
	input:
		hits_fasta="{run}/psiblast_out_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_hits.fasta",
		fasta_in=lambda wildcards: glob.glob("{input_dir}/{input_prefix}.fasta".format(
			input_dir=config["input_dir"],input_prefix=wildcards.input_prefix))
	output:
		post_search_msa = "{run}/clustalo_cluster{cluster}/{input_prefix}/db-{db_prefix}_query.msa.fasta"
	params:
		merged_input = "{run}/clustalo_cluster{cluster}/{input_prefix}/db-{db_prefix}_query_merged-input.fasta",
	threads: config["threads"]
	shell:
		"""
		cat {input.fasta_in} {input.hits_fasta} > {params.merged_input}
		clustalo --iter 1 --threads {threads} -i {params.merged_input} -o {output.post_search_msa} -v
		"""

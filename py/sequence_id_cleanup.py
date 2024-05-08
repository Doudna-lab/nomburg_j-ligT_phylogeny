# == Native Modules
import re
# == Installed Modules
from Bio import SeqIO
# == Project Modules


def main():
	# Snakemake Imports
	#   SMK Inputs
	hits_fasta = str(snakemake.input.hits_fasta)
	#   SMK Outputs
	clean_hits_fasta = str(snakemake.output.clean_hits_fasta)
	#   SMK Params
	target_string_cleanup = str(snakemake.params.target_string_cleanup)


	#DEBUG
	# hits_fasta = "/Users/bellieny/projects/nomburg_j-ligT_phylogeny/dump/ligT_phyrec_nucleoside_transporter_sequences/psiblast_out_clusternucleoside_transporter_sequences/merged/merge_db-nr_query_hits.fasta"
	# target_string_cleanup = "partial"

	with open(hits_fasta) as fasta_in_handle:
		records = SeqIO.parse(fasta_in_handle, "fasta")
		clean_records = []
		for record in records:
			if not re.search(r"{}".format(target_string_cleanup), record.description, re.IGNORECASE) or re.search(r"{}".format(target_string_cleanup), record.id, re.IGNORECASE):
				clean_records.append(record)

	with open(clean_hits_fasta, "w") as out_clean_fasta:
		SeqIO.write(clean_records, out_clean_fasta, "fasta")


if __name__ == "__main__":
	main()

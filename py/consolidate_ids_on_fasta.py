# == Native Modules
import re
# == Installed Modules
from Bio import SeqIO
# == Project Modules


def main():
	# Snakemake Imports
	#   SMK Inputs
	cluster_representatives = str(snakemake.input.cluster_representatives)
	fasta_in_list = list(snakemake.input.fasta_in)
	#   SMK Outputs
	merged_clustered_input = str(snakemake.output.merged_clustered_input)

	#DEBUG
	# fasta_in_list = ['/Users/bellieny/projects/nomburg_j-ligT_phylogeny/input_data/sequences_for_pdes/56/sequence_reps/NS1_protein__NP_694863__Human_erythrovirus_V9__72197.fasta',
	# 				 '/Users/bellieny/projects/nomburg_j-ligT_phylogeny/input_data/sequences_for_pdes/56/sequence_reps/NS1__YP_009058894__Bufavirus-3__1391667.fasta',
	# 				 '/Users/bellieny/projects/nomburg_j-ligT_phylogeny/input_data/sequences_for_pdes/55/all_members/CI__YP_001552416__Tobacco_vein_banding_mosaic_virus__33765.fasta']
	#
	# cluster_representatives = '/Users/bellieny/projects/nomburg_j-ligT_phylogeny/dump/ligT_phyrec_cl56/psiblast_out_cluster56/merged/merge_db-nr_query_hits.fasta'

	post_merge_records = []
	external_records = {}
	# Process Clustered Sequence Representatives
	# Evaluate the need to spike original sequences back in the FASTA pool
	with open(cluster_representatives, "r") as fasta_handle:
		spike_back_behavior = []
		for ref_record in SeqIO.parse(fasta_handle, "fasta"):
			trimmed_ref_id = ref_record.id.split("|")[0]
			for fasta_file in fasta_in_list:
				out_record = SeqIO.read(fasta_file, "fasta")
				if out_record.id not in set(external_records.keys()):
					external_records.setdefault(out_record.id, out_record)
				if re.search(trimmed_ref_id, out_record.id):
					if out_record.id not in spike_back_behavior:
						spike_back_behavior.append(out_record.id)

			post_merge_records.append(ref_record)

		for out_record_id in external_records:
			if out_record_id not in spike_back_behavior:
				post_merge_records.append(external_records[out_record_id])

	with open(merged_clustered_input, "w") as processed_msa_handle:
		SeqIO.write(post_merge_records, processed_msa_handle, "fasta")


if __name__ == "__main__":
	main()

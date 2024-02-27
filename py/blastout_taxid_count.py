# Native modules
import re
import time
# Installed modules
import urllib.request
import urllib.error
import pandas as pd
import yaml
from Bio import Entrez
from Bio import SeqIO

# DEBUG INPUTS
# blastout_path = "/groups/doudna/projects/daniel_projects/boger_r/phyrec_screening/psiblast_out/db-phyrec_db_2023-08_query-test.out"
# config_path = "/groups/doudna/team_resources/toolbox/phyrec_screening/config/phyrec_processing.yaml"
# # Load config files
# import yaml
# with open(config_path, "r") as f:
# 	config = yaml.load(f, Loader=yaml.FullLoader)
# blast_col_names = config["blast_custom_cols"]


def ncbi_fetch(acc_list, ncbi_db, file_format):
	Entrez.email = "thedoudnalab@gmail.com"
	record_list = []
	max_retries = 50
	for acc in acc_list:

		# Attempt the search with retries
		for _ in range(max_retries):
			try:
				handle = Entrez.efetch(db=ncbi_db, id=f"{acc}", rettype=file_format, retmode="text")
				record = SeqIO.read(handle, file_format)
				handle.close()
				break  # Break the loop if successful
			except urllib.error.HTTPError as e:
				if e.code == 429:  # HTTP 429: Too Many Requests
					print(f"Received HTTP 429 error. Retrying in 10 seconds...")
					time.sleep(10)
				else:
					continue  # Re-raise other HTTP errors
			except urllib.error.URLError:
				continue

		# Synchronize accessions orthography
		try:
			sync_acc = re.search(acc, record.id)
			record.id = sync_acc.group()
			print(sync_acc, record.id)
			# Put it back on the record
			record_list.append(record)
		except (AttributeError, NameError):
			pass
	return record_list


def export_records(rec_list, output_path, file_format,):
	with open(f"{output_path}", "w") as handle:
		SeqIO.write(rec_list, handle, file_format)


def attach_label_to_fasta(seqrecords_list, id_to_label_df):
	new_seqrecords_list = []
	match_dict = {}
	id_to_label_df.apply(lambda row: match_dict.setdefault(row["sacc"], row["staxid"]), axis=1).to_dict()
	for record in seqrecords_list:
		try:
			record.id = f"{record.id}|{match_dict[record.id]}"
		except KeyError:
			pass
		new_seqrecords_list.append(record)
	return new_seqrecords_list


# # Load config file
# config_path = "../config/phylogeny_processing.yaml"
# with open(config_path, "r") as f:
# 	config = yaml.load(f, Loader=yaml.FullLoader)


def main():
	# Snakemake Imports
	#   SMK Inputs
	blastout_path = str(snakemake.input.psiblast_out)
	blast_col_names = str(snakemake.params.blast_col_names)
	#   SMK Outputs
	output_taxcount_table = str(snakemake.output.taxid_counts)
	output_fasta_hits = str(snakemake.output.hits_fasta)

	# Import blastout table
	blast_col_names_list = blast_col_names.split(" ")
	blast_df = pd.read_csv(blastout_path,
	                       names=blast_col_names_list,
	                       index_col=False).convert_dtypes().infer_objects()
	# Process pd Dataframe and take the unique set of hit IDs
	uniq_blast_hits = set(blast_df['sacc'].dropna().tolist())
	id_to_taxid = blast_df[blast_df['sacc'].isin(uniq_blast_hits)][['sacc', 'staxid']].drop_duplicates()

	# Count occurrences of each TaxID in the relevant column
	taxid_counts_df = id_to_taxid['staxid'].value_counts().reset_index()
	taxid_counts_df.columns = ['sacc', 'staxid_count']
	# Export output
	taxid_counts_df.to_csv(output_taxcount_table, sep="\t", index=False)

	# Parse hit sequences to FASTA file
	hit_seq_record_list = ncbi_fetch(uniq_blast_hits, 'protein', 'fasta')
	hit_taxid_seq_record_list = attach_label_to_fasta(hit_seq_record_list, id_to_taxid)
	export_records(hit_taxid_seq_record_list, output_fasta_hits, 'fasta')
	#
	# # Start the stopwatch / counter
	# t1_start = process_time()
	#
	# # Setup a support variable to control the IDs search across the taxid reference
	# processed_line = 0
	# line_report_threshold = 1000000
	# processed_ids_count = 0
	# t1_loop = 0
	# blasthit2taxid = {}
	# # Search taxid match file -> don't process any line unnecessarily
	# with open(taxid_match_path, 'r') as tax_handle:
	# 	for line in tax_handle:
	# 		processed_line += 1
	# 		if line.split("\t")[0] in unique_blast_hits:
	# 			processed_ids_count += 1
	# 			blasthit2taxid.setdefault(line.split("\t")[0], line.split("\t")[1])
	# 		if processed_ids_count == len(unique_blast_hits):
	# 			break
	# 		if processed_line >= line_report_threshold:
	# 			line_report_threshold += 1000000
	# 			t1_current = process_time()
	# 			print(f"Took {t1_current - t1_loop} to Process 1M lines. Now {processed_line} processed lines in total")
	# 			t1_loop = process_time()
	#
	#
	# # Stop the stopwatch / counter
	# t1_stop = process_time()
	#
	# print("Elapsed time:", t1_stop, t1_start)
	#
	# print("Elapsed time during the whole program in seconds:",
	#       t1_stop - t1_start)


if __name__ == "__main__":
	main()

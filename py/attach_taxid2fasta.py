# == Native Modules

# == Installed Modules
from Bio import SeqIO
import pandas as pd
import re
# == Project Modules


def extract_txid_from_fasta(fasta_seq_records_path):
	txid_list = []
	with open(fasta_seq_records_path, "r") as fasta_handle:
		for record in SeqIO.parse(fasta_handle, "fasta"):
			try:
				fasta_id_txid = int(record.id.split("|")[1].strip())
				txid_list.append(fasta_id_txid)
			except IndexError:
				continue
	return txid_list


def import_taxdump(taxdump_lineage_path):
	taxdump_dict = {}
	with open(taxdump_lineage_path, "r") as txid_lineage_handle:
		for line in txid_lineage_handle:
			key_txid = int(line.strip().split("\t|\t")[0])
			string_txid = line.strip().split("\t|\t")[1].strip("\t|").strip()
			associated_txids = string_txid.split(" ")
			taxdump_dict.setdefault(key_txid, associated_txids)
	return taxdump_dict


def add_txid_lineage(taxid_dict, taxdump_lineage_dict, txid_list):
	# reference_dict; txid_lineage
	for key_txid in set(taxdump_lineage_dict.keys()):
		if key_txid in set(txid_list):
			associated_txids = taxdump_lineage_dict[key_txid]
			anchor_txid = ''

			for associated_txid in set(associated_txids):
				associated_txid = int(associated_txid)
				if associated_txid in set(taxid_dict.keys()):
					anchor_txid = str(associated_txid)

			if anchor_txid:
				for txid_index in range(associated_txids.index(anchor_txid), len(associated_txids)):
					integeger_anchor = int(anchor_txid)
					taxid_dict.setdefault(associated_txids[txid_index], taxid_dict[integeger_anchor])
					taxid_dict.setdefault(key_txid, taxid_dict[integeger_anchor])

	return taxid_dict


def attach_taxid2fasta(fasta_seq_records_path, taxid_to_gentype_dictionary, full_txid_dictionary):
	loop_fasta_seq_records = []
	count = 0
	with open(fasta_seq_records_path, "r") as fasta_handle:
		for record in SeqIO.parse(fasta_handle, "fasta"):
			try:
				fasta_id_txid = int(record.id.split("|")[1].strip())
				fasta_seq_id = record.id.split("|")[0].strip()
			except IndexError:
				loop_fasta_seq_records.append(record)
				continue
			try:
				txid_lineage_string = "|".join(map(str, full_txid_dictionary[fasta_id_txid]))
				record.id = f"{fasta_seq_id}|{txid_lineage_string}|{fasta_id_txid}"
			except KeyError:
				loop_fasta_seq_records.append(record)
				continue
			if fasta_id_txid in set(taxid_to_gentype_dictionary.keys()):
				record.id = f"{record.id}|{taxid_to_gentype_dictionary[fasta_id_txid]}"
				count += 1
				print(f"{count}: {record.id}")
			loop_fasta_seq_records.append(record)
		print(f"Found {count}")
	return loop_fasta_seq_records


def remove_fasta_duplicates(fasta_seq_records):
	unique_records = []
	check_duplicates = []
	for record in fasta_seq_records:
		if record.id not in set(check_duplicates):
			unique_records.append(record)
			check_duplicates.append(record.id)
	return unique_records


def filter_txid_lineage(fasta_seq_records, txid_lineage_filter, keep_sequences):
	filtered_seq_records = []
	for record in fasta_seq_records:
		if re.search(r"\|{}\|".format(txid_lineage_filter), record.id):
			filtered_seq_records.append(record)
		elif record.id in set(keep_sequences):
			filtered_seq_records.append(record)
	return filtered_seq_records


def main():
	# Snakemake Imports
	#   SMK Inputs
	id_match_table = str(snakemake.input.id_match_table)
	merged_msa_fasta = str(snakemake.input.merged_msa_fasta)
	txid_lineage = str(snakemake.input.txid_lineage)
	#   SMK Params
	match_colum = str(snakemake.params.match_colum)
	add_colum = str(snakemake.params.add_colum)
	filter_taxid = str(snakemake.params.filter_taxid)
	reference_prefixes = list(snakemake.params.reference_prefixes)
	#   SMK Outputs
	merged_msa_fasta_txid = str(snakemake.output.merged_msa_fasta_txid)

	#DEBUG INPUT
	# id_match_table = "/Users/bellieny/projects/nomburg_j-ligT_phylogeny/input_data/family_genome_types_taxID.tsv"
	# merged_msa_fasta = "/Users/bellieny/projects/nomburg_j-ligT_phylogeny/dump/ligT_phyrec/clustalo_clusterligT/merged/db-nr_query_merged-input.msa.fasta"
	# txid_lineage = "/Users/bellieny/projects/nomburg_j-ligT_phylogeny/input_data/taxidlineage.dmp"
	# match_colum = "family_taxonID"
	# add_colum = "family_genome_type"

	# Import Reference Table
	df_reference = pd.read_csv(id_match_table, sep="\t")
	reference_dict = {}
	df_reference.apply(
		lambda row: reference_dict.setdefault(
			row.loc[match_colum], row.loc[add_colum]),
		axis=1
	)

	txids_in_fasta = extract_txid_from_fasta(merged_msa_fasta)

	taxdump_dictionary = import_taxdump(txid_lineage)

	full_lineage_dict = add_txid_lineage(reference_dict, taxdump_dictionary, txids_in_fasta)

	txid_records = attach_taxid2fasta(merged_msa_fasta, full_lineage_dict, taxdump_dictionary)

	unique_txid_records = remove_fasta_duplicates(txid_records)

	filtered_txid_records = unique_txid_records

	if filter_taxid != "false":
		filter_taxid = int(filter_taxid)
		filtered_txid_records = filter_txid_lineage(unique_txid_records, filter_taxid, reference_prefixes)
		print(f"Number of records after filtering taxid: {len(filtered_txid_records)}")

	with open(merged_msa_fasta_txid, "w") as processed_msa_handle:
		SeqIO.write(filtered_txid_records, processed_msa_handle, "fasta")


if __name__ == "__main__":
	main()

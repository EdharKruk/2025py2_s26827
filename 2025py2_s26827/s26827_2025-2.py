#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Basic script to connect to NCBI and retrieve genetic sequence records for a given taxonomic ID.
"""

from Bio import Entrez, SeqIO
from io import StringIO
import os
import time
import pandas as pd
import matplotlib.pyplot as plt



class NCBIRetriever:
    def __init__(self, email, api_key):
        self.email = email
        self.api_key = api_key

        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid):
        print(f"Searching for records with taxID: {taxid}")
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")

            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print(f"No records found for {organism_name}")
                return None

            print(f"Found {count} records")

            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count

            return count

        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, start=0, max_records=10):
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            batch_size = min(max_records, 500)
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            records_text = handle.read()
            handle.close()

            gbio = StringIO(records_text)
            parsed_records = list(SeqIO.parse(gbio, "genbank"))
            return parsed_records

        except Exception as e:
            print(f"Error fetching records: {e}")
            return []

def main():
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")

    retriever = NCBIRetriever(email, api_key)

    taxid = input("Enter taxonomic ID (taxid) of the organism: ")
    min_length = int(input("Enter minimum sequence length: "))
    max_length = int(input("Enter maximum sequence length: "))

    count = retriever.search_taxid(taxid)
    if not count:
        print("No records found. Exiting.")
        return

    print("\nFetching sample records...")
    sample_records = retriever.fetch_records(start=0, max_records=20)

    filtered_records = [
        rec for rec in sample_records
        if min_length <= len(rec.seq) <= max_length
    ]

    print(f"Filtered to {len(filtered_records)} records within length range.")

    output_file = f"taxid_{taxid}_sample.gb"
    with open(output_file, "w") as f:
        SeqIO.write(filtered_records, f, "genbank")

    print(f"Saved filtered sample records to {output_file}")

    csv_data = []
    for record in filtered_records:
        csv_data.append({
            "Accession": record.id,
            "Length": len(record.seq),
            "Description": record.description
        })

    df = pd.DataFrame(csv_data)
    csv_file = f"taxid_{taxid}_report.csv"
    df.to_csv(csv_file, index=False)
    print(f"Saved CSV report to {csv_file}")
    

    df_sorted = df.sort_values(by="Length", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.plot(df_sorted["Accession"], df_sorted["Length"], marker='o')
    plt.xticks(rotation=90)
    plt.xlabel("Accession Number")
    plt.ylabel("Sequence Length")
    plt.title(f"Sequence Lengths for TaxID {taxid}")
    plt.tight_layout()

    chart_file = f"taxid_{taxid}_lengths.png"
    plt.savefig(chart_file)
    plt.close()

    print(f"Saved chart to {chart_file}")


    
if __name__ == "__main__":
    main()

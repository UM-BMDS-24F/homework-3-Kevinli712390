from Bio.Blast.Applications import NcbiblastpCommandline, NcbimakeblastdbCommandline
from Bio import SeqIO
import os

# Paths
human_fasta = "human.fa"
mouse_fasta = "mouse.fa"
output_file = "blast_results.txt"

# 1. Create a BLAST database for the mouse sequences
makeblastdb_cline = NcbimakeblastdbCommandline(dbtype="prot", input_file=mouse_fasta)
stdout, stderr = makeblastdb_cline()

# 2. Set up BLAST search for each human sequence
with open(output_file, "w") as out_handle:
    for record in SeqIO.parse(human_fasta, "fasta"):
        # Save human sequence to a temporary file
        with open("temp_query.fa", "w") as temp_query:
            SeqIO.write(record, temp_query, "fasta")
        
        # Run BLAST
        blastp_cline = NcbiblastpCommandline(
            query="temp_query.fa",
            db=mouse_fasta,
            evalue=0.001,
            outfmt=5,
            out="temp_result.xml"
        )
        stdout, stderr = blastp_cline()
        
        # Parse BLAST results (assuming XML output)
        from Bio.Blast import NCBIXML
        with open("temp_result.xml") as result_handle:
            blast_record = NCBIXML.read(result_handle)
            
            # Get the top hit
            if blast_record.alignments:
                best_alignment = blast_record.alignments[0]
                best_hsp = best_alignment.hsps[0]
                out_handle.write(f"Human ID: {record.id}\n")
                out_handle.write(f"Mouse ID: {best_alignment.hit_id}\n")
                out_handle.write(f"Alignment:\n{best_hsp.align_length}\n")
                out_handle.write(f"E-value: {best_hsp.expect}\n")
                out_handle.write(f"Bitscore: {best_hsp.bits}\n\n")

# Clean up temporary files
os.remove("temp_query.fa")
os.remove("temp_result.xml")

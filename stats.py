import time
from Bio.Blast import NCBIXML, NCBIWWW
from config import DATA_DIR
import os
from Bio import SeqIO
def truncate_sequence(seq, max_len=60):
    if len(seq) > max_len:
        part = max_len // 3
        return seq[:part] + '...' + seq[-part:]
    return seq
def run_blast_to_xml(query_seq, output_xml):

    try:
        print("Running BLASTp on NCBI... (this may take 2-5 minutes)")
        result_handle = NCBIWWW.qblast(
            program="blastp",
            database="swissprot",  # Smaller database for faster results
            sequence=query_seq,
            matrix_name="BLOSUM62",
            gapcosts="11 1",  # Match your aligner's gap penalties
            descriptions=100,
            alignments=100,
            hitlist_size=10
        )
        
        with open(output_xml, "w") as out_file:
            out_file.write(result_handle.read())
        print(f"BLAST results saved to {output_xml}")
        
    except Exception as e:
        print(f"BLAST failed: {e}. Retrying in 30 seconds...")
        time.sleep(30)
        run_blast_to_xml(query_seq, output_xml)  # Simple retry
def generate_blast_report(blast_xml_path, output_path):
    with open(blast_xml_path) as in_file:
        blast_records = list(NCBIXML.parse(in_file))
    
    if not blast_records:
        print("No BLAST records found.")
        return
    
    blast_record = blast_records[0]
    
    # Extract header information
    query_id = blast_record.query_id
    query_desc = blast_record.query
    query_len = blast_record.query_length
    database = blast_record.database
    
    # Extract parameters correctly
    params = getattr(blast_record, 'params', None)
    matrix = getattr(params, 'matrix', 'N/A') if params else 'N/A'
    expect = getattr(params, 'expect', 'N/A') if params else 'N/A'
    gap_open = getattr(params, 'gap_open', 'N/A') if params else 'N/A'  
    gap_extend = getattr(params, 'gap_extend', 'N/A') if params else 'N/A'
    
    # Prepare hits data with domain annotations
    hits_data = []
    for alignment in blast_record.alignments:
        hit_accession = alignment.accession
        hit_def = alignment.hit_def
        
        # Extract species and domains
        scientific_name = hit_def.split('[')[-1].split(']')[0] if '[' in hit_def else 'N/A'
        functional_notes = []
        if 'RING' in hit_def.upper():
            functional_notes.append('RING domain')
        if 'BRCT' in hit_def.upper():
            functional_notes.append('BRCT domain')
        domain_notes = ', '.join(functional_notes) if functional_notes else 'None'
        
        hsp = alignment.hsps[0]
        percent_identity = (hsp.identities / hsp.align_length) * 100
        
        hits_data.append({
            'accession': hit_accession,
            'description': hit_def.split(' [')[0],  
            'scientific_name': scientific_name,
            'bit_score': hsp.bits,
            'e_value': hsp.expect,
            'percent_identity': percent_identity,
            'align_length': hsp.align_length,
            'q_seq': truncate_sequence(hsp.query),
            'midline': truncate_sequence(hsp.match),
            'h_seq': truncate_sequence(hsp.sbjct),
            'q_start': hsp.query_start,
            'q_end': hsp.query_end,
            'h_start': hsp.sbjct_start,
            'h_end': hsp.sbjct_end,
            'domains': domain_notes
        })
    
    # Write report with improved formatting
    with open(output_path, 'w') as out_file:
        # Header
        out_file.write(f"BLAST Report\n")
        out_file.write(f"{'='*80}\n")
        out_file.write(f"Query ID: {query_id}\n")
        out_file.write(f"Description: {query_desc}\n")
        out_file.write(f"Length: {query_len} residues\n")
        out_file.write(f"Database: {database}\n\n")
        
        out_file.write(f"Parameters:\n")
        out_file.write(f"Matrix: {matrix}\n")
        out_file.write(f"Expected threshold: {expect}\n")
        out_file.write(f"Gap costs: open={gap_open}, extend={gap_extend}\n\n")
        
        # Hits table (wider columns)
        out_file.write(f"Significant hits ({len(hits_data)} found):\n")
        out_file.write("-" * 120 + "\n")
        out_file.write("No.  Accession   Description                                   Scientific Name              Bit Score    E-value    % Identity  Align Length\n")
        out_file.write("-" * 120 + "\n")
        
        for idx, hit in enumerate(hits_data, 1):
            out_file.write(f"{idx:<4} {hit['accession']:<12} {hit['description'][:50]:<50} {hit['scientific_name'][:25]:<25} "
                           f"{hit['bit_score']:>10.1f}  {hit['e_value']:>10.2g}  {hit['percent_identity']:>10.1f}%  {hit['align_length']:>12}\n")
        out_file.write("-" * 120 + "\n\n")
        
        # Alignment details with domain annotations
        for idx, hit in enumerate(hits_data, 1):
            out_file.write(f"Hit {idx}: {hit['accession']} ({hit['scientific_name']})\n")
            out_file.write(f"Bit score: {hit['bit_score']:.1f}, E-value: {hit['e_value']:.2g}, % Identity: {hit['percent_identity']:.1f}%\n")
            out_file.write(f"Conserved domains: {hit['domains']}\n")
            out_file.write("Alignment:\n")
            out_file.write(f"Query: {hit['q_start']} {hit['q_seq']} {hit['q_end']}\n")
            out_file.write(f"        {hit['midline']}\n")
            out_file.write(f"Hit:   {hit['h_start']} {hit['h_seq']} {hit['h_end']}\n\n")
        
        # Footer
        out_file.write(f"{'='*80}\n")
        out_file.write(f"Total hits found: {len(hits_data)}\n")
        out_file.write(f"Parameters used: matrix={matrix}, expect={expect}, gap_open={gap_open}, gap_extend={gap_extend}\n")


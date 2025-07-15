import argparse
import logging
import os.path
import time

from BuildPhIPSeqLibrary.config import OUTPUT_DIR, BARCODED_NUC_FILE, ORDER_FILE # done
from BuildPhIPSeqLibrary.construct_nucleotide_sequences import aa_to_nuc, get_edge_restrictions
from BuildPhIPSeqLibrary.output_to_order import transfer_to_order # done
from BuildPhIPSeqLibrary.read_input_files import get_input_files, read_file # done
from BuildPhIPSeqLibrary.sequence_ids import add_sequences_to_files_list # done
from BuildPhIPSeqLibrary.split_sequences_to_oligos import split_and_map_new_sequences # done

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Building new library main.")
    parser.add_argument('--overwrite', dest='overwrite', default=False, action='store_true')
    args = parser.parse_args()
    if args.overwrite and os.path.exists(OUTPUT_DIR):
        assert len([filename for filename in os.listdir(OUTPUT_DIR) if
                    filename != 'README.md' and not filename.startswith(
                        '.')]) == 0, f"""In order to overwrite you must empty the output dir {OUTPUT_DIR}
Consider running:
for filename in os.listdir('{OUTPUT_DIR}'):
    if filename != 'README.md':
        os.remove(os.path.join('{OUTPUT_DIR}', filename))""""""
        """
    
    # Make sure Output dir is exist, if not, create it.    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    get_edge_restrictions()     # From construct_nucleotide_sequences.py
    
    files = get_input_files()   # Get input files
    logging.basicConfig(level=logging.INFO)     # Enables the logging.INFO statement to output log information.
    
    # Start to iterate input files
    for filename in files:
        # os.path.basename(filename): Obtains the filename (excluding the path).
        # os.path.splitext(...) [0]: Remove the file extension, keep only the main file name, 
        # and store it in the file_end variable.
        file_end = os.path.splitext(os.path.basename(filename))[0]
        
        # Record the file name and time of the current process to facilitate debugging and tracking progress.
        logging.info(f"Working on file {file_end}, {time.ctime()}")
        logging.info(f"{file_end}: reading file, {time.ctime()}")
        
        # Returns a dictionary with the sequence ID as the key and the corresponding sequence as the value. 
        # dictionary with sequence_ID -> AA_sequence
        seq_id_to_sequences = read_file(filename)   # From read_input_files.py
        
        # Record the number of sequences read and begin to recognize new sequences.
        logging.info(f"{file_end}: Got {len(seq_id_to_sequences)} sequences, {time.ctime()}")
        logging.info(f"{file_end}: Identifying new sequences, {time.ctime()}")
        
        # This file adds a list of sequences from a newly read file to the sequences table.
        # It will return only newly-added sequences.
        seq_id_to_sequences = add_sequences_to_files_list(seq_id_to_sequences, filename)    # From sequence_ids.py
        
        # Records the number of sequences processed by add_sequences_to_files_list().
        logging.info(f"{file_end}: Got {len(seq_id_to_sequences)} new sequences, {time.ctime()}")
        
        # If there are no new sequences (seq_id_to_sequences is empty), log and go to next file.
        if len(seq_id_to_sequences) == 0:
            logging.info(f"{file_end}: Contributed no new sequences.")
            continue
        
        # Records are undergoing oligonucleotide (oligos) conversion.
        logging.info(f"{file_end}: Converting sequences to oligos, {time.ctime()}")
        all_oligos_aa_sequences, new_oligos_aa_sequences = split_and_map_new_sequences(seq_id_to_sequences) # From split_sequences_to_oligos.py
        
        # Record the number of oligos after conversion: new + all oligos
        logging.info(f"{file_end}: "
              f"Converted {len(new_oligos_aa_sequences)} new oligos "
              f"overall {len(all_oligos_aa_sequences)} oligos, {time.ctime()}")
        
        # Record barcoding oligos from oligos_sequence.csv to barcoded_nuc_file.csv
        # AA -> DNA nucleotide barcodes
        logging.info(f"{file_end}: Barcoding oligos, {time.ctime()}")
        oligo_barcoded_sequences = aa_to_nuc(new_oligos_aa_sequences) # From construct_nucleotide_sequences.py
        
        logging.info(
            f"{file_end}: Finished barcoding oligos. "
            f"Current number of oligos is {len(oligo_barcoded_sequences)}, {time.ctime()}")
        
    # 
    if len(files) > 0:
        logging.info(f"Finished creating sequences. Find them in {BARCODED_NUC_FILE}, {time.ctime()}")
        logging.info(f"Converting all sequences to order, {time.ctime()}")
        # Generate an order file 
        transfer_to_order(oligo_barcoded_sequences) # From output_to_order.py
        logging.info(f"Find file of order in {ORDER_FILE}, {time.ctime()}")
    else:
        logging.info("No new files")

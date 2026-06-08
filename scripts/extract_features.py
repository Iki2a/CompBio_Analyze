import os
import zipfile
import tempfile
import shutil
import pandas as pd
from Bio import SeqIO
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DATA_DIR = PROJECT_ROOT / "data" / "processed"
GENOMES_DIR = PROJECT_ROOT / "data" / "genomes"

def get_gc(seq):
    """Menghitung persentase GC secara mandiri tanpa bergantung pada versi Biopython"""
    if not seq:
        return None
    s = str(seq).upper()
    total = len(s)
    if total == 0:
        return None
    gc_count = s.count('G') + s.count('C')
    return (gc_count / total) * 100.0

def extract_features():
    selected_csv = PROCESSED_DATA_DIR / 'selected_genomes_for_download.csv'
    df = pd.read_csv(selected_csv)
    
    results = []
    
    for idx, row in df.iterrows():
        acc = str(row['Assembly_or_accession']).strip()
        # Mengembalikan nama spesies asli seperti di CSV TEMPURA
        species_original = row['Genus_and_species']
        species_safe = str(row['Genus_and_species']).replace(' ', '_').replace('"', '')
        topt_col = [c for c in df.columns if 'Topt_average' in c][0]
        topt = row[topt_col]
        
        # Pengelompokan
        if topt > 45:
            group = 'Termofil'
        else:
            group = 'Mesofil'
            
        zip_path = GENOMES_DIR / f"{acc}_{species_safe}.zip"
        
        if not zip_path.exists():
            print(f"[{idx+1}/{len(df)}] File tidak ditemukan: {zip_path}, melewati...")
            continue
            
        print(f"[{idx+1}/{len(df)}] Memproses {acc} ({species_original})...")
        
        temp_dir = tempfile.mkdtemp()
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
                
            # Mencari file FNA dan GFF3
            fna_file = None
            gff_file = None
            
            for root, dirs, files in os.walk(temp_dir):
                for file in files:
                    if file.endswith('.fna'):
                        fna_file = os.path.join(root, file)
                    if file.endswith('.gff'):
                        gff_file = os.path.join(root, file)
                        
            if not fna_file or not gff_file:
                print(f"  -> Melewati {acc}: File FASTA atau GFF tidak lengkap.")
                continue
                
            # Membaca sequence FASTA
            seq_dict = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))
            
            # Hitung GCw (Whole Genome)
            total_gc = 0
            total_len = 0
            for seq_record in seq_dict.values():
                seq_str = str(seq_record.seq).upper()
                total_len += len(seq_str)
                total_gc += seq_str.count('G') + seq_str.count('C')
            
            gc_w = (total_gc / total_len * 100) if total_len > 0 else None
            
            # Membaca Anotasi GFF3 untuk tRNA, 16S, 23S
            trna_seqs = ""
            rrna_16s_seqs = ""
            rrna_23s_seqs = ""
            
            with open(gff_file, 'r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('#'): continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9: continue
                    
                    seqid = parts[0]
                    feat_type = parts[2]
                    # Format GFF menggunakan 1-based index
                    start = int(parts[3]) - 1
                    end = int(parts[4])
                    strand = parts[6]
                    attributes = parts[8]
                    
                    if seqid not in seq_dict: continue
                    
                    if feat_type == 'tRNA':
                        seq_obj = seq_dict[seqid].seq[start:end]
                        if strand == '-': seq_obj = seq_obj.reverse_complement()
                        trna_seqs += str(seq_obj)
                        
                    elif feat_type == 'rRNA':
                        # Atribut GFF RefSeq biasanya menulis nama produk seperti '16S ribosomal RNA'
                        if '16S' in attributes:
                            seq_obj = seq_dict[seqid].seq[start:end]
                            if strand == '-': seq_obj = seq_obj.reverse_complement()
                            rrna_16s_seqs += str(seq_obj)
                        elif '23S' in attributes:
                            seq_obj = seq_dict[seqid].seq[start:end]
                            if strand == '-': seq_obj = seq_obj.reverse_complement()
                            rrna_23s_seqs += str(seq_obj)
                            
            gc_trna = get_gc(trna_seqs)
            gc_16s = get_gc(rrna_16s_seqs)
            gc_23s = get_gc(rrna_23s_seqs)
            
            results.append({
                'species': species_original,
                'Topt': topt,
                'group': group,
                'GCw': round(gc_w, 3) if gc_w is not None else None,
                'GC_tRNA': round(gc_trna, 3) if gc_trna is not None else None,
                'GC_16S': round(gc_16s, 3) if gc_16s is not None else None,
                'GC_23S': round(gc_23s, 3) if gc_23s is not None else None
            })
            
        finally:
            shutil.rmtree(temp_dir)
            
    out_df = pd.DataFrame(results)
    out_csv = PROCESSED_DATA_DIR / 'real_gc_data.csv'
    out_df.to_csv(out_csv, index=False)
    print(f"\nSelesai! Fitur GC berhasil diekstrak dan disimpan ke '{out_csv}'")

if __name__ == '__main__':
    extract_features()

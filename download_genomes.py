import pandas as pd
import os
import subprocess
import zipfile
import shutil

def main():
    # Buat direktori tujuan
    out_dir = "dataset_genom"
    os.makedirs(out_dir, exist_ok=True)
    
    # Baca CSV TEMPURA
    csv_file = '200617_TEMPURA.csv'
    print(f"Membaca {csv_file}...")
    
    # Memastikan encoding terbaca dengan baik
    try:
        df = pd.read_csv(csv_file, encoding='utf-8')
    except UnicodeDecodeError:
        df = pd.read_csv(csv_file, encoding='latin1')
        
    # Cari kolom Topt_average
    topt_col = [c for c in df.columns if 'Topt_average' in c][0]
    
    # Cleaning data: hanya ambil yang Bacteria dan punya Accession GCA/GCF
    df['Assembly_or_accession'] = df['Assembly_or_accession'].fillna('')
    df_valid = df[
        (df['Superkingdom'] == 'Bacteria') & 
        (df['Assembly_or_accession'].str.startswith('GC'))
    ].copy()
    
    # Pastikan Topt adalah numerik
    df_valid[topt_col] = pd.to_numeric(df_valid[topt_col], errors='coerce')
    df_valid = df_valid.dropna(subset=[topt_col])
    
    # Filter kelompok
    df_thermo = df_valid[df_valid[topt_col] > 45]
    df_meso = df_valid[(df_valid[topt_col] >= 15) & (df_valid[topt_col] <= 45)]
    
    print(f"Total Termofilik valid: {len(df_thermo)}")
    print(f"Total Mesofilik valid: {len(df_meso)}")
    
    # Ambil sampel (54 Termofilik, 37 Mesofilik)
    n_thermo = min(54, len(df_thermo))
    n_meso = min(37, len(df_meso))
    
    df_thermo_sample = df_thermo.sample(n=n_thermo, random_state=42)
    df_meso_sample = df_meso.sample(n=n_meso, random_state=42)
    
    final_df = pd.concat([df_thermo_sample, df_meso_sample])
    
    # Simpan daftar yang akan diunduh ke CSV baru agar mudah dilacak
    final_df.to_csv("selected_genomes_for_download.csv", index=False)
    print(f"\nDisimpan daftar 91 genom ke 'selected_genomes_for_download.csv'")
    
    # Download dengan datasets CLI
    print("\nMulai mengunduh FASTA dan GFF3 menggunakan NCBI Datasets CLI...")
    success_count = 0
    
    for idx, row in final_df.iterrows():
        acc = str(row['Assembly_or_accession']).strip()
        species = str(row['Genus_and_species']).replace(' ', '_').replace('"', '')
        zip_filename = os.path.join(out_dir, f"{acc}_{species}.zip")
        
        if os.path.exists(zip_filename):
            print(f"[{success_count+1}/91] Melewati {acc} ({species}) - sudah ada...")
            success_count += 1
            continue

        print(f"[{success_count+1}/91] Mengunduh {acc} ({species})...")
        
        cmd = [
            ".\\datasets.exe", "download", "genome", "accession", acc,
            "--include", "genome,gff3",
            "--filename", zip_filename
        ]
        
        # Eksekusi CLI
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if result.returncode == 0:
                success_count += 1
                print(f"  -> Sukses diunduh ke {zip_filename}")
            else:
                print(f"  -> Gagal mengunduh {acc}: {result.stderr.strip()}")
        except Exception as e:
            print(f"  -> Terjadi kesalahan saat memanggil CLI: {e}")
            
    print(f"\nSelesai! {success_count} genom berhasil diunduh dan disimpan di dalam folder '{out_dir}/'")

if __name__ == '__main__':
    main()

# Analisis GC Content dan Suhu Optimal Bakteri

Proyek ini menganalisis hubungan antara GC content dan suhu pertumbuhan optimal
(`Topt`) pada bakteri termofilik dan mesofilik. Data suhu berasal dari TEMPURA,
sedangkan data genom/anotasi berasal dari NCBI Datasets.

Repo ini sengaja menyimpan CSV, ZIP genom, figure PNG, dan `datasets.exe` agar
analisis bisa dijalankan ulang secara praktis tanpa perlu download ulang semua
data.

## Struktur Proyek

- `data/raw/200617_TEMPURA.csv`: dataset TEMPURA mentah.
- `data/processed/selected_genomes_for_download.csv`: 91 genom yang dipilih dari TEMPURA.
- `data/processed/real_gc_data.csv`: fitur GC hasil ekstraksi dari genom lokal.
- `data/external/gct/`: dataset eksternal `gcT` dari paket R `seqinr`.
- `data/genomes/`: ZIP genom lokal hasil download dari NCBI Datasets.
- `scripts/download_genomes.py`: memilih sampel dan mengunduh genom/anotasi.
- `scripts/extract_features.py`: mengekstrak `GCw`, `GC_tRNA`, `GC_16S`, dan `GC_23S`.
- `scripts/gc_analysis.py`: menjalankan statistik, model machine learning, dan membuat figure.
- `scripts/import_gct_dataset.py`: mengimpor dataset eksternal `gcT` tanpa runtime R.
- `scripts/gct_validation.py`: validasi eksternal GCw vs `Topt` menggunakan `gcT`.
- `notebooks/gc_analysis.ipynb`: versi notebook dari analisis.
- `figures/`: output visualisasi utama.
- `tools/datasets.exe`: NCBI Datasets CLI lokal.
- `docs/gemini.md`: catatan/ringkasan proyek.
- `docs/gct_external_dataset.md`: penjelasan dataset eksternal `gcT`.

## Setup

Gunakan Python 3.10+.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Jika memakai Jupyter:

```powershell
jupyter notebook notebooks/gc_analysis.ipynb
```

## Cara Menjalankan Pipeline

Untuk menjalankan analisis dari data yang sudah tersedia:

```powershell
python scripts/gc_analysis.py
```

Untuk mengekstrak ulang fitur GC dari ZIP genom lokal:

```powershell
python scripts/extract_features.py
```

Untuk memilih sampel dan mengunduh ulang genom memakai NCBI Datasets CLI:

```powershell
python scripts/download_genomes.py
```

`scripts/download_genomes.py` memakai `tools/datasets.exe`.

## Dataset Validasi Eksternal

Dataset `gcT` dari paket R `seqinr` digunakan sebagai validasi eksternal untuk
hubungan genomic GC content dan suhu pertumbuhan optimal. Dataset ini berasal
dari studi klasik Galtier & Lobry (1997) dan berbeda dari TEMPURA maupun fitur
yang diekstrak langsung dari genom lokal.

Penjelasan lebih lengkap untuk metode, hasil, dan cara membahas dataset ini di
paper tersedia di `docs/gct_external_dataset.md`.

Untuk mengimpor dataset:

```powershell
python scripts/import_gct_dataset.py
```

Output utama:

- `data/external/gct/gct_species.csv`: 772 spesies bakteri dengan `GCw` dan `Topt`.
- `data/external/gct/gc16S.csv`, `gctRNA.csv`, `gc23S.csv`, `gc5S.csv`: tabel GC stem RNA dari objek `gcT`.
- `data/external/gct/gct_metadata.md`: sumber dan catatan interpretasi.

Untuk menjalankan validasi eksternal:

```powershell
python scripts/gct_validation.py
```

Validasi ini membuat `figures/gct_external_gc_topt.png` dan mencetak korelasi
Pearson/Spearman antara `GCw` dan `Topt`. Data `gcT` tidak digabung ke
`data/processed/real_gc_data.csv`.

## Output Figure

- `figures/fig1_distribusi_gc.png`: distribusi GC content per kelompok bakteri.
- `figures/fig2_distribusi_topt.png`: distribusi `Topt` per kelompok.
- `figures/fig3_perbedaan_gc.png`: selisih rata-rata GC dan effect size.
- `figures/fig4_scatter.png`: scatter GC content vs `Topt`.
- `figures/fig5_heatmap.png`: heatmap korelasi Pearson.
- `figures/fig6_pgls.png`: perbandingan OLS dan simulasi PGLS.
- `figures/fig8_model_comparison.png`: feature importance RF dan prediksi vs aktual.
- `figures/fig9_roc_cm.png`: ROC curve dan confusion matrix.
- `figures/fig10_single_auc.png`: AUC per fitur tunggal.
- `figures/gct_external_gc_topt.png`: validasi eksternal GCw vs `Topt` pada dataset `gcT`.

Tidak ada `fig7` pada pipeline saat ini.

## Catatan Reproducibility dan Interpretasi

Fitur `GCw`, `GC_tRNA`, `GC_16S`, dan `GC_23S` diekstrak dari genom/anotasi
lokal. Namun `GC_core`, `GC1`, `GC2`, `GC3`, dan hasil PGLS pada
`scripts/gc_analysis.py` saat ini adalah estimasi/simulasi untuk eksplorasi, bukan fitur
hasil ekstraksi langsung dari genom atau analisis filogenetik penuh.

Tabel RNA dari `gcT` berisi stem GC pada level genus dari dataset klasik, jadi
tidak setara langsung dengan fitur `GC_tRNA`, `GC_16S`, dan `GC_23S` yang
diekstrak dari genom lokal di proyek ini.

Karena itu, klaim biologis utama sebaiknya difokuskan pada fitur yang benar-benar
diekstrak dari data genom lokal, terutama `GCw`, `GC_tRNA`, `GC_16S`, dan
`GC_23S`.

## Verifikasi Cepat

```powershell
python -m py_compile scripts/download_genomes.py scripts/extract_features.py scripts/gc_analysis.py scripts/import_gct_dataset.py scripts/gct_validation.py
python scripts/gc_analysis.py
python scripts/import_gct_dataset.py
python scripts/gct_validation.py
```

Setelah selesai, pastikan figure `fig1`, `fig2`, `fig3`, `fig4`, `fig5`,
`fig6`, `fig8`, `fig9`, dan `fig10` tetap terbentuk.

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
- `data/genomes/`: ZIP genom lokal hasil download dari NCBI Datasets.
- `scripts/download_genomes.py`: memilih sampel dan mengunduh genom/anotasi.
- `scripts/extract_features.py`: mengekstrak `GCw`, `GC_tRNA`, `GC_16S`, dan `GC_23S`.
- `scripts/gc_analysis.py`: menjalankan statistik, model machine learning, dan membuat figure.
- `notebooks/gc_analysis.ipynb`: versi notebook dari analisis.
- `figures/`: output visualisasi utama.
- `tools/datasets.exe`: NCBI Datasets CLI lokal.
- `docs/gemini.md`: catatan/ringkasan proyek.

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

Tidak ada `fig7` pada pipeline saat ini.

## Catatan Reproducibility dan Interpretasi

Fitur `GCw`, `GC_tRNA`, `GC_16S`, dan `GC_23S` diekstrak dari genom/anotasi
lokal. Namun `GC_core`, `GC1`, `GC2`, `GC3`, dan hasil PGLS pada
`scripts/gc_analysis.py` saat ini adalah estimasi/simulasi untuk eksplorasi, bukan fitur
hasil ekstraksi langsung dari genom atau analisis filogenetik penuh.

Karena itu, klaim biologis utama sebaiknya difokuskan pada fitur yang benar-benar
diekstrak dari data genom lokal, terutama `GCw`, `GC_tRNA`, `GC_16S`, dan
`GC_23S`.

## Verifikasi Cepat

```powershell
python -m py_compile scripts/download_genomes.py scripts/extract_features.py scripts/gc_analysis.py
python scripts/gc_analysis.py
```

Setelah selesai, pastikan figure `fig1`, `fig2`, `fig3`, `fig4`, `fig5`,
`fig6`, `fig8`, `fig9`, dan `fig10` tetap terbentuk.

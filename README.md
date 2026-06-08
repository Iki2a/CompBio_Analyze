# Analisis GC Content dan Suhu Optimal Bakteri

Proyek ini menganalisis hubungan antara GC content dan suhu pertumbuhan optimal
(`Topt`) pada bakteri termofilik dan mesofilik. Data suhu berasal dari TEMPURA,
sedangkan data genom/anotasi berasal dari NCBI Datasets.

Repo ini sengaja menyimpan CSV, ZIP genom, figure PNG, dan `datasets.exe` agar
analisis bisa dijalankan ulang secara praktis tanpa perlu download ulang semua
data.

## Struktur Proyek

- `200617_TEMPURA.csv`: dataset TEMPURA mentah.
- `selected_genomes_for_download.csv`: 91 genom yang dipilih dari TEMPURA.
- `dataset_genom/`: ZIP genom lokal hasil download dari NCBI Datasets.
- `real_gc_data.csv`: fitur GC hasil ekstraksi dari genom lokal.
- `download_genomes.py`: memilih sampel dan mengunduh genom/anotasi.
- `extract_features.py`: mengekstrak `GCw`, `GC_tRNA`, `GC_16S`, dan `GC_23S`.
- `gc_analysis.py`: menjalankan statistik, model machine learning, dan membuat figure.
- `gc_analysis.ipynb`: versi notebook dari analisis.
- `fig*.png`: output visualisasi utama.

## Setup

Gunakan Python 3.10+.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Jika memakai Jupyter:

```powershell
jupyter notebook gc_analysis.ipynb
```

## Cara Menjalankan Pipeline

Untuk menjalankan analisis dari data yang sudah tersedia:

```powershell
python gc_analysis.py
```

Untuk mengekstrak ulang fitur GC dari ZIP genom lokal:

```powershell
python extract_features.py
```

Untuk memilih sampel dan mengunduh ulang genom memakai NCBI Datasets CLI:

```powershell
python download_genomes.py
```

`download_genomes.py` memakai `datasets.exe` yang sudah ada di root repo.

## Output Figure

- `fig1_distribusi_gc.png`: distribusi GC content per kelompok bakteri.
- `fig2_distribusi_topt.png`: distribusi `Topt` per kelompok.
- `fig3_perbedaan_gc.png`: selisih rata-rata GC dan effect size.
- `fig4_scatter.png`: scatter GC content vs `Topt`.
- `fig5_heatmap.png`: heatmap korelasi Pearson.
- `fig6_pgls.png`: perbandingan OLS dan simulasi PGLS.
- `fig8_model_comparison.png`: feature importance RF dan prediksi vs aktual.
- `fig9_roc_cm.png`: ROC curve dan confusion matrix.
- `fig10_single_auc.png`: AUC per fitur tunggal.

Tidak ada `fig7` pada pipeline saat ini.

## Catatan Reproducibility dan Interpretasi

Fitur `GCw`, `GC_tRNA`, `GC_16S`, dan `GC_23S` diekstrak dari genom/anotasi
lokal. Namun `GC_core`, `GC1`, `GC2`, `GC3`, dan hasil PGLS pada
`gc_analysis.py` saat ini adalah estimasi/simulasi untuk eksplorasi, bukan fitur
hasil ekstraksi langsung dari genom atau analisis filogenetik penuh.

Karena itu, klaim biologis utama sebaiknya difokuskan pada fitur yang benar-benar
diekstrak dari data genom lokal, terutama `GCw`, `GC_tRNA`, `GC_16S`, dan
`GC_23S`.

## Verifikasi Cepat

```powershell
python -m py_compile download_genomes.py extract_features.py gc_analysis.py
python gc_analysis.py
```

Setelah selesai, pastikan figure `fig1`, `fig2`, `fig3`, `fig4`, `fig5`,
`fig6`, `fig8`, `fig9`, dan `fig10` tetap terbentuk.

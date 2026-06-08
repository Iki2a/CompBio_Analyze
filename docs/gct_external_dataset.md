# Dataset Eksternal gcT

## Latar Belakang

Analisis utama proyek ini memakai TEMPURA sebagai sumber suhu pertumbuhan
optimal (`Topt`) dan genom/anotasi NCBI sebagai sumber fitur GC. Untuk
memperkuat interpretasi, proyek ini menambahkan dataset eksternal `gcT` sebagai
pembanding independen. Tujuannya bukan memperbesar data utama, tetapi menguji
apakah pola hubungan antara GC content dan `Topt` tetap terlihat pada dataset
klasik yang dikurasi dari sumber berbeda.

Dataset `gcT` relevan karena berasal dari studi Galtier dan Lobry (1997) tentang
hubungan genomic G+C content, struktur sekunder RNA, dan suhu pertumbuhan
optimal pada prokariota. Dalam proyek ini, `gcT` dipakai sebagai validasi
eksternal untuk pertanyaan sederhana: apakah genomic GC content saja cukup kuat
untuk menjelaskan variasi `Topt` di luar dataset TEMPURA.

## Dataset Apa

`gcT` adalah dataset bawaan paket R `seqinr`. Objek `gcT` berisi beberapa tabel:

- `species`: 772 spesies bakteri dengan genomic GC content dan `Topt`.
- `genus`: ringkasan pada level genus.
- `details`: rentang nilai GC dan `Topt` dari sumber asli.
- `gc16S`, `gctRNA`, `gc23S`, `gc5S`: GC pada stem RNA untuk beberapa kelompok RNA.

Di repo ini, hasil impor disimpan di `data/external/gct/`:

- `gct_species.csv`: tabel utama untuk validasi eksternal.
- `gc16S.csv`, `gctRNA.csv`, `gc23S.csv`, `gc5S.csv`: tabel pendukung RNA stem GC.
- `gct_metadata.md`: ringkasan sumber, citation, dan batas interpretasi.

Kolom utama di `gct_species.csv`:

- `species`: nama spesies.
- `genus`: genus dari spesies tersebut.
- `species_epithet`: epithet spesies.
- `GCw`: genomic GC content dalam persen.
- `Topt`: suhu pertumbuhan optimal.
- `thermal_group`: label sederhana, `Termofil` jika `Topt > 45`, selain itu `Non-termofil`.

Ringkasan data `gct_species.csv` setelah impor:

- Total baris: 772 spesies.
- Non-termofil: 718 spesies.
- Termofil: 54 spesies.
- Rata-rata `GCw`: 51.09%.
- Rata-rata `Topt`: 35.12 degC.
- Rentang `Topt`: 12.00 sampai 105.00 degC.

## Cara Dipakai dalam Paper

Dataset ini paling tepat ditulis sebagai dataset validasi eksternal. Contoh
narasi metode:

> Sebagai validasi eksternal, kami menggunakan dataset `gcT` dari paket R
> `seqinr`, yang berasal dari kompilasi Galtier dan Lobry (1997). Dataset ini
> menyediakan genomic G+C content dan suhu pertumbuhan optimal untuk 772 spesies
> bakteri. Dataset `gcT` tidak digabungkan dengan dataset utama, melainkan
> dianalisis secara terpisah untuk mengevaluasi apakah hubungan GC content dan
> `Topt` konsisten pada sumber data independen.

Contoh narasi hasil:

> Pada dataset eksternal `gcT`, hubungan antara genomic GC content dan `Topt`
> menunjukkan korelasi negatif yang lemah tetapi signifikan (Pearson r = -0.1891,
> p = 1.20e-07; Spearman rho = -0.3406, p = 2.04e-22). Hasil ini mendukung
> interpretasi bahwa genomic GC content secara keseluruhan bukan prediktor tunggal
> yang kuat untuk suhu pertumbuhan optimal.

Contoh narasi diskusi:

> Perbedaan antara hasil dataset utama dan `gcT` menggarisbawahi pentingnya
> membedakan genomic GC content dari GC pada RNA struktural. Pada dataset utama,
> fitur RNA seperti `GC_16S` dan `GC_23S` memperlihatkan hubungan yang lebih kuat
> dengan `Topt`, sedangkan genomic `GCw` sendiri cenderung lemah. Dataset `gcT`
> memberikan pembanding independen yang memperkuat kehati-hatian dalam memakai
> whole-genome GC sebagai biomarker termoadaptasi.

Citation utama:

Galtier, N. & Lobry, J.R. (1997). Relationships between genomic G+C content, RNA
secondary structures, and optimal growth temperature in prokaryotes. Journal of
Molecular Evolution 44:632-636.

Sumber teknis:

`gcT` dataset dari paket R `seqinr` version 4.2-36.

## Perbedaan dengan Dataset Utama

Dataset utama proyek ini:

- Berasal dari TEMPURA untuk data suhu.
- Menggunakan genom/anotasi NCBI lokal untuk menghitung fitur `GCw`, `GC_tRNA`,
  `GC_16S`, dan `GC_23S`.
- Berisi 85 genom hasil ekstraksi fitur setelah filtering dan ketersediaan file.
- Cocok untuk analisis fitur genom aktual dan model prediksi.

Dataset `gcT`:

- Berasal dari studi klasik dan paket `seqinr`.
- Sudah berupa tabel agregat, bukan hasil ekstraksi ulang dari ZIP genom lokal.
- Berisi lebih banyak spesies untuk analisis genomic GC vs `Topt`.
- Tabel RNA di dalamnya adalah stem GC pada level genus, sehingga tidak setara
  langsung dengan fitur RNA hasil ekstraksi genom pada dataset utama.
- Cocok sebagai validasi eksternal, bukan sebagai data training gabungan.

Karena sumber, level agregasi, dan cara pengukuran berbeda, `gcT` tidak dicampur
ke `data/processed/real_gc_data.csv`.

## Hasil Perbandingan

Validasi dijalankan dengan:

```powershell
python scripts/gct_validation.py
```

Hasil pada `gct_species.csv`:

- Pearson r = -0.1891, p = 1.20e-07.
- Spearman rho = -0.3406, p = 2.04e-22.
- Plot tersimpan sebagai `figures/gct_external_gc_topt.png`.

Interpretasinya: pada dataset eksternal `gcT`, genomic GC content memiliki
hubungan negatif yang lemah dengan `Topt`. Ini sejalan dengan kehati-hatian pada
analisis utama bahwa whole-genome GC (`GCw`) saja bukan indikator kuat
termofili. Dalam dataset utama, fitur RNA struktural tetap lebih informatif
dibandingkan `GCw`.

## Hal yang Biasanya Dibahas tentang Dataset

Beberapa poin yang sebaiknya dibahas saat memakai dataset ini:

- **Sumber dan provenance:** jelaskan bahwa `gcT` berasal dari paket `seqinr` dan
  studi Galtier dan Lobry, bukan dari TEMPURA.
- **Unit analisis:** `gct_species.csv` berada pada level spesies, sedangkan tabel
  RNA pendukung berada pada level genus.
- **Definisi fitur:** `GCw` adalah genomic GC content; tabel RNA memakai stem GC,
  bukan hasil ekstraksi fitur dengan pipeline lokal.
- **Ketidakseimbangan kelas:** dataset didominasi non-termofil, dengan 718
  non-termofil dan 54 termofil.
- **Kegunaan utama:** validasi eksternal tren, bukan penggabungan data utama.
- **Batas interpretasi:** korelasi signifikan tidak berarti efek biologis kuat;
  ukuran efek dan sumber variasi tetap harus dibahas.
- **Reproducibility:** dataset diimpor dengan `scripts/import_gct_dataset.py`,
  sehingga proses transformasi dari objek `gcT` ke CSV dapat diulang.


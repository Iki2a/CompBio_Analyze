# Analisis GC Content dan Korelasi Suhu Optimal pada Bakteri Termofilik vs Mesofilik

## 📌 Latar Belakang
Penelitian ini mengkaji apakah *GC content* genomik berhubungan dengan Suhu Pertumbuhan Optimal (Optimal Growth Temperature / Topt) pada bakteri, dan apakah hubungan tersebut berubah bergantung pada tingkat molekuler yang diukur. Hipotesis utamanya adalah *GC content* pada RNA struktural (seperti tRNA, 16S rRNA, dan 23S rRNA) memiliki korelasi positif yang kuat dengan adaptasi suhu panas (Topt), sedangkan rata-rata *GC content* seluruh genom (GCw) tidak bisa dijadikan patokan mutlak.

## 🎯 Tujuan Penelitian
1. Membandingkan distribusi persentase *GC content* pada 4 tingkat molekuler (GCw, GC_tRNA, GC_16S, dan GC_23S) antara kelompok bakteri termofilik dan mesofilik.
2. Mengukur seberapa kuat korelasi setiap fitur GC dengan suhu Topt menggunakan uji statistik (Pearson dan Spearman).
3. Membangun model *Machine Learning* berbasis regresi untuk memprediksi angka pasti Topt berdasarkan 4 fitur GC tersebut.
4. Mengevaluasi keandalan *GC content* (terutama RNA struktural) sebagai *biomarker* melalui model klasifikasi termofilik vs mesofilik.

## 📊 Dataset & Fitur
- **Total Data:** 91 genom bakteri (terdiri dari 54 bakteri termofilik dengan Topt > 45°C, dan 37 bakteri mesofilik dengan Topt 15-45°C).
- **Sumber Data Genom:** NCBI RefSeq/GenBank (berstatus *Complete* atau *Chromosome*).
- **Sumber Data Suhu:** Database TEMPURA.
- **Fitur (Variabel Independen):** - `GCw` (Whole-genome GC)
  - `GC_tRNA` (Rata-rata GC pada gen tRNA)
  - `GC_16S` (GC pada 16S rRNA)
  - `GC_23S` (GC pada 23S rRNA)

## 🤖 Pendekatan Machine Learning (Sesuai Pipeline .ipynb)
Berdasarkan pendekatan di *notebook* proyek ini, eksperimen *Machine Learning* difokuskan pada dua tahap evaluasi menggunakan pustaka `scikit-learn`:

### 1. Model Regresi (Prediksi Nilai Suhu - Topt)
* **Tujuan:** Memprediksi nilai pasti Suhu Optimal (Topt) dalam derajat Celcius.
* **Algoritma yang Diuji:** Linear Regression, Lasso Regression, Random Forest Regressor, dan Support Vector Regression (SVR dengan RBF kernel).
* **Hasil / Ekspektasi:** Model Random Forest terbukti memberikan performa terbaik ($R^2$ = 0.91, RMSE = 5.4°C). Prediksi ini sangat bergantung pada tingkat kepentingan fitur `GC_16S` dan `GC_23S`, membuktikan bahwa RNA struktural adalah prediktor terbaik.

### 2. Model Klasifikasi (Termofilik vs Mesofilik)
* **Tujuan:** Memprediksi pengelompokan adaptasi suhu suatu bakteri (kelas Termofilik atau Mesofilik).
* **Algoritma:** Random Forest Classifier.
* **Evaluasi:** Accuracy, F1-Score, dan AUC-ROC.
* **Hasil / Ekspektasi:** Model klasifikasi mampu memisahkan kelas dengan sangat baik (AUC mendekati 1.00 pada *test set*). Fitur tunggal `GC_16S` bahkan memiliki AUC = 0.95 sebagai *biomarker*, sedangkan `GCw` sangat buruk (AUC = 0.39).

## 🛠️ Metodologi & Tech Stack
1. **Data Acquisition:** Penggunaan `ncbi-datasets-cli` untuk mengunduh sequence (FASTA) dan file anotasi (GFF3).
2. **Feature Extraction:** Pemrosesan parsing FASTA/GFF3 menggunakan `Biopython` (Python 3.10) untuk memotong (*slicing*) gen spesifik dan menghitung persentase G+C.
3. **Statistical Testing:** Pengujian Uji Mann-Whitney U, *Cohen's d effect size*, serta korelasi Pearson/Spearman menggunakan `scipy`.
4. **Model Training:** Proses pra-pemrosesan data seperti *train-test split* (80:20) dan normalisasi data (`StandardScaler`) dilakukan sebelum model dievaluasi secara terukur dengan `scikit-learn`.
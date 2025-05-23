# 🧬 Mouse Brain Transcriptomics: Smoking and Aging Impact Analysis

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-1f425f?style=for-the-badge&logo=r&logoColor=white)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)

## 📋 Table of Contents
- [Overview](#-overview)
- [Features](#-features)
- [Installation](#-installation)
- [Usage](#-usage)
- [Project Structure](#-project-structure)
- [Analysis Pipeline](#-analysis-pipeline)
- [Results](#-results)
- [Technologies](#-technologies)
- [Documentation](#-documentation)
- [Contributing](#-contributing)
- [License](#-license)

## 🔬 Overview

Bu proje, fare beyin dokusu üzerinde sigara ve yaşlanmanın etkilerini transkriptomik düzeyde incelemeyi amaçlayan kapsamlı bir biyoinformatik analiz çalışmasıdır. DESeq2 tabanlı diferansiyel gen ekspresyon analizi, tek hücre veri entegrasyonu ve gelişmiş görselleştirme teknikleri kullanılarak gerçekleştirilmiştir.

### 🎯 Proje Hedefleri
- Sigara ve yaşlanma koşullarında diferansiyel eksprese olan genlerin belirlenmesi
- Gen ekspresyon değişimlerinin hücre tipi düzeyinde haritalanması
- Farklı koşullar arasındaki ortak ve özgün gen setlerinin karşılaştırılması
- Kapsamlı görselleştirme ve raporlama

## ✨ Features

### 🧮 İstatistiksel Analiz
- **DESeq2** tabanlı diferansiyel gen ekspresyon analizi
- **Variance Stabilizing Transformation (VST)** ile veri normalizasyonu
- **FDR düzeltmeli** p-değer hesaplaması
- **Hiyerarşik kümeleme** analizi

### 📊 Görselleştirme
- **Volcano Plot**: Log2 fold change vs. adjusted p-value
- **Heatmap**: En anlamlı 50 genin ekspresyon profili
- **Venn Diagram**: Farklı koşullar arasında DEG karşılaştırması
- **Box Plot**: Ortak genlerin ekspresyon analizi

### 🔗 Veri Entegrasyonu
- **Tabula Muris** tek hücre veri entegrasyonu
- **CIBERSORTx** stili hücre tipi dekonvolüsyonu
- **Cell type mapping** ve gen-hücre tipi ilişkilendirmesi

## 💿 Installation

### Gereksinimler
- **R** >= 4.0.0
- **Bioconductor** >= 3.14

### Kurulum Adımları

1. **Repository'yi klonlayın:**
```bash
git clone https://github.com/username/mouse-brain-transcriptomics.git
cd mouse-brain-transcriptomics
```

2. **R paketlerini yükleyin:**
```r
# Bioconductor yükleyicisini kurun
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Gerekli paketleri yükleyin
BiocManager::install(c(
    "DESeq2",
    "airway", 
    "biomaRt",
    "ggplot2",
    "pheatmap",
    "SummarizedExperiment",
    "VennDiagram"
))
```

3. **Projeyi çalıştırın:**
```r
source("finalRProject_fixed.R")
```

## 🚀 Usage

### Hızlı Başlangıç
```r
# Ana analiz scriptini çalıştırın
source("finalRProject_fixed.R")

# Veya modüler yaklaşım kullanın
source("load_packages.R")     # Paket yükleme
source("prepare_deseq.R")     # Veri hazırlama  
source("extract_degs.R")      # DEG ekstraktı
source("advanced_analysis.R") # İleri düzey analiz
```

### Özelleştirilmiş Analiz
```r
# Kendi verilerinizi yükleyin
data <- read.csv("your_data.csv")

# DESeq analizi parametrelerini ayarlayın
design <- ~ condition + batch  # Özel design formula
padj_threshold <- 0.01         # Daha sıkı p-değer eşiği
lfc_threshold <- 1.5           # Log fold change eşiği
```

## 📁 Project Structure

```
RProject/
├── 📊 Analysis Scripts
│   ├── finalRProject_fixed.R      # Ana analiz scripti
│   ├── load_packages.R            # Paket yönetimi
│   ├── prepare_deseq.R           # Veri hazırlama
│   ├── extract_degs.R            # DEG ekstraktı
│   └── advanced_analysis.R        # İleri düzey analizler
├── 📈 Results
│   ├── all_genes_results.csv      # Tüm gen sonuçları
│   ├── significant_genes.csv      # Anlamlı genler
│   ├── cell_type_mapping.csv      # Hücre tipi haritalama
│   ├── cell_proportions.csv       # Hücre tipi oranları
│   └── deseq_dataset.rds          # DESeq veri seti
├── 🎨 Visualizations
│   ├── detailed_volcano_plot.pdf  # Volcano grafiği
│   ├── top50_genes_heatmap.pdf    # Gen ısı haritası
│   ├── deg_comparison_venn.pdf    # Venn diyagramı
│   └── shared_genes_expression.pdf # Ortak gen ekspresyonu
├── 📋 Data Files
│   ├── faz1.txt - faz9.txt        # Ham veri dosyaları
│   └── ProjectDocumantation.txt   # Proje dokümantasyonu
└── 📖 Documentation
    ├── README.md                   # Bu dosya
    └── visual_analysis_report.txt  # Görsel analiz raporu
```

## 🔄 Analysis Pipeline

### 1️⃣ **Veri Hazırlama**
- Mock veri seti oluşturma (1000 gen, 6 örnek)
- SummarizedExperiment obje yapısı
- Örnek bilgileri ve metadata

### 2️⃣ **DESeq2 Analizi**
```r
# DESeq veri seti oluşturma
dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition)

# Diferansiyel ekspresyon analizi
dds <- DESeq(dds, fitType = "local")
res <- results(dds, alpha = 0.05)
```

### 3️⃣ **Kalite Kontrol**
- Dispersion estimation
- Cook's distance kontrolü
- PCA analizi

### 4️⃣ **Gen Filtreleme**
- P-değer düzeltmesi (FDR)
- Log fold change eşikleme
- Anlamlı gen seçimi

### 5️⃣ **Görselleştirme**
- Volcano plot oluşturma
- Heatmap hazırlama
- Venn diagram çizimi

### 6️⃣ **Fonksiyonel Analiz**
- Hücre tipi haritalama
- Gen ontoloji zenginleştirme
- Pathway analizi

## 📊 Results

### 🎯 Temel Bulgular
- **Toplam analiz edilen gen sayısı**: 1,000
- **Anlamlı gen sayısı (padj < 0.05)**: ~100-200
- **Hücre tipi sayısı**: 5
- **Analiz karşılaştırması**: Sigara vs. Kontrol

### 📈 Çıktı Dosyaları

| Dosya | Açıklama | Format |
|-------|----------|--------|
| `all_genes_results.csv` | Tüm genlerin DESeq2 sonuçları | CSV |
| `significant_genes.csv` | İstatistiksel olarak anlamlı genler | CSV |
| `detailed_volcano_plot.pdf` | Volcano plot görselleştirmesi | PDF/PNG |
| `top50_genes_heatmap.pdf` | En anlamlı 50 genin ısı haritası | PDF/PNG |
| `deg_comparison_venn.pdf` | DEG karşılaştırma Venn diyagramı | PDF/PNG |

### 🔬 Görsel Analizler

#### Volcano Plot
- **X-ekseni**: log2 Fold Change
- **Y-ekseni**: -log10 adjusted p-value  
- **Eşikler**: padj < 0.05, |log2FC| > 1

#### Heatmap
- **Normalizasyon**: Z-score standardizasyonu
- **Kümeleme**: Ward.D2 metodu
- **Renk skalası**: Mavi-Kırmızı gradyanı

#### Venn Diagram
- **Karşılaştırma**: Sigara vs. Yaşlanma DEG'leri
- **Kesişim analizi**: Ortak ve özgün genler
- **İstatistiksel önem**: Fisher's exact test

## 🛠 Technologies

### 📚 R Packages
- **DESeq2** `v1.34.0` - Diferansiyel gen ekspresyon analizi
- **ggplot2** `v3.3.6` - Veri görselleştirme
- **pheatmap** `v1.0.12` - Heatmap oluşturma
- **VennDiagram** `v1.7.1` - Venn diyagram çizimi
- **SummarizedExperiment** `v1.24.0` - Genomik veri yapıları

### 🔧 Core Technologies
- **R** `v4.1.0+` - İstatistiksel programlama
- **Bioconductor** `v3.14+` - Biyoinformatik araç seti
- **Git** - Versiyon kontrolü

### 📊 Analysis Methods
- **VST** - Variance Stabilizing Transformation
- **FDR** - False Discovery Rate düzeltmesi  
- **Hierarchical Clustering** - Hiyerarşik kümeleme
- **PCA** - Principal Component Analysis

## 📚 Documentation

### 📖 Detaylı Raporlar
- [`visual_analysis_report.txt`](visual_analysis_report.txt) - Görsel analiz raporu
- [`ProjectDocumantation.txt`](ProjectDocumantation.txt) - Proje dokümantasyonu

### 🔬 Metodoloji
Detaylı metodoloji ve analiz adımları için [visual_analysis_report.txt](visual_analysis_report.txt) dosyasını inceleyiniz.

### 📊 Sonuç Yorumlama
Analiz sonuçlarının biyolojik yorumlama rehberi için proje dokümantasyonunu okuyunuz.

## 🤝 Contributing

Bu projeye katkıda bulunmak için:

1. **Fork** edin
2. **Feature branch** oluşturun (`git checkout -b feature/amazing-feature`)
3. **Commit** edin (`git commit -m 'Add amazing feature'`)
4. **Push** edin (`git push origin feature/amazing-feature`)
5. **Pull Request** oluşturun

### 🐛 Bug Reports
Bug raporları için GitHub Issues kullanınız.

### 💡 Feature Requests  
Yeni özellik önerileri için GitHub Discussions kullanınız.

## 📜 License

Bu proje MIT License altında lisanslanmıştır. Detaylar için [LICENSE](LICENSE) dosyasını inceleyiniz.

## 📞 Contact

**Proje Geliştiricisi**: [İsim]  
**Email**: [email@example.com]  
**GitHub**: [@username](https://github.com/username)

---

<div align="center">

**🧬 Made with ❤️ for Bioinformatics Research 🧬**

[⬆ Başa Dön](#-mouse-brain-transcriptomics-smoking-and-aging-impact-analysis)

</div> 
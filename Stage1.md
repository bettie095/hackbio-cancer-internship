<!--StartFragment-->


### **Development of a RNA-Seq Based Prognostic Signature**

### **in Lung Adenocarcinoma**

Authors (@slack): Manal Agdada (@Manal), Elizabeth Abodunrin (@Bettie\_95), Hagar Haitham Elazab (@HBONH33), Rahma Nabil Sallam (@rahmanabil2002)

Github Repo: https://github.com/manal-agdada/Stage-1---HackBio-internship/blob/main/stage1.md

LinkedIn post: https://www.linkedin.com/posts/rahma-nabil-27a080225_hackbio-research-rnaseq-activity-7237856063700537345-D6ol?utm_source=share&utm_medium=member_android

Lung adenocarcinoma is a major subtype of non-small lung cancer (NSCLC) with a high mortality rate \[1]. Previously, microarrays have been widely used for gene expression analysis and the development of prognostic signatures in lung adenocarcinoma. However, microarray-derived signatures are limited by biases, as they do not account for non-coding RNA and low abundance transcripts \[2]. To address these limitations, Sudhanshu et al. used RNA-seq data from a lung adenocarcinoma cohort to develop a robust prognostic gene signature. This signature could be directly incorporated into RNA-seq-based clinical tests to improve the prediction of patient outcomes and to guide personalized treatment decisions \[3].

**Development of a prognostic model for lung adenocarcinoma**

The authors analyzed RNA-seq data from The Cancer Genome Atlas (TCGA) patients with lung adenocarcinoma to develop a prognostic signature that classifies lung adenocarcinoma patients into high- and low-risk groups. In particular, their analysis delved into utilizing Cox Proportional Hazard Regression Models to identify four genes (RHOV, CD109, LINC00941 and FRRS1) that make the prognostic signature that helps to predict overall and metastasis-free survival of lung adenocarcinoma patients. 

The built model was used to calculate risk scores for each patient included in the dataset to stratify patients into groups based on the risk score and expression of these four genes. The derived signature was validated in independent cohorts, showing its ability to stratify patients into high- and low-risk groups.

Survival analysis of the overall and metastasis-free survival in independent cohorts showed a higher risk of mortality in high-risk patients than low-risk patients (Figure 1).

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdcoYINxNNYqEuMR7k3bFTlwisCjq5r4vONiT6kr_ZdMTdAuSWYKd1K4U64NaynAkHN4OEKfpyvAR9AHHMqFUEHMk1wFv4eQAPgl0tLofEvs0t508LqtS1xszK1EASa6LqLyQr5wB4Zpgxj0SoQJXX3_Lm8?key=YLu_DE8OGL6Caiap1UzJTw)

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeqxTo6nI_fS__4Th_rGQPSVkletiiA8LDbWY8ofSEJOfEdswxPYbon9jdcwjrOrr0ceSBIOzkjTDP9vfX7dWzIN_jeyz7wmuWsE5ACzGJQzOtSKb1UYSVlQOe7ka_TJMcKJqLbEtCnPkJGntxp-dKy15D1?key=YLu_DE8OGL6Caiap1UzJTw)


### **Possible clinical applications of the four-gene signature**

The four-gene signature could be integrated into clinical sequencing programs that are already applied in clinical practice to guide treatment decisions, especially for early-stage lung adenocarcinoma patients that have worse outcomes than other patients. This could help clinicians in determining who might benefit from standard chemotherapy or more aggressive treatments based on their risk category.


### **Conclusion**

The four-gene prognostic signature identified in lung adenocarcinoma has been shown to significantly predict survival outcomes and could be integrated into RNA-seq-based clinical tests. Given that early stage lung adenocarcinoma often has poorer survival rates than other cancers, more aggressive and efficient treatment approaches are needed. This signature holds the potential to personalize therapy and improve patient outcomes in both early and advanced stages of the disease.

### **References** 

1\. Pao W, Girard N. New driver mutations in non-small-cell lung cancer. Lancet Oncol. 2011;12(2):175–180.

2\. Zheng Y, Bueno R. Commercially available prognostic molecular models in early-stage lung cancer: A review of the Pervenio Lung RS and Myriad myPlan Lung Cancer tests. Expert Rev Mol Diagn. 2015;15(5):589–596.

3\. Sudhanshu Shukla, Joseph R Evans, Rohit Malik, Felix Y Feng, Saravana M Dhanasekaran, Xuhong Cao, Guoan Chen, David G Beer, Hui Jiang, Arul M Chinnaiyan, Development of a RNA-Seq Based Prognostic Signature in Lung Adenocarcinoma, JNCI: Journal of the National Cancer Institute, 109(1) 2017





<!--EndFragment-->

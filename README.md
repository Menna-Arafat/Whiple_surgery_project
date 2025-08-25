# Proteomic Profiling of Pancreatic Cancer Patients: Time Series Study Before and After Whipple Surgery
In this study, I applied a network-based analysis to characterize temporal protein dynamics and functional organization.
Implemented Weighted Gene Co-expression Network Analysis (WGCNA) to identify co-expressed protein modules, correlating module eigengenes with sample metadata to detect modules associated with temporal progression.
Modules showing significant changes, including the blue module at 7 days post-operation, were analyzed using generalized additive models, and hub proteins with high intramodular connectivity (datKME > 0.85) were identified as central drivers.
Functional interpretation was achieved by integrating the blue module with protein–protein interaction networks from STRING, followed by enrichment analysis using clusterProfiler and modularity clustering with the Louvain algorithm, revealing pathway modules with crosstalk and coordinated biological activity.

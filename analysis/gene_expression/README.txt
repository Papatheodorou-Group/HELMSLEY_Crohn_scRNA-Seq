step1: RUN_clustering_ccRegress.sh
computes clustering solutions for different HVG selection and clustering criteria 

step2: RUN_DEA_ccRegress.sh
run DEA for the selected solutions (best silhouette score and min number of clusters)

step3: RUN_filtering.sh
remove low quality clusters, find doublets and merge objects

step4: RUN_Crohn_filtering.sh
remove low quality clusters, find doublets and merge objects

step5: RUN_filtering_2ndRound.sh
remove doublets - all of them or by cluster

step6: RUN_reclustering.sh
retrieve external cell annotation (from scANVI), recluster cells and evaluate clustering solutions

step7: RUN_DEA_reclustering.sh
run DEA for the subclusters

step8: RUN_reannotation.sh
annotation curation based on GCA information, clustering, and cluster markers

step9: RUN_subclustering.sh
re-cluster selected clusters into subclusters

step10: RUN_annotation_final.sh



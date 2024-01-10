from sklearn.cluster import DBSCAN
import numpy as np

def extract_cluster_labels_dbscan(vae_mu, data, nameCluster, eps, min_samples=30):
    dbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(vae_mu)
    cluster_labels = dbscan.labels_
    
    # Number of clusters in labels, ignoring noise if present.
    n_clusters = len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)
    n_noise = list(cluster_labels).count(-1)
    
    print('Estimated number of clusters: %d' % n_clusters)
    print('Estimated number of noise points: %d' % n_noise)

    for cluster_id in range(n_clusters):
        indices = np.where(cluster_labels == cluster_id)[0]
        sequence_names = data["name"][indices]
        
        if len(indices) > 30:  # Write clusters with more than 30 reads (30x)
            # Add 1 to cluster_id for filename
            with open(f'{nameCluster}_{cluster_id + 1}-.txt', 'w') as f:
                for name in sequence_names:
                    f.write(name + '\n')

    # Saving noise points
    noise_indices = np.where(cluster_labels == -1)[0]
    noise_sequence_names = data["name"][noise_indices]
    with open(nameCluster + 'noise.txt', 'w') as f:
        for name in noise_sequence_names:
            f.write(name + '\n')

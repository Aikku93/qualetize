/************************************************/
#include <math.h>
/************************************************/

//! Insert element to start of linked list
static void InsertAtListHead(uint32_t Idx, uint32_t *Head, uint32_t *List) {
	List[Idx] = *Head;
	*Head = Idx;
}

//! Add/remove a cluster in a distortion linked list (Head = Most distorted)
static void InsertToDistortedClusterList(TCluster_t *Clusters, uint32_t Idx, uint32_t *Head) {
	uint32_t Prev = CLUSTER_END_OF_LIST;
	uint32_t Next = *Head;
	float Dist = Clusters[Idx].TotalDist;
	if(Dist != 0.0f) {
		while(Next != CLUSTER_END_OF_LIST && Dist < Clusters[Next].TotalDist) {
			Prev = Next;
			Next = Clusters[Next].NextCluster;
		}
		Clusters[Idx].NextCluster = Next;
		if(Prev != CLUSTER_END_OF_LIST) Clusters[Prev].NextCluster = Idx;
		else *Head = Idx;
	}
}
static uint32_t PopDistortedClusterList(TCluster_t *Clusters, uint32_t *Head) {
	uint32_t x = *Head;
	*Head = Clusters[x].NextCluster;
	return x;
}

//! Add/remove a cluster in the "empty clusters" linked list
static void InsertToEmptyClusterList(TCluster_t *Clusters, uint32_t Idx, uint32_t *Head) {
	Clusters[Idx].NextCluster = *Head;
	*Head = Idx;
}
static uint32_t PopEmptyClusterList(TCluster_t *Clusters, uint32_t *Head) {
	uint32_t x = *Head;
	*Head = Clusters[x].NextCluster;
	return x;
}

/************************************************/

//! Clear training data
static inline void ClearTraining(TCluster_t *x, uint32_t nDims) {
	x->nPoints = 0;
	x->TrainWeight = 0.0f;
	TCluster_ClearTrainingVector(x, nDims);
}

//! Train cluster on data
static inline void TrainCluster(TCluster_t *x, const TClusterData_t *Data, float w, uint32_t nDims) {
	TCluster_AddToTraining(x, Data, w, nDims);
	x->TrainWeight += w;
	x->nPoints++;
}

//! Calculate cluster centroid and recalculate distortions
static uint8_t ResolveCluster(
	TCluster_t *Cluster,
	const TClusterData_t *Data,
	const uint32_t *ClusterListIndices,
	const float *Weights,
	uint32_t nDims
) {
	//! Resolve for the centroid
	if(Cluster->nPoints == 0 || Cluster->TrainWeight == 0.0f) return 0;
	TCluster_ResolveCentroid(Cluster, nDims);

	//! Recalculate distortion
	uint32_t Next = Cluster->FirstDataIdx;
	uint32_t PeakDistIdx = CLUSTER_END_OF_LIST;
	float    PeakDistVal = 0.0f;
	float    TotalDist   = 0.0f;
	while(Next != CLUSTER_END_OF_LIST) {
		float w = Weights ? Weights[Next] : 1.0f;
		float Dist = w * TCluster_Dist2ToCentroid(Cluster, Data + Next*nDims, nDims);
		if(Dist > PeakDistVal) {
			PeakDistIdx = Next;
			PeakDistVal = Dist;
		}
		TotalDist += Dist;
		Next = ClusterListIndices[Next];
	}
	Cluster->MaxDistIdx = PeakDistIdx;
	Cluster->TotalDist  = TotalDist;
	return 1;
}

/************************************************/

//! Apply cluster analysis to the specified data
static inline uint32_t TClusterize_Process(
	TCluster_t *Clusters,
	const TClusterData_t *Data,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	uint32_t nPasses,
	const float *Weights,
	uint32_t nDims
) {
	uint32_t n, k;
	if(!nClusters || !nDataPoints || !nPasses || !nDims) return 0;

	//! Do a quick pass for the first cluster
	uint32_t nCurrentClusters = 1;
	uint32_t DistClusterHead  = CLUSTER_END_OF_LIST;
	uint32_t EmptyClusterHead = CLUSTER_END_OF_LIST;
	Clusters[0].FirstDataIdx  = CLUSTER_END_OF_LIST;
	ClearTraining(&Clusters[0], nDims);
	for(n=0;n<nDataPoints;n++) {
		float w = Weights ? Weights[n] : 1.0f;
		InsertAtListHead(n, &Clusters[0].FirstDataIdx, ClusterListIndices);
		TrainCluster(&Clusters[0], &Data[n], w, nDims);
	}
	ResolveCluster(&Clusters[0], Data, ClusterListIndices, Weights, nDims);
	InsertToDistortedClusterList(Clusters, 0, &DistClusterHead);

	//! Begin creating additional clusters
	while(nCurrentClusters < nClusters && DistClusterHead != CLUSTER_END_OF_LIST) {
		//! Create new clusters from the most distorted data points
		//! Note that we are splitting out the most distorted point
		//! of the most distorted cluster, NOT the most distorted
		//! point general. This is an important distinction.
		float MaxDist = Clusters[DistClusterHead].TotalDist;
		do {
			//! Stop splitting once we've reached clusters with low distortion.
			//! The idea is to only split highly distorted clusters, so that
			//! we avoid early splitting of low-distortion clusters, as this
			//! would give sub-optimal results.
			uint32_t SrcCluster = PopDistortedClusterList(Clusters, &DistClusterHead);
			if(Clusters[SrcCluster].TotalDist < 0.5f*MaxDist) break;

			uint32_t DstCluster = nCurrentClusters++;
			TCluster_SetCentroidToData(&Clusters[DstCluster], Data + Clusters[SrcCluster].MaxDistIdx*nDims, nDims);
		} while(nCurrentClusters < nClusters && DistClusterHead != CLUSTER_END_OF_LIST);

		//! Begin refinement loop
		uint32_t Pass;
		float LastPassDist = INFINITY;
		for(Pass=0;Pass<nPasses;Pass++) {
			//! Assign data points to clusters
			for(k=0;k<nCurrentClusters;k++) {
				Clusters[k].FirstDataIdx = CLUSTER_END_OF_LIST;
				ClearTraining(&Clusters[k], nDims);
			}
			for(n=0;n<nDataPoints;n++) {
				uint32_t BestIdx  = 0;
				float    BestDist = INFINITY;
				for(k=0;k<nCurrentClusters;k++) {
					float Dist = TCluster_Dist2ToCentroid(&Clusters[k], Data + n*nDims, nDims);
					if(Dist < BestDist) {
						BestIdx  = k;
						BestDist = Dist;
					}
				}
				float w = Weights ? Weights[n] : 1.0f;
				InsertAtListHead(n, &Clusters[BestIdx].FirstDataIdx, ClusterListIndices);
				TrainCluster(&Clusters[BestIdx], Data + n*nDims, w, nDims);
			}

			//! Resolve all clusters
			float ThisPassDist = 0.0f;
			DistClusterHead  = CLUSTER_END_OF_LIST;
			EmptyClusterHead = CLUSTER_END_OF_LIST;
			for(k=0;k<nCurrentClusters;k++) {
				if(ResolveCluster(&Clusters[k], Data, ClusterListIndices, Weights, nDims)) {
					//! If the cluster resolves, update the distortion linked list
					InsertToDistortedClusterList(Clusters, k, &DistClusterHead);
					ThisPassDist += Clusters[k].TotalDist;
				} else {
					//! No resolve - append to empty-cluster linked list
					InsertToEmptyClusterList(Clusters, k, &EmptyClusterHead);
				}
			}

			//! Split the most distorted data points into empty clusters
			while(DistClusterHead != CLUSTER_END_OF_LIST && EmptyClusterHead != CLUSTER_END_OF_LIST) {
				uint32_t SrcCluster = PopDistortedClusterList(Clusters, &DistClusterHead);
				uint32_t DstCluster = PopEmptyClusterList    (Clusters, &EmptyClusterHead);
				TCluster_SetCentroidToData(&Clusters[DstCluster], Data + Clusters[SrcCluster].MaxDistIdx*nDims, nDims);
			}

			//! Early exit on convergence
			if(ThisPassDist == 0.0f || ThisPassDist >= LastPassDist) break;

			//! Set up for next iteration
			LastPassDist = ThisPassDist;
		}
	}
	return nCurrentClusters;
}

/************************************************/

//! Get cluster indices from clusters and data lists
static inline void TClusterize_GetClusterIndices_u8(
	uint8_t *DataClusterIndices,
	const TCluster_t *Clusters,
	uint32_t nClusters,
	const uint32_t *ClusterListIndices
) {
	uint32_t k;
	for(k=0;k<nClusters;k++) {
		uint32_t Next = Clusters[k].FirstDataIdx;
		while(Next != CLUSTER_END_OF_LIST) {
			DataClusterIndices[Next] = (uint8_t)k;
			Next = ClusterListIndices[Next];
		}
	}
}

/************************************************/
//! EOF
/************************************************/

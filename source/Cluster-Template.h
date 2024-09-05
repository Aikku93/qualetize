/************************************************/
#include <math.h>
/************************************************/

//! Insert element to start of linked list
static void InsertAtListHead(uint32_t Idx, uint32_t *Head, uint32_t *List) {
	List[Idx] = *Head;
	*Head = Idx;
}

//! Remove element from linked list, given Prev and Next
//! ie. this joins Prev and Next together in the list
static void RemoveFromList(uint32_t Prev, uint32_t Next, uint32_t *Head, uint32_t *List) {
	if(Prev != CLUSTER_END_OF_LIST) List[Prev] = Next;
	else *Head = Next;
}

//! Add/remove a cluster in a distortion linked list (Head = Most distorted)
static void InsertToDistortedClusterList(TCluster_t *Clusters, uint32_t Idx, uint32_t *Head) {
	uint32_t Prev = CLUSTER_END_OF_LIST;
	uint32_t Next = *Head;
	float Dist = Clusters[Idx].MaxDistVal;
	if(Dist != 0.0f) {
		while(Next != CLUSTER_END_OF_LIST && Dist < Clusters[Next].MaxDistVal) {
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

//! ADd/remove a cluster in the "empty clusters" linked list
//! NOTE: We don't use PrevCluster for this.
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
	TCluster_ClearTrainingVector(x, nDims);
}

//! Train cluster on data
static inline void TrainCluster(TCluster_t *x, const TClusterData_t *Data, uint32_t nDims) {
	TCluster_AddToTraining(x, Data, nDims);
	x->nPoints++;
}

//! Calculate cluster centroid and recalculate distortions
static uint8_t ResolveCluster(
	TCluster_t *Cluster,
	const TClusterData_t *Data,
	const uint32_t *ClusterListIndices,
	uint32_t nDims
) {
	//! Resolve for the centroid
	uint32_t nPoints = Cluster->nPoints;
	if(nPoints == 0) return 0;
	TCluster_ResolveCentroid(Cluster, nDims);

	//! Recalculate distortion
	//! Note that MaxDistVal stores the average between the average
	//! distortion, and the peak distortion. This generally reduces
	//! the chances of overfitting, while still giving a clearer idea
	//! of which cluster needs splitting most urgently.
	uint32_t Prev = CLUSTER_END_OF_LIST, Next = Cluster->FirstDataIdx;
	uint32_t PeakDistIdx = CLUSTER_END_OF_LIST, PeakDistIdxPrev = CLUSTER_END_OF_LIST;
	float    PeakDistVal = 0.0f;
	float    TotalDist   = 0.0f;
	while(Next != CLUSTER_END_OF_LIST) {
		float Dist = TCluster_Dist2ToCentroid(Cluster, Data + Next*nDims, nDims);
		if(Dist > PeakDistVal) {
			PeakDistIdx = Next, PeakDistIdxPrev = Prev;
			PeakDistVal = Dist;
		}
		TotalDist += Dist;
		Prev = Next;
		Next = ClusterListIndices[Next];
	}
	Cluster->MaxDistIdx = PeakDistIdx, Cluster->MaxDistIdxPrev = PeakDistIdxPrev;
	Cluster->MaxDistVal = PeakDistVal + (TotalDist / (float)nPoints);
	Cluster->TotalDist  = TotalDist;
	return 1;
}

//! Split cluster into two separate clusters
static void SplitCluster(
	uint32_t SrcClusterIdx,
	uint32_t DstClusterIdx,
	TCluster_t *Clusters,
	const TClusterData_t *Data,
	uint32_t *ClusterListIndices,
	uint32_t *DistClusterHead,
	uint32_t nDims
) {
	TCluster_t *SrcCluster = &Clusters[SrcClusterIdx];
	TCluster_t *DstCluster = &Clusters[DstClusterIdx];

	//! Create a new cluster from the most distorted point - this helps
	//! us make it out of a local optimum into a better cluster fit.
	TCluster_SetCentroidToData(DstCluster, Data + SrcCluster->MaxDistIdx*nDims, nDims);
	DstCluster->FirstDataIdx = CLUSTER_END_OF_LIST;
	ClearTraining(DstCluster, nDims);

	//! Now re-cluster the data points that fell into the source cluster
	uint32_t Prev = CLUSTER_END_OF_LIST;
	uint32_t Next = SrcCluster->FirstDataIdx;
	ClearTraining(SrcCluster, nDims);
	while(Next != CLUSTER_END_OF_LIST) {
		uint32_t Idx = Next;
		Next = ClusterListIndices[Idx];

		//! Decide whether to use the original cluster or the newly-created one
		const TClusterData_t *ThisData = Data + Idx*nDims;
		float DistSrc = TCluster_Dist2ToCentroid(SrcCluster, ThisData, nDims);
		float DistDst = TCluster_Dist2ToCentroid(DstCluster, ThisData, nDims);
		if(DistSrc < DistDst) {
			Prev = Idx;
			TrainCluster(SrcCluster, ThisData, nDims);
		} else {
			//! Remove from SrcCluster list before inserting to DstCluster list
			RemoveFromList(Prev, Next, &SrcCluster->FirstDataIdx, ClusterListIndices);
			InsertAtListHead(Idx, &DstCluster->FirstDataIdx, ClusterListIndices);
			TrainCluster(DstCluster, ThisData, nDims);
		}
	}

	//! Resolve the new cluster centroids
	//! NOTE: Do not re-insert these clusters into the distortion
	//! list, or we end up "prematurely optimizing" the result,
	//! causing technically lower RMSE, at the cost of nasty
	//! overfitting artifacts (where the algorithm will assign
	//! very few points to a cluster, even if that increases the
	//! error for other clusters significantly).
	(void)DistClusterHead;
	if(ResolveCluster(SrcCluster, Data, ClusterListIndices, nDims)) {
		//InsertToDistortedClusterList(Clusters, SrcClusterIdx, DistClusterHead);
	}
	if(ResolveCluster(DstCluster, Data, ClusterListIndices, nDims)) {
		//InsertToDistortedClusterList(Clusters, DstClusterIdx, DistClusterHead);
	}
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
	float AvgErrorThreshold,
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
		InsertAtListHead(n, &Clusters[0].FirstDataIdx, ClusterListIndices);
		TrainCluster(&Clusters[0], &Data[n], nDims);
	}
	ResolveCluster(&Clusters[0], Data, ClusterListIndices, nDims);
	InsertToDistortedClusterList(Clusters, 0, &DistClusterHead);

	//! Now begin splitting the most distorted clusters until we
	//! have the desired number of output clusters. Note that
	//! we only enter this loop if we had any distortion when
	//! we calculated the first cluster (otherwise we've already
	//! achieved perfect results and don't need to split).
	float LastSplitDist = Clusters[0].TotalDist;
	if(LastSplitDist != 0.0f) while(nCurrentClusters < nClusters) {
		//! Split the most distorted clusters out
		//! We can split as many times as we want here, with
		//! N=1..nClusters-1. In general, N=1 gives the best
		//! results, but can take absurd amounts of time when
		//! nClusters is large, whereas N=nClusters-1 still
		//! gives good results, just somewhat sub-optimal as
		//! compared to N=1 (-1dB PSNR ish), but is faster.
		uint32_t nNewClustersRem = nCurrentClusters;
		if(nCurrentClusters+nNewClustersRem > nClusters) nNewClustersRem = nClusters - nCurrentClusters;
		while(nNewClustersRem) {
			//! Find a distorted cluster to split from
			uint32_t SrcCluster;
			if(DistClusterHead != CLUSTER_END_OF_LIST) {
				SrcCluster = PopDistortedClusterList(Clusters, &DistClusterHead);
			} else break;

			//! If we have no empty clusters, then generate a new one
			uint32_t DstCluster;
			if(EmptyClusterHead != CLUSTER_END_OF_LIST) {
				DstCluster = PopEmptyClusterList(Clusters, &EmptyClusterHead);
			} else {
				DstCluster = nCurrentClusters;
				nCurrentClusters++;
				nNewClustersRem--;
			}

			//! Actually split the cluster
			SplitCluster(
				SrcCluster,
				DstCluster,
				Clusters,
				Data,
				ClusterListIndices,
				&DistClusterHead,
				nDims
			);
		}

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
				InsertAtListHead(n, &Clusters[BestIdx].FirstDataIdx, ClusterListIndices);
				TrainCluster(&Clusters[BestIdx], Data + n*nDims, nDims);
			}

			//! Resolve all clusters
			DistClusterHead  = CLUSTER_END_OF_LIST, DistClusterHead = CLUSTER_END_OF_LIST;
			EmptyClusterHead = CLUSTER_END_OF_LIST;
			for(k=0;k<nCurrentClusters;k++) {
				if(ResolveCluster(&Clusters[k], Data, ClusterListIndices, nDims)) {
					//! If the cluster resolves, update the distortion linked list
					InsertToDistortedClusterList(Clusters, k, &DistClusterHead);
				} else {
					//! No resolve - append to empty-cluster linked list
					InsertToEmptyClusterList(Clusters, k, &EmptyClusterHead);
				}
			}

			//! Split the most distorted clusters into any empty ones
			while(DistClusterHead != CLUSTER_END_OF_LIST && EmptyClusterHead != CLUSTER_END_OF_LIST) {
				uint32_t SrcCluster = PopDistortedClusterList(Clusters, &DistClusterHead);
				uint32_t DstCluster = PopEmptyClusterList    (Clusters, &EmptyClusterHead);
				SplitCluster(
					SrcCluster,
					DstCluster,
					Clusters,
					Data,
					ClusterListIndices,
					&DistClusterHead,
					nDims
				);
			}

			//! Calculate global distortion and early exit
			float ThisPassDist = 0.0f;
			for(k=0;k<nCurrentClusters;k++) ThisPassDist += Clusters[k].TotalDist;
			ThisPassDist /= (float)nDataPoints;
			if(ThisPassDist == 0.0f || ThisPassDist >= LastPassDist) {
				LastPassDist = ThisPassDist;
				break;
			}
			LastPassDist = ThisPassDist;
		}

		//! Early exit?
		if(
			LastPassDist <= (AvgErrorThreshold * (float)nDataPoints) ||
			LastPassDist >= LastSplitDist
		) break;
		LastSplitDist = LastPassDist;
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

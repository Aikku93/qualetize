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
	x->TrainWeight = 0.0f;
	TCluster_ClearTrainingVector(x, nDims);
}

//! Train cluster on data
static inline void TrainCluster(TCluster_t *x, const TClusterData_t *Data, float w, uint32_t nDims) {
	TCluster_AddToTraining(x, Data, w, nDims);
	x->TrainWeight += w;
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
	if(Cluster->TrainWeight == 0.0f) return 0;
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

//! Split cluster in two
static void SplitCluster(
	TCluster_t *DstCluster,
	TCluster_t *SrcCluster,
	const TClusterData_t *Data,
	const uint32_t *ClusterListIndices,
	const float *Weights,
	uint32_t nDims
) {
	//! First, set the destination cluster to the most distorted point
	TCluster_SetCentroidToData(DstCluster, Data + SrcCluster->MaxDistIdx*nDims, nDims);

	//! Now find the point furthest away from that for the source cluster
	uint32_t Next = SrcCluster->FirstDataIdx;
	uint32_t MaxDistIdx = Next;
	float    MaxDistVal = -INFINITY;
	while(Next != CLUSTER_END_OF_LIST) {
		float w = Weights ? Weights[Next] : 1.0f;
		float d = w * TCluster_Dist2ToCentroid(DstCluster, Data + Next*nDims, nDims);
		if(d > MaxDistVal) {
			MaxDistIdx = Next;
			MaxDistVal = d;
		}
		Next = ClusterListIndices[Next];
	}
	TCluster_SetCentroidToData(SrcCluster, Data + MaxDistIdx*nDims, nDims);
}

/************************************************/

//! Create initial clusters by splitting out the most
//! distorted point into a new cluster and reclustering,
//! repeating until we have as many clusters as needed.
//! The output is usable, but not particularly great.
//! Note: 1/Variance is stored to InvVariancePtr only
//! if the distortion is not 0.0!
static uint32_t CreateClusters(
	TCluster_t *Clusters,
	const TClusterData_t *Data,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	uint32_t nFinalPasses,
	const float *Weights,
	float *InvVariancePtr,
	uint32_t nDims
) {
	uint32_t n, k;

	//! Do a quick pass for the first cluster
	uint32_t nCurrentClusters = 1;
	uint32_t DistClusterHead  = CLUSTER_END_OF_LIST;
	uint32_t EmptyClusterHead = CLUSTER_END_OF_LIST;
	Clusters[0].FirstDataIdx  = CLUSTER_END_OF_LIST;
	ClearTraining(&Clusters[0], nDims);
	for(n=0;n<nDataPoints;n++) {
		float w = Weights ? Weights[n] : 1.0f;
		InsertAtListHead(n, &Clusters[0].FirstDataIdx, ClusterListIndices);
		TrainCluster(&Clusters[0], Data + n*nDims, w, nDims);
	}
	ResolveCluster(&Clusters[0], Data, ClusterListIndices, Weights, nDims);
	InsertToDistortedClusterList(Clusters, 0, &DistClusterHead);

	//! Because we averaged all the data together at a single point,
	//! the variance of the dataset is equal to TotalDist/TrainWeight.
	if(Clusters[0].TotalDist != 0.0f) {
		*InvVariancePtr = Clusters[0].TrainWeight / Clusters[0].TotalDist;
	}

	//! Begin creating additional clusters
	float LastPassDist = Clusters[0].TotalDist;
	while(nCurrentClusters < nClusters && DistClusterHead != CLUSTER_END_OF_LIST) {
		//! Create new cluster from the most distorted data point.
		//! Note that we are splitting out the most distorted point
		//! of the most distorted cluster, NOT the most distorted
		//! point general. This is an important distinction.
		//! Also note the commented-out code: This allows splitting
		//! multiple clusters at once, but because of the very low
		//! number of iterations we do, the results are quite bad.
		do {
			uint32_t SrcCluster = PopDistortedClusterList(Clusters, &DistClusterHead);
			uint32_t DstCluster = nCurrentClusters++;
			SplitCluster(&Clusters[DstCluster], &Clusters[SrcCluster], Data, ClusterListIndices, Weights, nDims);
		} while(/*nCurrentClusters < nClusters && DistClusterHead != CLUSTER_END_OF_LIST*/0);

		//! Run at least two k-means passes to improve results a bit
		uint32_t Pass, nPasses = 2;
		if(nCurrentClusters == nClusters && nFinalPasses > nPasses) {
			nPasses = nFinalPasses;
		}
		for(Pass=0;Pass<nPasses;Pass++) {
			//! Assign points to clusters
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

			//! Resolve clusters
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

			//! Terminate if we exceeded the maximum number of passes.
			//! This happens if the last iteration had empty clusters
			//! so that we can assign /something/ to them.
			if(nPasses >= nFinalPasses) break;

			//! Split the most distorted data points into empty clusters
			if(DistClusterHead != CLUSTER_END_OF_LIST && EmptyClusterHead != CLUSTER_END_OF_LIST) {
				//! Do at least one more pass to assign to this cluster
				nPasses++;
			}
			while(DistClusterHead != CLUSTER_END_OF_LIST && EmptyClusterHead != CLUSTER_END_OF_LIST) {
				uint32_t SrcCluster = PopDistortedClusterList(Clusters, &DistClusterHead);
				uint32_t DstCluster = PopEmptyClusterList    (Clusters, &EmptyClusterHead);
				SplitCluster(&Clusters[DstCluster], &Clusters[SrcCluster], Data, ClusterListIndices, Weights, nDims);
			}

			//! Early exit on convergence
			if(ThisPassDist == 0.0f || ThisPassDist == LastPassDist) break;
			LastPassDist = ThisPassDist;
		}
	}
	return nCurrentClusters;
}

/************************************************/

//! Refine clusters using equilibrium k-means
static void RefineClusters(
	TCluster_t *Clusters,
	const TClusterData_t *Data,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	uint32_t nPasses,
	const float *Weights,
	float alpha,
	uint32_t nDims
) {
	uint32_t n, k;
	uint32_t Pass;
	float LastPassDist = INFINITY;
	for(Pass=0;Pass<nPasses;Pass++) {
		for(k=0;k<nClusters;k++) {
			Clusters[k].FirstDataIdx = CLUSTER_END_OF_LIST;
			ClearTraining(&Clusters[k], nDims);
		}
		for(n=0;n<nDataPoints;n++) {
			const TClusterData_t *p = Data + n*nDims;
			float w_n = Weights ? Weights[n] : 1.0f;

			//! Find minimum distance to improve the LogSumExp-style summation.
			//! Note: We also store the hard assignment at this point.
			uint32_t MinIdx  = 0;
			float    MinDist = INFINITY;
			for(k=0;k<nClusters;k++) {
				float d_kn = TCluster_Dist2ToCentroid(&Clusters[k], p, nDims);
				if(d_kn < MinDist) {
					MinIdx  = k;
					MinDist = d_kn;
				}
			}
			InsertAtListHead(n, &Clusters[MinIdx].FirstDataIdx, ClusterListIndices);

			//! Let SoftminDistW = Sum[E^(-a*d_{i,n}), {i,K}]
			//! Let SoftminDist  = Sum[d_{i,n}*e^(-a*d_{i,n}), {i,K}] / SoftminDistW
			float WeightSum    = 0.0f;
			float SoftminDist  = 0.0f;
			float SoftminDistW = 0.0f;
			for(k=0;k<nClusters;k++) {
				float d_kn = TCluster_Dist2ToCentroid(&Clusters[k], p, nDims) - MinDist;
				float w = w_n * expf(-alpha * d_kn);
				SoftminDist  += w * d_kn;
				SoftminDistW += w;
			}
			if(SoftminDistW < 1.0e-10f) {
				//! Numerical collapse - perform hard assign only
				TrainCluster(&Clusters[MinIdx], p, w_n, nDims);
				continue;
			}
			float InvSoftminDistW = 1.0f / SoftminDistW;
			SoftminDist *= InvSoftminDistW;

			//! w_{k,n} = e^(-a*d_{k,n}) / SoftminDistW * (1 - a*(d_{k,n} - SoftminDist))
			for(k=0;k<nClusters;k++) {
				float d_kn = TCluster_Dist2ToCentroid(&Clusters[k], p, nDims) - MinDist;
				float w  = w_n * expf(-alpha * d_kn) * InvSoftminDistW;
				      w *= 1.0f - alpha*(d_kn - SoftminDist);
				TrainCluster(&Clusters[k], p, w, nDims);
				WeightSum += w;
			}
		}

		//! c_k = Sum[w_{k,n}*x_n, {n,N}] / Sum[w_{k,n}, {n,N}]
		float ThisPassDist = 0.0f;
		for(k=0;k<nClusters;k++) {
			ResolveCluster(&Clusters[k], Data, ClusterListIndices, Weights, nDims);
			ThisPassDist += Clusters[k].TotalDist;
		}

		//! Early exit on convergence
		if(ThisPassDist == 0.0f || ThisPassDist == LastPassDist) break;
		LastPassDist = ThisPassDist;
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
	float    Sharpness,
	const float *Weights,
	uint32_t nDims
) {
	if(!nClusters || !nDataPoints || !nDims) return 0;

	//! Build initial clusters
	float InvVariance = 0.0f; //! <- Shouldn't be needed, but gcc complains
	uint32_t nOutputClusters = CreateClusters(
		Clusters,
		Data,
		nClusters,
		nDataPoints,
		ClusterListIndices,
		(Sharpness <= 0.0f) ? nPasses : 0,
		Weights,
		&InvVariance,
		nDims
	);
	if(nOutputClusters < nClusters || Sharpness <= 0.0f) {
		//! Early convergence achieved, or no EKM desired - early exit
		return nOutputClusters;
	}

	//! Refine clusters using EKM
	RefineClusters(
		Clusters,
		Data,
		nClusters,
		nDataPoints,
		ClusterListIndices,
		nPasses,
		Weights,
		Sharpness * InvVariance,
		nDims
	);
	return nClusters;
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

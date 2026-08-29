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

//! Create initial clusters using a variation of median-cut
//! The variation is to split along the center of the principal axis
//! rather than the median of the axis of highest variance.
static uint32_t CreateInitialClusters(
	TCluster_t *Clusters,
	const TClusterData_t *Data,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	const float *Weights,
	uint32_t nDims
) {
	uint32_t n, k, This;

	//! Assign all points to cluster 0 to start with
	uint32_t nCurrentClusters = 1;
	Clusters[0].FirstDataIdx = CLUSTER_END_OF_LIST;
	ClearTraining(&Clusters[0], nDims);
	for(n=0;n<nDataPoints;n++) {
		float w = Weights ? Weights[n] : 1.0f;
		TrainCluster(&Clusters[0], Data + n*nDims, w, nDims);
		InsertAtListHead(n, &Clusters[0].FirstDataIdx, ClusterListIndices);
	}
	ResolveCluster(&Clusters[0], Data, ClusterListIndices, Weights, nDims);

	//! Now perform cuts along the principal axes of maximum variance
	if(Clusters[0].TotalDist != 0.0f) while(nCurrentClusters < nClusters) {
		//! Find the principal axis of each cluster, and save
		//! the cluster of highest variance for splitting
		uint32_t MaxDistClusterIdx = 0;
		float MaxDistClusterAvg = 0;
		float MaxDistClusterVar = 0;
		for(k=0;k<nCurrentClusters;k++) {
			//! Find principal axis
			TCluster_EigenAxis(&Clusters[k], Data, ClusterListIndices, Weights, nDims);

			//! Now find mean and variance along this axis
			float Avg = TCluster_EigenProjectCentroid(&Clusters[k], nDims);
			float Var = 0;
			for(This=Clusters[k].FirstDataIdx;This!=CLUSTER_END_OF_LIST;This=ClusterListIndices[This]) {
				float w = Weights ? Weights[This] : 1.0f;
				float p = TCluster_EigenProject(&Clusters[k], Data + This*nDims, nDims);
				float d = p - Avg;
				Var += w * d*d;
			}
			if(Var > MaxDistClusterVar) {
				MaxDistClusterIdx = k;
				MaxDistClusterAvg = Avg;
				MaxDistClusterVar = Var;
			}
		}
		if(MaxDistClusterVar == 0.0f) {
			//! No distortion remaining
			break;
		}

		//! Split cluster through the center of its principal axis
		uint32_t SrcCluster = MaxDistClusterIdx;
		uint32_t NewCluster = nCurrentClusters;
		uint32_t Next = Clusters[SrcCluster].FirstDataIdx;
		Clusters[SrcCluster].FirstDataIdx = CLUSTER_END_OF_LIST;
		Clusters[NewCluster].FirstDataIdx = CLUSTER_END_OF_LIST;
		ClearTraining(&Clusters[SrcCluster], nDims);
		ClearTraining(&Clusters[NewCluster], nDims);
		for(This=Next;This!=CLUSTER_END_OF_LIST;This=Next) {
			float p = TCluster_EigenProject(&Clusters[SrcCluster], Data + This*nDims, nDims);
			float w = Weights ? Weights[This] : 1.0f;
			uint32_t DstCluster = (p < MaxDistClusterAvg) ? SrcCluster : NewCluster;
			Next = ClusterListIndices[This];
			TrainCluster(&Clusters[DstCluster], Data + This*nDims, w, nDims);
			InsertAtListHead(This, &Clusters[DstCluster].FirstDataIdx, ClusterListIndices);
		}
		ResolveCluster(&Clusters[SrcCluster], Data, ClusterListIndices, Weights, nDims);
		ResolveCluster(&Clusters[NewCluster], Data, ClusterListIndices, Weights, nDims);
		nCurrentClusters++;
	}
	return nCurrentClusters;
}

/************************************************/

//! Create initial clusters using median-cut, then refine
//! the results using standard k-means clustering.
static uint32_t CreateClusters(
	TCluster_t *Clusters,
	const TClusterData_t *Data,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	uint32_t nPasses,
	const float *Weights,
	float *InvVariancePtr,
	uint32_t nDims
) {
	uint32_t n, k;

	//! Create initial clusters using median-cut
	uint32_t nMedianCutClusters = CreateInitialClusters(Clusters, Data, nClusters, nDataPoints, ClusterListIndices, Weights, nDims);
	if(nMedianCutClusters < nClusters) return nMedianCutClusters;

	//! Calculate initial distortions
	//! Note that median-cut guarantees that no empty clusters
	uint32_t DistClusterHead  = CLUSTER_END_OF_LIST;
	uint32_t EmptyClusterHead = CLUSTER_END_OF_LIST;
	float LastPassDist = 0.0f; {
		float w = 0.0f;
		for(k=0;k<nClusters;k++) {
			InsertToDistortedClusterList(Clusters, k, &DistClusterHead);
			LastPassDist += Clusters[k].TotalDist;
			w += Clusters[k].TrainWeight;
		}
		if(LastPassDist == 0.0f) {
			//! No distortion remaining
			*InvVariancePtr = 0;
			return nClusters;
		}
		*InvVariancePtr = w / LastPassDist;
	}

	//! Begin k-means passes
	uint8_t LastPassSameDist = 0;
	while(nPasses || (EmptyClusterHead != CLUSTER_END_OF_LIST && DistClusterHead != CLUSTER_END_OF_LIST)) {
		//! Fill empty clusters from the most distorted data points.
		//! Note that we are splitting out the most distorted point
		//! of the most distorted clusters, NOT the most distorted
		//! points general. This is an important distinction.
		while(EmptyClusterHead != CLUSTER_END_OF_LIST && DistClusterHead != CLUSTER_END_OF_LIST) {
			uint32_t SrcCluster = PopDistortedClusterList(Clusters, &DistClusterHead);
			uint32_t DstCluster = PopEmptyClusterList    (Clusters, &EmptyClusterHead);
			SplitCluster(&Clusters[DstCluster], &Clusters[SrcCluster], Data, ClusterListIndices, Weights, nDims);
		}

		//! Re-cluster the data points
		if(nPasses != 0) {
			//! Assign points to clusters
			for(k=0;k<nClusters;k++) {
				Clusters[k].FirstDataIdx = CLUSTER_END_OF_LIST;
				ClearTraining(&Clusters[k], nDims);
			}
			for(n=0;n<nDataPoints;n++) {
				uint32_t BestIdx  = 0;
				float    BestDist = INFINITY;
				for(k=0;k<nClusters;k++) {
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
			nPasses--;
		}

		//! Resolve new clusters
		float ThisPassDist = 0.0f; {
			DistClusterHead  = CLUSTER_END_OF_LIST;
			EmptyClusterHead = CLUSTER_END_OF_LIST;
			for(k=0;k<nClusters;k++) {
				if(ResolveCluster(&Clusters[k], Data, ClusterListIndices, Weights, nDims)) {
					//! If the cluster resolves, update the distortion linked list
					InsertToDistortedClusterList(Clusters, k, &DistClusterHead);
					ThisPassDist += Clusters[k].TotalDist;
				} else {
					//! No resolve - append to empty-cluster linked list
					InsertToEmptyClusterList(Clusters, k, &EmptyClusterHead);
				}
			}
		}

		//! Early exit on convergence
		if(ThisPassDist == 0.0f) break;
		if(ThisPassDist == LastPassDist) {
			//! If distortion matches the last pass, we try once
			//! more just in case this was due to rounding error
			if(LastPassSameDist++ != 0) break;
		} else LastPassSameDist = 0;
		LastPassDist = ThisPassDist;
	}
	return nClusters;
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
	uint8_t LastPassSameDist = 0;
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
			}
		}

		//! c_k = Sum[w_{k,n}*x_n, {n,N}] / Sum[w_{k,n}, {n,N}]
		float ThisPassDist = 0.0f;
		for(k=0;k<nClusters;k++) {
			ResolveCluster(&Clusters[k], Data, ClusterListIndices, Weights, nDims);
			ThisPassDist += Clusters[k].TotalDist;
		}

		//! Early exit on convergence
		if(ThisPassDist == 0.0f) break;
		if(ThisPassDist == LastPassDist) {
			//! If distortion matches the last pass, we try once
			//! more just in case this was due to rounding error
			if(LastPassSameDist++ != 0) break;
		} else LastPassSameDist = 0;
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
	//! Note that EKM appears to work best if we feed it
	//! the median-cut data directly rather than running
	//! k-means on it first!
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
	if(InvVariance == 0.0f || Sharpness <= 0.0f) {
		//! No distortion remaining, or no EKM desired - early exit
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

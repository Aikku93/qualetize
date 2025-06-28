/************************************************/
#include <math.h>
/************************************************/
#include "Cluster.h"
/************************************************/

typedef        float     TClusterData_t;
typedef struct Cluster_t TCluster_t;

/************************************************/

static float TCluster_Dist2ToCentroid(TCluster_t *x, const TClusterData_t *Data, uint32_t nDims) {
	uint32_t n;
	float Dist2 = 0.0f;
	for(n=0;n<nDims;n++) {
		float d = Data[n] - x->Centroid[n];
		Dist2 += d*d;
	}
	return Dist2;
}

static void TCluster_ClearTrainingVector(TCluster_t *x, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) x->Training[n] = 0.0f;
}

static void TCluster_AddToTraining(TCluster_t *x, const TClusterData_t *Data, float w, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) {
		x->Training[n] += Data[n]*w;
	}
}

static void TCluster_ResolveCentroid(TCluster_t *x, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) {
		x->Centroid[n] = x->Training[n] / x->TrainWeight;
	}
}

static void TCluster_SetCentroidToData(TCluster_t *x, const TClusterData_t *Data, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) {
		x->Centroid[n] = Data[n];
	}
}

/************************************************/
#include "Cluster-Template.h"
/************************************************/

//! Apply cluster analysis to the specified data
uint32_t Clusterize_Process(
	struct Cluster_t *Clusters,
	const float *Data,
	uint32_t nDims,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	uint32_t nPasses,
	float    SplitRatio,
	const float *Weights
) {
	return TClusterize_Process(
		Clusters,
		Data,
		nClusters,
		nDataPoints,
		ClusterListIndices,
		nPasses,
		SplitRatio,
		Weights,
		nDims
	);
}

/************************************************/

//! Get cluster indices from clusters and data lists
void Clusterize_GetClusterIndices_u8(
	uint8_t *DataClusterIndices,
	const struct Cluster_t *Clusters,
	uint32_t nClusters,
	const uint32_t *ClusterListIndices
) {
	TClusterize_GetClusterIndices_u8(
		DataClusterIndices,
		Clusters,
		nClusters,
		ClusterListIndices
	);
}

/************************************************/
//! EOF
/************************************************/

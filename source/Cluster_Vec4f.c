/************************************************/
#include <math.h>
/************************************************/
#include "Cluster.h"
#include "Vec4f.h"
/************************************************/

typedef        Vec4f_t         TClusterData_t;
typedef struct Cluster_Vec4f_t TCluster_t;

/************************************************/

static float TCluster_Dist2ToCentroid(TCluster_t *x, const TClusterData_t *Data, uint32_t nDims) {
	(void)nDims;
	return Vec4f_Dist2(Data, &x->Centroid);
}

static void TCluster_ClearTrainingVector(TCluster_t *x, uint32_t nDims) {
	(void)nDims;
	x->Training = VEC4F_EMPTY;
}

static void TCluster_AddToTraining(TCluster_t *x, const TClusterData_t *Data, float w, uint32_t nDims) {
	(void)nDims;
	Vec4f_t t = Vec4f_Muli(Data, w);
	x->Training = Vec4f_Add(&x->Training, &t);
}

static void TCluster_ResolveCentroid(TCluster_t *x, uint32_t nDims) {
	(void)nDims;
	x->Centroid = Vec4f_Divi(&x->Training, x->TrainWeight);
}

static void TCluster_SetCentroidToData(TCluster_t *x, const TClusterData_t *Data, uint32_t nDims) {
	(void)nDims;
	x->Centroid = *Data;
}

/************************************************/
#include "Cluster-Template.h"
/************************************************/

//! Apply cluster analysis to the specified data
uint32_t Clusterize_Vec4f_Process(
	struct Cluster_Vec4f_t *Clusters,
	const Vec4f_t *Data,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	uint32_t nPasses,
	const float *Weights
) {
	return TClusterize_Process(
		Clusters,
		Data,
		nClusters,
		nDataPoints,
		ClusterListIndices,
		nPasses,
		Weights,
		1
	);
}

/************************************************/

//! Get cluster indices from clusters and data lists
void Clusterize_Vec4f_GetClusterIndices_u8(
	uint8_t *DataClusterIndices,
	const struct Cluster_Vec4f_t *Clusters,
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

/************************************************/
#include <math.h>
/************************************************/
#include "Cluster.h"
/************************************************/

typedef        float     TClusterData_t;
typedef        float*    TClusterVector_t;
typedef struct Cluster_t TCluster_t;

/************************************************/
#if CLUSTER_USE_COMPENSATED_SUMMATION
/************************************************/

//! Kahan summation
static inline float KahanSum(float Sum, float Input, float *c) {
	float y = Input - *c;
	float t = Sum + y;
	*c = t - Sum - y;
	return t;
}

/************************************************/
#endif
/************************************************/

//! Standard vector operations
static void VecClear(TClusterVector_t *x, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) (*x)[n] = 0.0f;
}
static void VecCopy(TClusterVector_t *Dst, const TClusterVector_t *Src, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) (*Dst)[n] = (*Src)[n];
}
static void VecSet(TClusterVector_t *Dst, const TClusterData_t *Src, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) (*Dst)[n] = Src[n];
}
#if CLUSTER_USE_COMPENSATED_SUMMATION
static void VecSum(TClusterVector_t *Sum, const TClusterData_t *Input, TClusterVector_t *c, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) (*Sum)[n] = KahanSum((*Sum)[n], Input[n], &(*c)[n]);
}
#else
static void VecSum(TClusterVector_t *Sum, const TClusterData_t *Input, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) (*Sum)[n] += Input[n];
}
#endif
static void VecDivi(TClusterVector_t *x, float y, uint32_t nDims) {
	uint32_t n;
	for(n=0;n<nDims;n++) (*x)[n] /= y;
}
static float VecDist2(const TClusterData_t *a, const TClusterVector_t *b, uint32_t nDims) {
	uint32_t n;
	float Dist = 0.0f;
	for(n=0;n<nDims;n++) {
		float d = a[n] - (*b)[n];
		Dist += d*d;
	}
	return Dist;
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
	float AvgErrorThreshold
) {
	return TClusterize_Process(
		Clusters,
		Data,
		nClusters,
		nDataPoints,
		ClusterListIndices,
		nPasses,
		AvgErrorThreshold,
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

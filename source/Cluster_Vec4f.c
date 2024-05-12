/************************************************/
#include <math.h>
/************************************************/
#include "Cluster.h"
#include "Vec4f.h"
/************************************************/

typedef        Vec4f_t         TClusterData_t;
typedef        Vec4f_t         TClusterVector_t;
typedef struct Cluster_Vec4f_t TCluster_t;

/************************************************/
#if CLUSTER_USE_COMPENSATED_SUMMATION
/************************************************/

//! Kahan summation
static inline Vec4f_t Vec4f_KahanSum(const Vec4f_t *Sum, const Vec4f_t *Input, Vec4f_t *c) {
	Vec4f_t y = Vec4f_Sub(Input, c);
	Vec4f_t t = Vec4f_Add(Sum, &y);
	       *c = Vec4f_Sub(&t, Sum);
	       *c = Vec4f_Sub(c, &y);
	return t;
}

/************************************************/
#endif
/************************************************/

//! Standard vector operations
static void VecClear(TClusterVector_t *x, uint32_t nDims) {
	(void)nDims;
	*x = VEC4F_EMPTY;
}
static void VecCopy(TClusterVector_t *Dst, const TClusterVector_t *Src, uint32_t nDims) {
	(void)nDims;
	*Dst = *Src;
}
static void VecSet(TClusterVector_t *Dst, const TClusterData_t *Src, uint32_t nDims) {
	(void)nDims;
	*Dst = *Src;
}
#if CLUSTER_USE_COMPENSATED_SUMMATION
static void VecSum(TClusterVector_t *Sum, const TClusterData_t *Input, TClusterVector_t *c, uint32_t nDims) {
	(void)nDims;
	*Sum = Vec4f_KahanSum(Sum, Input, c);
}
#else
static void VecSum(TClusterVector_t *Sum, const TClusterData_t *Input, uint32_t nDims) {
	(void)nDims;
	*Sum = Vec4f_Add(Sum, Input);
}
#endif
static void VecDivi(TClusterVector_t *x, float y, uint32_t nDims) {
	(void)nDims;
	*x = Vec4f_Divi(x, y);
}
static float VecDist2(const TClusterData_t *a, const TClusterVector_t *b, uint32_t nDims) {
	(void)nDims;
	return Vec4f_Dist2(a, b);
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

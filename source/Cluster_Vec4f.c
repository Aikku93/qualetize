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

static void TCluster_EigenAxis(TCluster_t *x, const TClusterData_t *Data, const uint32_t *ClusterListIndices, const float *Weights, uint32_t nDims) {
#define FLATTEN_LTRI(i,j) ((i)*((i)+1)/2+(j))
	(void)nDims;
	uint32_t i, j, k, m, Next;

	//! Calculate covariance matrix
	//! Because the matrix is symmetric, we keep only the
	//! lower-triangular part, which is N*(N+1)/2 entries
	float K[FLATTEN_LTRI(4, 0)];
	for(m=0;m<FLATTEN_LTRI(4,0);m++) K[m] = 0.0f;
	for(Next=x->FirstDataIdx;Next!=CLUSTER_END_OF_LIST;Next=ClusterListIndices[Next]) {
		float w = Weights ? Weights[Next] : 1.0f;
		const Vec4f_t p = Vec4f_Sub(&Data[Next], &x->Centroid);
		for(i=0;i<4;i++) for(j=0;j<=i;j++) {
			float pi = p.f32[i];
			float pj = p.f32[j];
			K[FLATTEN_LTRI(i,j)] += w * (pi*pj);
		}
	}

	//! Use power iteration to find the principal axis
	Vec4f_t *Axis = &x->Axis;
	Vec4f_t  Temp;
	*Axis = Vec4f_Broadcast(1.0f);
	for(k=0;k<EIGEN_ITERS;k++) {
		Temp = Vec4f_Broadcast(0.0f);
		for(i=0;i<4;i++) for(j=0;j<4;j++) {
			Temp.f32[i] += K[(i >= j) ? FLATTEN_LTRI(i,j) : FLATTEN_LTRI(j,i)]*Axis->f32[j];
		}
		float Norm2 = Vec4f_Length2(&Temp);
		if(Norm2 == 0.0f) break;
		float InvNorm = 1.0f / sqrtf(Norm2);
		*Axis = Vec4f_Muli(&Temp, InvNorm);
	}
#undef FLATTEN_LTRI
}

static float TCluster_EigenProject(TCluster_t *x, const TClusterData_t *Data, uint32_t nDims) {
	(void)nDims;
	return Vec4f_Dot(&x->Axis, Data);
}

static float TCluster_EigenProjectCentroid(TCluster_t *x, uint32_t nDims) {
	(void)nDims;
	return Vec4f_Dot(&x->Axis, &x->Centroid);
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
	float    Sharpness,
	const float *Weights
) {
	return TClusterize_Process(
		Clusters,
		Data,
		nClusters,
		nDataPoints,
		ClusterListIndices,
		nPasses,
		Sharpness,
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

/************************************************/
#include <math.h>
#include <stdlib.h>
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

static void TCluster_EigenAxis(TCluster_t *x, const TClusterData_t *Data, const uint32_t *ClusterListIndices, const float *Weights, uint32_t nDims) {
#define FLATTEN_LTRI(i,j) ((i)*((i)+1)/2+(j))
	uint32_t i, j, k, m, Next;

	//! Calculate covariance matrix
	//! Because the matrix is symmetric, we keep only the
	//! lower-triangular part, which is N*(N+1)/2 entries
	float *K = (float*)alloca(sizeof(float[FLATTEN_LTRI(nDims, 0)]));
	for(m=0;m<FLATTEN_LTRI(nDims,0);m++) K[m] = 0.0f;
	for(Next=x->FirstDataIdx;Next!=CLUSTER_END_OF_LIST;Next=ClusterListIndices[Next]) {
		float w = Weights ? Weights[Next] : 1.0f;
		const float *p = Data + Next*nDims;
		for(i=0;i<nDims;i++) for(j=0;j<=i;j++) {
			float pi = p[i] - x->Centroid[i];
			float pj = p[j] - x->Centroid[j];
			K[FLATTEN_LTRI(i,j)] += w * (pi*pj);
		}
	}

	//! Use power iteration to find the principal axis
	float *Axis = x->Axis;
	float *Temp = (float*)alloca(sizeof(float[nDims]));
	for(m=0;m<nDims;m++) Axis[m] = 1.0f;
	for(k=0;k<EIGEN_ITERS;k++) {
		for(m=0;m<nDims;m++) Temp[m] = 0;
		for(i=0;i<nDims;i++) for(j=0;j<nDims;j++) {
			Temp[i] += K[(i >= j) ? FLATTEN_LTRI(i,j) : FLATTEN_LTRI(j,i)]*Axis[j];
		}
		float Norm2 = 0;
		for(m=0;m<nDims;m++) Norm2 += Temp[m]*Temp[m];
		if(Norm2 == 0.0f) break;
		float InvNorm = 1.0f / sqrtf(Norm2);
		for(m=0;m<nDims;m++) Axis[m] = Temp[m]*InvNorm;
	}
#undef FLATTEN_LTRI
}

static float TCluster_EigenProject(TCluster_t *x, const TClusterData_t *Data, uint32_t nDims) {
	uint32_t n;
	float p = 0.0f;
	for(n=0;n<nDims;n++) {
		p += x->Axis[n] * Data[n];
	}
	return p;
}

static float TCluster_EigenProjectCentroid(TCluster_t *x, uint32_t nDims) {
	uint32_t n;
	float p = 0.0f;
	for(n=0;n<nDims;n++) {
		p += x->Axis[n] * x->Centroid[n];
	}
	return p;
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

/************************************************/
#pragma once
/************************************************/
#include <stdint.h>
/************************************************/
#include "Vec4f.h"
/************************************************/

//! "End of list" signal in linked lists
#define CLUSTER_END_OF_LIST (~0u)

/************************************************/

//! Generic cluster
//! NOTE: Centroid, Training, and cTraining must be assigned
//! by the caller; the functions assume that these pointers
//! have already been set. Each vector will store N values,
//! where N is the number of dimensions.
struct Cluster_t {
	uint32_t nPoints;        //! Number of points assigned to cluster
	uint32_t NextCluster;    //! Next cluster in list
	uint32_t FirstDataIdx;   //! First data point assigned to this cluster
	uint32_t MaxDistIdx;     //! Index of the most distorted data point
	float    TotalDist;      //! Sum of distortion of all data points in this cluster
	float   *Centroid;       //! Centroid of cluster
	float   *Training;       //! Training summation
};

//! Specialized Vec4f_t cluster
//! Because SSE-enhanced Vec4f_t requires alignment, then so too
//! does Cluster_Vec4f_t.
struct
#ifdef __SSE__
__attribute__((aligned(16)))
#endif
Cluster_Vec4f_t {
	uint32_t nPoints;
	uint32_t NextCluster;
	uint32_t FirstDataIdx;
	uint32_t MaxDistIdx;
	float    TotalDist;
	Vec4f_t  Centroid;
	Vec4f_t  Training;
};

/************************************************/

//! Apply cluster analysis to the specified data
//! Returns the number of used clusters.
//! Notes:
//!  -Data[] must be composed of sequential vectors of length nDims;
//!   eg. `float Data[nDataPoints][nDims]`.
//!  -ClusterListIndices[] is used to accelerate searches for all data
//!   points in the specified cluster, and must be nDataPoints in size.
uint32_t Clusterize_Process(
	struct Cluster_t *Clusters,
	const float *Data,
	uint32_t nDims,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	uint32_t nPasses
);
uint32_t Clusterize_Vec4f_Process(
	struct Cluster_Vec4f_t *Clusters,
	const Vec4f_t *Data,
	uint32_t nClusters,
	uint32_t nDataPoints,
	uint32_t *ClusterListIndices,
	uint32_t nPasses
);

//! Get cluster indices from clusters and data lists (uint8_t)
//! Notes:
//!  -ClusterListIndices[] should be input from Clusterize_Process().
void Clusterize_GetClusterIndices_u8(
	uint8_t *DataClusterIndices,
	const struct Cluster_t *Clusters,
	uint32_t nClusters,
	const uint32_t *ClusterListIndices
);
void Clusterize_Vec4f_GetClusterIndices_u8(
	uint8_t *DataClusterIndices,
	const struct Cluster_Vec4f_t *Clusters,
	uint32_t nClusters,
	const uint32_t *ClusterListIndices
);

/************************************************/
//! EOF
/************************************************/

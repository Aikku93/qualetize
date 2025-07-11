/************************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
/************************************************/
#include "Bitmap.h"
#include "Cluster.h"
#include "Qualetize.h"
#include "Qualetize-Colourspace.h"
#include "Vec4f.h"
/************************************************/
#define CLAMP(x, Min, Max) ((x) < (Min) ? (Min) : (x) > (Max) ? (Max) : (x))
/************************************************/

//! Tiles are clustered in eight dimensions: Four channels each
//! for the mean, and another four for the variances.
#define TILE_CLUSTER_DIMENSIONS 8

//! Default number of clustering passes
//! When using these values, these refer to the number of passes over
//! all the clusters, and the actual number of passes per cluster will
//! be equal to DEFAULT_PASSES/nClusters; this keeps the running time
//! more manageable.
#define DEFAULT_TILECLUSTER_PASSES   1024
#define DEFAULT_COLOURCLUSTER_PASSES 1024

//! Number of passes to use for palette sorting
#define PALETTE_SORT_PASSES 32

/************************************************/

//! Convert BGRA8 to native Vec4f type in RGBA format
static Vec4f_t BGRA8_To_Vec4fRGBA(const BGRA8_t *x) {
	Vec4f_t y;
	y.f32[0] = (float)x->r;
	y.f32[1] = (float)x->g;
	y.f32[2] = (float)x->b;
	y.f32[3] = (float)x->a;
	return Vec4f_Divi(&y, 255.0f);
}

//! Convert from Vec4f RGBA to BGRA8
static BGRA8_t Vec4fRGBA_To_BGRA8(const Vec4f_t *x) {
	Vec4f_t yf = Vec4f_Clamp(&yf, 0.0f, 1.0f);
	        yf = Vec4f_Muli(x, 255.0f);
	        yf = Vec4f_Round(&yf);
	BGRA8_t y;
	y.r = (uint8_t)yf.f32[0];
	y.g = (uint8_t)yf.f32[1];
	y.b = (uint8_t)yf.f32[2];
	y.a = (uint8_t)yf.f32[3];
	return y;
}

/************************************************/

//! Calculate the colour value of a tile
uint8_t CalculateTileColourValue(
	float *TileValue,
	float *TileWeight,
	const Vec4f_t *InputPxData,
	uint32_t TileX,
	uint32_t TileY,
	uint32_t InputWidth,
	uint32_t InputHeight,
	const struct QualetizePlan_t *Plan
) {
	//! Convert colours to Lab space, find the mean, and check opacity
	//! CIELab and Oklab work fine here, but have different strengths.
	//! In general, Oklab seems to preserve smooth gradients better,
	//! but CIELab preserves colour information a lot better.
	uint32_t x, y;
	uint32_t nSampledPoints = 0;
	Vec4f_t  Mean  = VEC4F_EMPTY, Var = VEC4F_EMPTY;
	Vec4f_t wMean  = VEC4F_EMPTY;
	float   wMeanW = 0.0f;
	for(y=0;y<Plan->TileHeight;y++) for(x=0;x<Plan->TileWidth;x++) {
		uint32_t vx = TileX*Plan->TileWidth  + x;
		uint32_t vy = TileY*Plan->TileHeight + y;
		if(vx < InputWidth && vy < InputHeight) {
			uint32_t PxOffs = vy*InputWidth + vx;
			Vec4f_t Px = ConvertToColourspace(&InputPxData[PxOffs], COLOURSPACE_CIELAB);
			if(!Plan->FirstColourIsTransparent || Px.f32[3] != 0.0f) {
				float w = Px.f32[0] * Px.f32[3];
				Vec4f_t Px2 = Vec4f_Mul(&Px, &Px);
				Vec4f_t wPx = Vec4f_Muli(&Px, w);
				Mean  = Vec4f_Add(&Mean,  &Px);
				Var   = Vec4f_Add(&Var,   &Px2);
				wMean = Vec4f_Add(&wMean, &wPx);
				wMeanW += w;
				nSampledPoints++;
			}
		}
	}

	//! If tile is fully transparent AND we have a forced transparent
	//! palette colour, then we can skip processing this tile
	if(Plan->FirstColourIsTransparent && !nSampledPoints) return 0;

	//! Normalize means and variance, calculate standard deviation
	Vec4f_t Mean2, Dev;
	static const Vec4f_t Zero = VEC4F_EMPTY;
	Mean  = Vec4f_Divi(&Mean, (float)nSampledPoints);
	Var   = Vec4f_Divi(&Var,  (float)nSampledPoints);
	wMean = (wMeanW != 0.0f) ? Vec4f_Divi(&wMean, wMeanW) : Mean;
	Mean2 = Vec4f_Mul(&Mean, &Mean);
	Var   = Vec4f_Sub(&Var, &Mean2);
	Var   = Vec4f_Max(&Var, &Zero); //! <- Protect against round-off error before square root
	Dev   = Vec4f_Sqrt(&Var);

	//! Calculate saturation of the weighted-mean colours.
	//! The colours were weighted by luma and alpha, such that we
	//! give more importance to brighter and more opaque colours
	//! than darker and/or translucent colours.
	//! We use the saturation as a cluster weight, because we want
	//! less-saturated clusters to get fewer palettes, as we can
	//! instead re-distribute them to more colourful tiles, which
	//! should result in less averaging of very different colours,
	//! and thus reduce the amount of desaturation that happens.
	//! Also note that we weight the saturation again, by the luma
	//! component - a colourful-but-dark tile is unimportant.
	float Luma2   = wMean.f32[0]*wMean.f32[0];
	float Chroma2 = wMean.f32[1]*wMean.f32[1] + wMean.f32[2]*wMean.f32[2];
	float Satur   = (Luma2 > 0.0f) ? sqrtf(Chroma2 / (Luma2 + Chroma2)) : 0.0f;

	//! Fill the tile data.
	//! Note that we further weight by the average alpha (not the
	//! weighted average alpha!), since low-alpha areas can be
	//! banded together into fewer palettes with few issues.
	//! Finally, we combine saturation and the amount of deviation
	//! into a single weight quantity: we want saturated tiles to
	//! split into more palettes as needed, and we also want tiles
	//! with high variance to split into more palettes; for tiles
	//! that are saturated and have high variance, we want these
	//! to have higher priority for splitting.
	//! Note that this latter saturation is itself weighted by the
	//! weighted average luma, such that brighter/more prominent
	//! colours are favoured with higher weights, because we want
	//! to penalize the opposite condition (colourful and/or varied
	//! /dark/ colours, where we can afford to ignore distortion).
	TileValue[0] = 0.25f*wMean.f32[0];
	TileValue[1] = wMean.f32[1];
	TileValue[2] = wMean.f32[2];
	TileValue[3] = wMean.f32[3];
	TileValue[4] = 0.25f*Dev.f32[0];
	TileValue[5] = Dev.f32[1];
	TileValue[6] = Dev.f32[2];
	TileValue[7] = Dev.f32[3];
	*TileWeight = Mean.f32[3]*((1.0f - wMean.f32[0]) + wMean.f32[0]*(1.0f+Satur));
	return 1;
}

/************************************************/

//! Compare hue (assumes that the hue value is stored in the first element)
static int CompareHueAscending(const void *a_, const void *b_) {
	const struct Cluster_Vec4f_t *a = (const struct Cluster_Vec4f_t*)a_;
	const struct Cluster_Vec4f_t *b = (const struct Cluster_Vec4f_t*)b_;
	return (a->Centroid.f32[0] < b->Centroid.f32[0]) ? (-1) : (+1);
}

//! Compare luma
static inline float LumaFromRGB(const Vec4f_t *x) {
	//! This just uses the luma calculation from YCbCr, in linear RGB
	float r = (x->f32[0] > 0.0f) ? powf(x->f32[0], 2.4f) : 0.0f;
	float g = (x->f32[1] > 0.0f) ? powf(x->f32[1], 2.4f) : 0.0f;
	float b = (x->f32[2] > 0.0f) ? powf(x->f32[2], 2.4f) : 0.0f;
	return cbrtf(0.2126f*r + 0.71520f*g + 0.0722f*b) * x->f32[3];
}
static int CompareLumaAscending(const void *a_, const void *b_) {
	const Vec4f_t *a = (const Vec4f_t*)a_;
	const Vec4f_t *b = (const Vec4f_t*)b_;
	return (LumaFromRGB(a) < LumaFromRGB(b)) ? (-1) : (+1);
}
static int CompareLumaDescending(const void *a_, const void *b_) {
	const Vec4f_t *a = (const Vec4f_t*)a_;
	const Vec4f_t *b = (const Vec4f_t*)b_;
	return (LumaFromRGB(a) > LumaFromRGB(b)) ? (-1) : (+1);
}

//! Apply qualetize operation to input image using the specified plan
uint8_t Qualetize(
	      uint8_t  *OutputPxData,
	      BGRA8_t  *OutputPalette,
	const void     *InputBitmap,
	const BGRA8_t  *InputPalette,
	      uint32_t  InputWidth,
	      uint32_t  InputHeight,
	const struct QualetizePlan_t *Plan,
	Vec4f_t *RMSE
) {
	uint32_t nTilesX     = (InputWidth  - 1) / Plan->TileWidth  + 1;
	uint32_t nTilesY     = (InputHeight - 1) / Plan->TileHeight + 1;
	uint32_t nTilesTotal = nTilesX * nTilesY;
	uint32_t nPxPerTile  = Plan->TileWidth * Plan->TileHeight;
	uint32_t nPxTotal    = nTilesTotal * nPxPerTile;

	//! Original image data, used for final colour matching.
	//! Contains ImageWidth*ImageHeight elements.
	Vec4f_t *InputPxData;

	//! Pixels for tiles.
	//! Contains ImageWidth*ImageHeight elements.
	//! Each tile takes TileWidth*TileHeight elements, and the tiles are then stored
	//! sequentially, in raster-scan order for simplicity.
	Vec4f_t *TilePxData;

	//! Palette data.
	//! Contains nPaletteColours*nTilePalettes elements.
	//! This is needed to contain the perceptual-space colours before writing
	//! to the output BGRA8 palette destination.
	Vec4f_t *PaletteData;

	//! Palette indices for tiles.
	//! Contains (ImageWidth/TileWidth)*(ImageHeight/TileHeight) elements.
	//! This data is used in the final clustering step to build the final pixel
	//! palette indices by combining with the clustered colour from the palette.
	uint8_t *TilePaletteIndices;

	//! Colour data for tile clustering.
	//! Contains (ImageWidth/TileWidth)*(ImageHeight/TileHeight) elements.
	//! These values are used for clustering the tiles into their appropriate
	//! palettes. Note that these are NOT Vec4f_t, but rather float[N] to
	//! instead store different attributes to optimize.
	float *TileColourValues;

	//! Weights for tile clustering.
	//! Contains (ImageWidth/TileWidth)*(ImageHeight/TileHeight) elements.
	float *TileClusterWeights;

	//! Pixel data for cluster processing (used to store lists of colours to cluster).
	//! Contains up to ImageWidth*ImageHeight elements.
	//! The reason we need so many elements is because if all pixels of the image
	//! have been assigned to a single palette, then we need to bring all of the
	//! original pixels in for clustering.
	Vec4f_t *ColourClusterBuffer;

	//! Tile clusters for clustering tiles together.
	//! Contains nTilePalettes elements, plus extra memory for storing the cluster
	//! vectors {Centroid,Training}.
	//! ---
	//! This pointer is aliased by ColourClusters, explained below:
	//! ---
	//! Colour clusters for clustering tile colours together into a palette.
	//! Contains nPaletteColours elements.
	struct Cluster_t       *TileClusters;
	struct Cluster_Vec4f_t *ColourClusters;

	//! Cluster indices for linked-list offsets
	//! Contains up to ImageWidth*ImageHeight elements.
	uint32_t *ClusterListIndices;

	//! Allocate buffers
	void *AllocBuffer; {
		//! Get buffer offsets and allocation size
		uint32_t AllocSize = 0;
		uint32_t TotalPalCols = Plan->nPaletteColours * Plan->nTilePalettes;
		uint32_t ClustersDataSize; {
			//! NOTE: We have to manually allocate space for Centroid[], Training[]
			uint32_t ColourClusterSize = Plan->nPaletteColours * sizeof(struct Cluster_Vec4f_t);
			uint32_t TilesClusterSize  = Plan->nTilePalettes   * sizeof(struct Cluster_t);
			         TilesClusterSize += Plan->nTilePalettes   * sizeof(float)*TILE_CLUSTER_DIMENSIONS;
			         TilesClusterSize += Plan->nTilePalettes   * sizeof(float)*TILE_CLUSTER_DIMENSIONS;
			if(ColourClusterSize > TilesClusterSize) {
				ClustersDataSize = ColourClusterSize;
			} else {
				ClustersDataSize = TilesClusterSize;
			}
		}
#ifdef __SSE__
# define CREATE_BUFFER_ALIGN(x) (((x) + 15) &~ 15) //! <- Align to 16 bytes
#else
# define CREATE_BUFFER_ALIGN(x) (((x) +  3) &~  3) //! <- Align to 4 bytes (shouldn't be needed, though)
#endif
#define CREATE_BUFFER(Name, Sz) uint32_t Name##_Offs = AllocSize; AllocSize += (uint32_t)CREATE_BUFFER_ALIGN(Sz)
		CREATE_BUFFER(InputPxData,         nPxTotal     * sizeof(Vec4f_t));
		CREATE_BUFFER(TilePxData,          nPxTotal     * sizeof(Vec4f_t));
		CREATE_BUFFER(PaletteData,         TotalPalCols * sizeof(Vec4f_t));
		CREATE_BUFFER(TilePaletteIndices,  nTilesTotal  * sizeof(uint8_t));
		CREATE_BUFFER(TileColourValues,    nTilesTotal  * sizeof(float)*TILE_CLUSTER_DIMENSIONS);
		CREATE_BUFFER(TileClusterWeights,  nTilesTotal  * sizeof(float));
		CREATE_BUFFER(ColourClusterBuffer, nPxTotal     * sizeof(Vec4f_t));
		CREATE_BUFFER(ClustersData,        ClustersDataSize);
		CREATE_BUFFER(ClusterListIndices,  nPxTotal     * sizeof(uint32_t));
#undef CREATE_BUFFER
#undef CREATE_BUFFER_ALIGN

		//! Allocate buffer space
#ifdef __SSE__
		AllocSize += 16-1; //! <- For forced alignment to 128-bit lanes (16 bytes)
#endif
		AllocBuffer = malloc(AllocSize);
		if(!AllocBuffer) return 0;
		uintptr_t AllocBufferBase = (uintptr_t)AllocBuffer;
#ifdef __SSE__
		AllocBufferBase = (AllocBufferBase+15) &~ 15; //! <- Align to 16 bytes
#endif

		//! Assign pointers
#define ASSIGN_BUFFER_FROM(Name, Type, Offs) Name = (Type*)(AllocBufferBase + (Offs))
#define ASSIGN_BUFFER(Name, Type) ASSIGN_BUFFER_FROM(Name, Type, Name##_Offs)
		ASSIGN_BUFFER(InputPxData, Vec4f_t);
		ASSIGN_BUFFER(TilePxData,  Vec4f_t);
		ASSIGN_BUFFER(PaletteData, Vec4f_t);
		ASSIGN_BUFFER(TilePaletteIndices,  uint8_t);
		ASSIGN_BUFFER(TileColourValues,    float);
		ASSIGN_BUFFER(TileClusterWeights,  float);
		ASSIGN_BUFFER(ColourClusterBuffer, Vec4f_t);
		ASSIGN_BUFFER_FROM(TileClusters,   struct Cluster_t,       ClustersData_Offs);
		ASSIGN_BUFFER_FROM(ColourClusters, struct Cluster_Vec4f_t, ClustersData_Offs);
		ASSIGN_BUFFER(ClusterListIndices,  uint32_t);
#undef ASSIGN_BUFFER
#undef ASSIGN_BUFFER_FROM
	}

	//! Read the image data into the native format and write tile pixels
	//! NOTE: BlankTilesIndices aliases near the end of TileColourValues[];
	//! this is fine, however, as that memory won't be reached.
	uint32_t *BlankTileIndices = (uint32_t*)(TileColourValues + nTilesTotal*TILE_CLUSTER_DIMENSIONS);
	const uint32_t *BlankTileIndicesEnd = BlankTileIndices;
	uint32_t nNonBlankTiles = 0; {
		//! Convert raw pixels
		uint32_t x, y;
		if(InputPalette) {
			const uint8_t *Src = (const uint8_t*)InputBitmap;
			for(y=0;y<InputHeight;y++) for(x=0;x<InputWidth;x++) {
				uint8_t Idx = Src[y*InputWidth+x];
				BGRA8_t Col = InputPalette[Idx];
				Col.a = 255;
				if(
					(Plan->FirstColourIsTransparent && Idx == 0) ||
					(
						Col.r == Plan->TransparentColour.r &&
						Col.g == Plan->TransparentColour.g &&
						Col.b == Plan->TransparentColour.b
					)
				) {
					Col = (BGRA8_t){0,0,0,0};
				}
				InputPxData[y*InputWidth+x] = BGRA8_To_Vec4fRGBA(&Col);
			}
		} else {
			const BGRA8_t *Src = (const BGRA8_t*)InputBitmap;
			for(y=0;y<InputHeight;y++) for(x=0;x<InputWidth;x++) {
				BGRA8_t Col = Src[y*InputWidth+x];
				if(Plan->TransparentColour.a != 0) {
					if(
						Col.r == Plan->TransparentColour.r &&
						Col.g == Plan->TransparentColour.g &&
						Col.b == Plan->TransparentColour.b
					) Col = (BGRA8_t){0,0,0,0};
				}
				InputPxData[y*InputWidth+x] = BGRA8_To_Vec4fRGBA(&Col);
			}
		}

		//! If alpha is not pre-multiplied, apply it now
		if(!Plan->PremultipliedAlpha) {
			uint32_t n;
			for(n=0;n<nPxTotal;n++) {
				float a = InputPxData[n].f32[3];
				InputPxData[n].f32[0] *= a;
				InputPxData[n].f32[1] *= a;
				InputPxData[n].f32[2] *= a;
			}
		}

		//! Chop up the image into tiles and store their pixels
		uint32_t tx, ty;
		for(ty=0;ty<nTilesY;ty++) for(tx=0;tx<nTilesX;tx++) {
			uint32_t TileIdx = ty*nTilesX + tx;
			Vec4f_t *ThisTilePxData = TilePxData + TileIdx*nPxPerTile;

			//! Calculate this tile's value (for tile clustering)
			if(CalculateTileColourValue(
				TileColourValues + nNonBlankTiles*TILE_CLUSTER_DIMENSIONS,
				&TileClusterWeights[TileIdx],
				InputPxData,
				tx,
				ty,
				InputWidth,
				InputHeight,
				Plan
			)) {
				nNonBlankTiles++;
			} else {
				*--BlankTileIndices = TileIdx;
			}

			//! Store tile pixels
			uint32_t nPxWritten = 0;
			for(y=0;y<Plan->TileHeight;y++) for(x=0;x<Plan->TileWidth;x++) {
				uint32_t vx = tx*Plan->TileWidth  + x;
				uint32_t vy = ty*Plan->TileHeight + y;
				if(vx < InputWidth && vy < InputHeight) {
					uint32_t PxOffs = vy*InputWidth + vx;
					Vec4f_t Px = ConvertToColourspace(&InputPxData[PxOffs], Plan->Colourspace);
					ThisTilePxData[nPxWritten++] = Px;
					InputPxData[PxOffs] = Px;
				}
			}
			if(nPxWritten < nPxPerTile) ThisTilePxData[nPxWritten].f32[0] = INFINITY;
		}
	}

	//! Cluster tiles together by palette
	if(Plan->nTilePalettes > 1) {
		//! Set up the cluster pointers
		uint32_t n, k;
		float *NextPtr = (float*)(TileClusters + Plan->nTilePalettes);
		for(k=0;k<Plan->nTilePalettes;k++) {
			TileClusters[k].Centroid = NextPtr, NextPtr += TILE_CLUSTER_DIMENSIONS;
			TileClusters[k].Training = NextPtr, NextPtr += TILE_CLUSTER_DIMENSIONS;
		}

		//! Renormalize tile values
		float ValueMins[TILE_CLUSTER_DIMENSIONS];
		float ValueMaxs[TILE_CLUSTER_DIMENSIONS];
		for(k=0;k<TILE_CLUSTER_DIMENSIONS;k++) ValueMins[k] = +INFINITY;
		for(k=0;k<TILE_CLUSTER_DIMENSIONS;k++) ValueMaxs[k] = -INFINITY;
		for(n=0;n<nNonBlankTiles;n++) for(k=0;k<TILE_CLUSTER_DIMENSIONS;k++) {
			float v = TileColourValues[n*TILE_CLUSTER_DIMENSIONS+k];
			if(v < ValueMins[k]) ValueMins[k] = v;
			if(v > ValueMaxs[k]) ValueMaxs[k] = v;
		}
		for(k=0;k<TILE_CLUSTER_DIMENSIONS;k++) {
			//! I have no idea why, but adding a large bias of 1.0
			//! improves results much better than adding a tiny bias.
			float Range = ValueMaxs[k] - ValueMins[k];
			float Scale = 1.0f / (1.0f + Range);
			for(n=0;n<nNonBlankTiles;n++) {
				//! For clustering purposes, we don't need to center before scaling
				float *vp = &TileColourValues[n*TILE_CLUSTER_DIMENSIONS+k];
				//*vp = (*vp - Mean) * Scale + Mean;
				*vp *= Scale;
			}
		}

		//! Perform actual clustering now
		uint32_t TileClusterPasses = Plan->nTileClusterPasses;
		if(!TileClusterPasses) TileClusterPasses = DEFAULT_TILECLUSTER_PASSES / Plan->nTilePalettes;
		Clusterize_Process(
			TileClusters,
			TileColourValues,
			TILE_CLUSTER_DIMENSIONS,
			Plan->nTilePalettes,
			nNonBlankTiles,
			ClusterListIndices,
			TileClusterPasses,
			0.0f,
			TileClusterWeights
		);
		Clusterize_GetClusterIndices_u8(
			TilePaletteIndices,
			TileClusters,
			Plan->nTilePalettes,
			ClusterListIndices
		);

		//! Expand palette indices so that empty tiles get mapped as well,
		//! since we removed them during clustering to avoid empty palettes.
		//! Note that this loop works backwards, since we have to expand
		//! the data out, and going forwards would overwrite the data.
		if(BlankTileIndices < BlankTileIndicesEnd) {
			uint32_t NextTileIdx = nTilesTotal-1;
			      uint8_t *TileIndexDst = TilePaletteIndices + nTilesTotal;
			const uint8_t *TileIndexSrc = TilePaletteIndices + nNonBlankTiles;
			do {
				uint32_t NextBlankTileIdx = *BlankTileIndices++;

				//! Write non-empty tile indices
				while(NextTileIdx > NextBlankTileIdx) {
					*--TileIndexDst = *--TileIndexSrc;
					NextTileIdx--;
				}

				//! And now set the index of the empty tile
				*--TileIndexDst = 0;
				NextTileIdx--;
			} while(BlankTileIndices < BlankTileIndicesEnd);
		}
	} else {
		//! If we only requested one palette, then all tiles belong
		//! to the same palette index
		uint32_t n;
		for(n=0;n<nTilesTotal;n++) TilePaletteIndices[n] = 0;
	}

	//! And now cluster the palette colours
	uint32_t PalIdx;
	uint32_t nOutputColours = Plan->nPaletteColours;
	float SplitRatio = Plan->SplitRatio;
	if(SplitRatio < 0.0f) {
		if(nOutputColours <= 16) SplitRatio = 0.0f;
		else SplitRatio = 1.0f - powf(2.0f, 1.0f - (float)nOutputColours/16.0f);
	}
	if(Plan->FirstColourIsTransparent) nOutputColours--;
	uint32_t ColourClusterPasses = Plan->nColourClusterPasses;
	if(!ColourClusterPasses) ColourClusterPasses = DEFAULT_COLOURCLUSTER_PASSES / nOutputColours;
	for(PalIdx=0;PalIdx<Plan->nTilePalettes;PalIdx++) {
		//! Read pixels of all tiles falling into this palette
		uint32_t n;
		uint32_t DataCnt = 0;
		const uint8_t *PalIndices = TilePaletteIndices;
		for(n=0;n<nTilesTotal;n++) {
			if(PalIndices[n] == PalIdx) {
				uint32_t k;
				const Vec4f_t *Src = TilePxData + n*nPxPerTile;
				for(k=0;k<nPxPerTile;k++) {
					Vec4f_t Px = Src[k];

					//! Stop after hitting the last pixel stored
					if(Px.f32[0] == INFINITY) break;

					//! If we have a transparent palette colour, then
					//! don't add pixels with alpha=0.
					if(!Plan->FirstColourIsTransparent || Px.f32[3] != 0.0f) {
						ColourClusterBuffer[DataCnt++] = Px;
					}
				}
			}
		}

		//! Perform cluster analysis
		//! If a palette is unused, mark with magenta, since this
		//! really should never happen.
		Vec4f_t UnusedColour = (Vec4f_t){{1.0f,0.0f,1.0f,1.0f}};
		UnusedColour = ConvertToColourspace(&UnusedColour, Plan->Colourspace);
		if(DataCnt != 0) {
			uint32_t nOutputClusters = Clusterize_Vec4f_Process(
				ColourClusters,
				ColourClusterBuffer,
				nOutputColours,
				DataCnt,
				ClusterListIndices,
				ColourClusterPasses,
				SplitRatio,
				NULL
			);
			for(n=nOutputClusters;n<nOutputColours;n++) {
				ColourClusters[n].Centroid = UnusedColour;
			}
		} else {
			for(n=0;n<nOutputColours;n++) {
				ColourClusters[n].Centroid = UnusedColour;
			}
		}

		//! Extract palette colours from the centroids
		Vec4f_t *ThisPalette = PaletteData + PalIdx*Plan->nPaletteColours;
		if(Plan->FirstColourIsTransparent) *ThisPalette++ = VEC4F_EMPTY;
		for(n=0;n<nOutputColours;n++) ThisPalette[n] = ColourClusters[n].Centroid;

		//! Now sort the palette colours
		//! This step is not needed, but it can help when
		//! running palette effects on the image
		if(DataCnt != 0) {
			uint32_t k;

			//! Generate hue values for each palette entry
			//! Note that we quantize the hue to a standard hexagon,
			//! plus corners, giving us 12 sides (-6 .. +6). We then
			//! "unpack" the hue theta as cos/sin for clustering, as
			//! this avoids having to deal with wraparound distance.
			for(n=0;n<nOutputColours;n++) {
				const float Scale = (float)(6.0/M_PI);
				Vec4f_t x = ConvertFromColourspace(&ThisPalette[n], Plan->Colourspace);
				x = ConvertToColourspace(&x, COLOURSPACE_OKLAB);
				float h = atan2f(x.f32[2], x.f32[1]);
				h = roundf(h * Scale) / Scale;
				ColourClusterBuffer[n] = VEC4F_EMPTY;
				ColourClusterBuffer[n].f32[0] = cosf(h);
				ColourClusterBuffer[n].f32[1] = sinf(h);
			}

			//! Now separate the palette colours into clusters
			uint32_t nPalClusters = Clusterize_Vec4f_Process(
				ColourClusters,
				ColourClusterBuffer,
				nOutputColours,
				nOutputColours,
				ClusterListIndices,
				PALETTE_SORT_PASSES,
				0.0f,
				NULL
			);

			//! Sort the clusters by hue
			//! NOTE: We abuse the ::NextCluster member to store the sort order
			for(k=0;k<nPalClusters;k++) ColourClusters[k].NextCluster = k;
			qsort(
				ColourClusters,
				nPalClusters,
				sizeof(struct Cluster_Vec4f_t),
				CompareHueAscending
			);

			//! Now store the colours of each cluster sequentially,
			//! sorting each cluster's colours by their luma value
			uint8_t  IsAscending = 1;
			uint32_t nSingleColours = 0;
			Vec4f_t *NextColour       = ThisPalette;
			Vec4f_t *NextSingleColour = ThisPalette + nOutputColours;
			for(n=0;n<nOutputColours;n++) {
				ColourClusterBuffer[n] = ConvertFromColourspace(&ThisPalette[n], Plan->Colourspace);
			}
			for(k=0;k<nPalClusters;k++) {
				uint32_t Next = ColourClusters[ColourClusters[k].NextCluster].FirstDataIdx;
				if(Next == CLUSTER_END_OF_LIST) continue; //! Empty cluster
				if(ClusterListIndices[Next] != CLUSTER_END_OF_LIST) { //! Have more than one colour?
					uint32_t nClusterColours = 0;
					do {
						NextColour[nClusterColours++] = ColourClusterBuffer[Next];
						Next = ClusterListIndices[Next];
					} while(Next != CLUSTER_END_OF_LIST);
					qsort(
						NextColour,
						nClusterColours,
						sizeof(Vec4f_t),
						IsAscending ? CompareLumaAscending : CompareLumaDescending
					);
					NextColour += nClusterColours;
					IsAscending ^= 1;
				} else {
					//! Put single colours at the back
					*--NextSingleColour = ColourClusterBuffer[Next];
					nSingleColours++;
				}
			}
			if(nSingleColours > 0) qsort(
				NextSingleColour,
				nSingleColours,
				sizeof(Vec4f_t),
				IsAscending ? CompareLumaAscending : CompareLumaDescending
			);
			for(n=0;n<nOutputColours;n++) {
				ThisPalette[n] = ConvertToColourspace(&ThisPalette[n], Plan->Colourspace);
			}
		}
	}

	//! Finally, convert the palette from perceptual space to RGB space
	for(PalIdx=0;PalIdx<Plan->nTilePalettes;PalIdx++) {
		uint32_t n, k;
		BGRA8_t *Dst = OutputPalette + PalIdx*Plan->nPaletteColours;
		Vec4f_t *Src = PaletteData   + PalIdx*Plan->nPaletteColours;
		Vec4f_t vZeros = Vec4f_Broadcast(0.0f);
		Vec4f_t vOnes  = Vec4f_Broadcast(1.0f);
		for(n=0;n<Plan->nPaletteColours;n++) {
			//! Convert back to RGB, undo pre-multiplied alpha as needed, and quantize
			Vec4f_t x = ConvertFromColourspace(&Src[n], Plan->Colourspace);
			if(!Plan->PremultipliedAlpha) {
				float a = x.f32[3];
				if(a != 0.0f) {
					x.f32[0] /= a;
					x.f32[1] /= a;
					x.f32[2] /= a;
				}
			}
			x = Vec4f_Max(&x, &vZeros);
			x = Vec4f_Min(&x, &vOnes);
			Vec4f_t xq = Vec4f_Quantize(&x, &Plan->ColourDepth);
			for(k=0;k<n;k++) {
				//! If we already quantized to this exact colour,
				//! split to floor/ceiling modes to get a better
				//! chance at error dithering
				Vec4f_t xk = BGRA8_To_Vec4fRGBA(&Dst[k]);
				if(Vec4f_Dist2(&xq, &xk) == 0.0f) {
					xq = Vec4f_QuantizeFloor(&x, &Plan->ColourDepth);
					Dst[k] = Vec4fRGBA_To_BGRA8(&xq);
					xq = Vec4f_QuantizeCeil(&x, &Plan->ColourDepth);
				}
			}
			Dst[n] = Vec4fRGBA_To_BGRA8(&xq);

			//! Pre-multiply by alpha again and convert to final colourspace
			if(!Plan->PremultipliedAlpha) xq = Vec4f_Muli(&xq, xq.f32[3]);
			Src[n] = ConvertToColourspace(&xq, Plan->Colourspace);
		}
	}

	//! Apply post-dithering and final colour matching
	DitherPaletteImage(
		OutputPxData,
		InputPxData,
		PaletteData,
		TilePaletteIndices,
		InputWidth,
		InputHeight,
		Plan->DitherType,
		Plan->DitherLevel,
		Plan->TileWidth,
		Plan->TileHeight,
		Plan->nPaletteColours,
		RMSE
	);

	//! Clean up
	free(AllocBuffer);
	return 1;
}

/************************************************/
//! EOF
/************************************************/

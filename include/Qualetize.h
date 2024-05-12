/************************************************/
#pragma once
/************************************************/
#include <stdint.h>
/************************************************/
#include "Bitmap.h"
#include "Vec4f.h"
/************************************************/

//! Dithering modes available
//! NOTE: Ordered dithering gives consistent tiled results, but Floyd-Steinberg can look nicer.
#define DITHER_NONE           (   0) //! No dithering
#define DITHER_ORDERED(n)     (   n) //! Ordered dithering (Kernel size: (2^n) x (2^n))
#define DITHER_CHECKER        (0xFF) //! Checkerboard pattern
#define DITHER_FLOYDSTEINBERG (0xFE) //! Floyd-Steinberg (diffusion)
#define DITHER_ATKINSON       (0xFD) //! Atkinson (diffusion)

/************************************************/

//! Colourspaces available
#define COLOURSPACE_SRGB       0 //! sRGB
#define COLOURSPACE_RGB_LINEAR 1 //! RGB (linear light)
#define COLOURSPACE_YCBCR      2 //! YCbCr (using ITU-R BT.709 constants)
#define COLOURSPACE_YCOCG      3 //! YCoCg
#define COLOURSPACE_CIELAB     4 //! CIELAB (D65 lighting)
#define COLOURSPACE_ICTCP      5 //! ICtCp
#define COLOURSPACE_OKLAB      6 //! OkLab
#define COLOURSPACE_RGB_PSY    7 //! RGB + Psyopt
#define COLOURSPACE_YCBCR_PSY  8 //! YCbCr + Psyopt
#define COLOURSPACE_YCOCG_PSY  9 //! YCoCg + Psyopt

/************************************************/

//! Plan for qualetize operation
struct QualetizePlan_t {
	uint16_t TileWidth;                //! Tile dimensions (in pixels)
	uint16_t TileHeight;
	uint16_t nPaletteColours;          //! Colours per palette
	uint16_t nTilePalettes;            //! Number of palettes to generate
	uint8_t  FirstColourIsTransparent; //! 0 = All palette entries used,  1 = First palette entry is always (0,0,0,0)
	uint8_t  PremultipliedAlpha;       //! 0 = Multiply colours by alpha, 1 = Colours are pre-multiplied by alpha
	uint8_t  DitherInputType;          //! Dithering to apply on conversion to target RGB range
	uint8_t  DitherOutputType;         //! Dithering to apply on final output
	uint16_t DitherInputLevel;         //! Dithering level to use on conversion to target RGB range (0 = None, 8000h = 100%)
	uint16_t DitherOutputLevel;        //! Dithering level to use on final output (0 = None, 8000h = 100%)
	uint8_t  Colourspace;              //! Colourspace to use during processing
	uint8_t  r1[3];
	uint32_t nTileClusterPasses;       //! Number of clustering passes to apply on tiles
	uint32_t nColourClusterPasses;     //! Number of clustering passes to apply on tile colours
	Vec4f_t  ColourDepth;              //! RGBA levels for output image
};

/************************************************/

//! Apply qualetize operation to input image using the specified plan
//! Returns 0 on failure, or 1 on sucess.
//! Notes:
//!  -If InputPalette is NULL, then InputBitmap is assumed to be BGRA8_t.
//!   Otherwise (InputPalette != NULL), InputBitmap is assumed to be
//!   uint8_t, using InputPalette for its colours.
//!  -RMSE is optional and can be NULL.
uint8_t Qualetize(
	      uint8_t  *OutputPxData,
	      BGRA8_t  *OutputPalette,
	const void     *InputBitmap,
	const BGRA8_t  *InputPalette,
	      uint32_t  InputWidth,
	      uint32_t  InputHeight,
	const struct QualetizePlan_t *Plan,
	Vec4f_t *RMSE
);

/************************************************/

//! Dither image data, given the target bit depth
//! Notes:
//!  -If using Floyd-Steinberg dithering, then if the function
//!   can't allocate the memory for it, it will disable dithering.
void DitherImage(
	      Vec4f_t *DstPx,
	const Vec4f_t *SrcPx,
	uint32_t Width,
	uint32_t Height,
	uint8_t  DitherType,
	float    DitherLevel,
	const Vec4f_t *ColourDepth
);

//! Dither palettized, tiled image data
//! Notes:
//!  -If using dithering, then if the function can't allocate the
//!   memory for the dither buffers, it will disable dithering.
//!  -RMSE is optional, and can be set to NULL.
void DitherPaletteImage(
	      uint8_t *DstPx,
	const Vec4f_t *SrcPx,
	const Vec4f_t *Palette,
	const uint8_t *TilePalettes,
	uint32_t Width,
	uint32_t Height,
	uint8_t  DitherType,
	float    DitherLevel,
	uint32_t TileWidth,
	uint32_t TileHeight,
	uint32_t nPaletteColours,
	uint32_t nTilePalettes,
	Vec4f_t *RMSE
);

/************************************************/
//! EOF
/************************************************/


/************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/************************************************/
#include "Bitmap.h"
#include "Qualetize.h"
/************************************************/

//! When not zero, the PSNR for each channel will be displayed
#define MEASURE_PSNR 1

/************************************************/

//! strcmp() implementation that ACTUALLY returns the difference between
//! characters instead of just the signs. Blame the C standard -_-
static int mystrcmp(const char *s1, const char *s2, const char **s1End) {
	while(*s1 != '\0' && *s1 == *s2) s1++, s2++;
	if(s1End) *s1End = s1;
	return (int)(*s1) - (int)(*s2);
}

/************************************************/

//! Convert symbolic colourspace name to pretty string
static const char *ColourspaceNameString(uint8_t Colourspace) {
	switch(Colourspace) {
		case COLOURSPACE_SRGB:      return "sRGB";
		case COLOURSPACE_YCBCR:     return "YCbCr";
		case COLOURSPACE_YCOCG:     return "YCoCg";
		case COLOURSPACE_CIELAB:    return "CIELAB";
		case COLOURSPACE_ICTCP:     return "ICtCp";
		case COLOURSPACE_OKLAB:     return "OkLab";
		case COLOURSPACE_RGB_PSY:   return "RGB + Psyopt";
		case COLOURSPACE_YCBCR_PSY: return "YCbCr + Psyopt";
		case COLOURSPACE_YCOCG_PSY: return "YCoCg + Psyopt";
		default: return "[Unknown colourspace]";
	}
}

/************************************************/

static int ParseColourspace(const char *s) {
	     if(!strcmp(s, "srgb"))      return COLOURSPACE_SRGB;
	else if(!strcmp(s, "ycbcr"))     return COLOURSPACE_YCBCR;
	else if(!strcmp(s, "ycocg"))     return COLOURSPACE_YCOCG;
	else if(!strcmp(s, "cielab"))    return COLOURSPACE_CIELAB;
	else if(!strcmp(s, "ictcp"))     return COLOURSPACE_ICTCP;
	else if(!strcmp(s, "oklab"))     return COLOURSPACE_OKLAB;
	else if(!strcmp(s, "rgb-psy"))   return COLOURSPACE_RGB_PSY;
	else if(!strcmp(s, "ycbcr-psy")) return COLOURSPACE_YCBCR_PSY;
	else if(!strcmp(s, "ycocg-psy")) return COLOURSPACE_YCOCG_PSY;
	else return -1;
}

static int ParseDitherMode(const char *s, uint8_t *ModePtr, uint16_t *LevelPtr) {
	int   r; //! <- Will contain output from mystrcmp()
	int   Mode;
	float Level;
	const char *MatchEnd;
#define DITHERMODE_MATCH(Name) \
	(((r = mystrcmp(s, Name, &MatchEnd)) || 1) && (r == '\0' || r == ','))
	     if(DITHERMODE_MATCH("none" ))    Mode = DITHER_NONE,           Level = 0.0f;
	else if(DITHERMODE_MATCH("floyd"))    Mode = DITHER_FLOYDSTEINBERG, Level = 1.0f;
	else if(DITHERMODE_MATCH("atkinson")) Mode = DITHER_ATKINSON,       Level = 1.0f;
	else if(DITHERMODE_MATCH("checker"))  Mode = DITHER_CHECKER,        Level = 1.0f;
	else if(DITHERMODE_MATCH("ord2" ))    Mode = DITHER_ORDERED(1),     Level = 1.0f;
	else if(DITHERMODE_MATCH("ord4" ))    Mode = DITHER_ORDERED(2),     Level = 1.0f;
	else if(DITHERMODE_MATCH("ord8" ))    Mode = DITHER_ORDERED(3),     Level = 1.0f;
	else if(DITHERMODE_MATCH("ord16"))    Mode = DITHER_ORDERED(4),     Level = 1.0f;
	else if(DITHERMODE_MATCH("ord32"))    Mode = DITHER_ORDERED(5),     Level = 1.0f;
	else if(DITHERMODE_MATCH("ord64"))    Mode = DITHER_ORDERED(6),     Level = 1.0f;
	else return -1;
#undef DITHERMODE_MATCH
	if(Mode != DITHER_NONE && r == ',') {
		Level = atof(MatchEnd+1); //! <- Skip over the ',' symbol
	}
	Level *= 32768.0f;
	if(Level <     0.0f) Level =     0.0f;
	if(Level > 65535.0f) Level = 65535.0f;
	*ModePtr  = (uint8_t )Mode;
	*LevelPtr = (uint16_t)Level;
	return Mode;
}

int main(int argc, const char *argv[]) {
	//! Check arguments
	if(argc < 3) {
		printf(
			"qualetize - Tiled colour-quantization tool\n"
			"Usage:\n"
			" qualetize Input.bmp Output.bmp [options]\n"
			"Options:\n"
			"  -tw:8                - Set tile width\n"
			"  -th:8                - Set tile height\n"
			"  -npal:16             - Set number of palettes available\n"
			"  -cols:16             - Set number of colours per palette\n"
			"                         Note that this value times the number of palettes must\n"
			"                         be less than or equal to 256.\n"
			"  -rgba:5551           - Set RGBA bit depth\n"
			"                         RGBA = 8888 is standard for BMP (ie. 24-bit colour,\n"
			"                         plus 8-bit alpha), but for the targets this tool is\n"
			"                         intended for, RGBA = 5551 is the norm.\n"
			"  -premulalpha:n       - Alpha is pre-multiplied (y/n)\n"
			"                         While most formats generally pre-multiply the colours\n"
			"                         by the alpha value, 32-bit BMP files generally do not.\n"
			"                         Note that if this option is set to `y`, then output\n"
			"                         colours in the palette will also be pre-multiplied.\n"
			"  -colspace:oklab      - Set colourspace\n"
			"                         Different colourspaces may give better/worse results\n"
			"                         depending on the input image, and it may be necessary\n"
			"                         to experiment to find the optimal one.\n"
			"  -ditherin:none       - Set dither mode, level for input\n"
			"                         This option will dither the input image (in sRGB space\n"
			"                         to get the \"correct\" output colours), and use this as\n"
			"                         the input for generating the palette data. This option\n"
			"                         may generate unused palette entries corresponding to\n"
			"                         dithered colours that weren't used during output, but\n"
			"                         can help remove tile blocking/banding artifacts.\n"
			"  -ditherout:floyd,1.0 - Set dither mode, level for output\n"
			"                         This can reduce some of the banding artifacts caused\n"
			"                         when the colours per palette is very small, at the\n"
			"                         expense of added \"noise\".\n"
			"  -tilepasses:0        - Set tile cluster passes (0 = default)\n"
			"  -colourpasses:0      - Set colour cluster passes (0 = default)\n"
			"                         Most of the processing time will be spent in the loop\n"
			"                         that clusters the colours together. If processing is\n"
			"                         taking excessive amounts of time, this option may be\n"
			"                         adjusted (eg. for 256-colour palettes, it is suggested\n"
			"                         to set this value to something like 4, whereas for\n"
			"                         16-colour palettes, a value between 32 and 64 should\n"
			"                         result in convergence in most cases, while still being\n"
			"                         of reasonable performance).\n"
			"  -col0isclear:y       - First colour of every palette is transparent (y/n)\n"
			"                         Note that this affects both input AND output images.\n"
			"                         To set transparency in a direct-colour input bitmap,\n"
			"                         an alpha channel must be used (ie. 32-bit input);\n"
			"                         translucent alpha values are supported by this tool.\n"
			"Colourspaces available:\n"
			"  srgb\n"
			"  rgb-psy      (Psy = Non-linear light, weighted components)\n"
			"  ycbcr[-psy]  (Psy = Non-linear luma, weighted chroma)\n"
			"  ycocg[-psy]  (Psy = Non-linear luma)\n"
			"  cielab\n"
			"    NOTE: CIELAB has poor performance in most cases.\n"
			"  ictcp\n"
			"  oklab\n"
			"Dither modes available (and default level):\n"
			"  none         - No dithering\n"
			"  floyd,1.0    - Floyd-Steinberg\n"
			"  atkinson,1.0 - Atkinson diffusion\n"
			"  checker,1.0  - Chekerboard dithering\n"
			"  ord2,1.0     - 2x2 ordered dithering\n"
			"  ord4,1.0     - 4x4 ordered dithering\n"
			"  ord8,1.0     - 8x8 ordered dithering\n"
			"  ord16,1.0    - 16x16 ordered dithering\n"
			"  ord32,1.0    - 32x32 ordered dithering\n"
			"  ord64,1.0    - 64x64 ordered dithering\n"
			"Note that ordered dithering on the output uses an approximation to the gradient\n"
			"slope only, as the palette is very unlikely to have a uniform spread.\n"
		);
		return 1;
	}

	//! Parse arguments
	struct QualetizePlan_t Plan;
	Plan.TileWidth            = 8;
	Plan.TileHeight           = 8;
	Plan.nPaletteColours      = 16;
	Plan.nTilePalettes        = 16;
	Plan.FirstColourIsTransparent = 1;
	Plan.PremultipliedAlpha   = 0;
	Plan.DitherInputType      = DITHER_NONE;
	Plan.DitherInputLevel     = 0;
	Plan.DitherOutputType     = DITHER_FLOYDSTEINBERG;
	Plan.DitherOutputLevel    = 0x8000;
	Plan.Colourspace          = COLOURSPACE_OKLAB;
	Plan.nTileClusterPasses   = 0;
	Plan.nColourClusterPasses = 0;
	Plan.ColourDepth          = (Vec4f_t){{31,31,31,1}};
	{
		int argi;
		for(argi=3;argi<argc;argi++) {
			uint8_t ArgOk = 0;

			const char *ArgStr;
#define ARGMATCH(Input, Target) \
	ArgStr = Input + strlen(Target); \
	if(!memcmp(Input, Target, strlen(Target)))
			ARGMATCH(argv[argi], "-tw:")   ArgOk = 1, Plan.TileWidth       = atoi(ArgStr);
			ARGMATCH(argv[argi], "-th:")   ArgOk = 1, Plan.TileHeight      = atoi(ArgStr);
			ARGMATCH(argv[argi], "-npal:") ArgOk = 1, Plan.nTilePalettes   = atoi(ArgStr);
			ARGMATCH(argv[argi], "-cols:") ArgOk = 1, Plan.nPaletteColours = atoi(ArgStr);
			ARGMATCH(argv[argi], "-rgba:") {
#define READCHANNEL(Target) (((Target = (int)*ArgStr++ - '0') || 1) && (Target > 0 && Target <= 8))
				int r, g, b, a;
				const char *Depth = ArgStr;
				if(
					READCHANNEL(r) &&
					READCHANNEL(g) &&
					READCHANNEL(b) &&
					READCHANNEL(a)
				) {
					Plan.ColourDepth.f32[0] = (float)((1<<r) - 1);
					Plan.ColourDepth.f32[1] = (float)((1<<g) - 1);
					Plan.ColourDepth.f32[2] = (float)((1<<b) - 1);
					Plan.ColourDepth.f32[3] = (float)((1<<a) - 1);
				} else printf("WARNING: Unrecognized RGBA depth sequence: %s\n", Depth);
				ArgOk = 1;
#undef READCHANNEL
			}
			ARGMATCH(argv[argi], "-premulalpha:") ArgOk = 1, Plan.PremultipliedAlpha = (ArgStr[0] == 'y') ? 1 : 0;
			ARGMATCH(argv[argi], "-colspace:") {
				int c = ParseColourspace(ArgStr);
				if(c != -1) Plan.Colourspace = (uint8_t)c;
				else printf("WARNING: Unrecognized colourspace: %s\n", ArgStr);
				ArgOk = 1;
			}
			ARGMATCH(argv[argi], "-ditherin:") {
				if(ParseDitherMode(ArgStr, &Plan.DitherInputType, &Plan.DitherInputLevel) == -1) {
					printf("WARNING: Unrecognized input dither mode: %s\n", ArgStr);
				}
				ArgOk = 1;
			}
			ARGMATCH(argv[argi], "-ditherout:") {
				if(ParseDitherMode(ArgStr, &Plan.DitherOutputType, &Plan.DitherOutputLevel) == -1) {
					printf("WARNING: Unrecognized output dither mode: %s\n", ArgStr);
				}
				ArgOk = 1;
			}
			ARGMATCH(argv[argi], "-tilepasses:")   ArgOk = 1, Plan.nTileClusterPasses       = atoi(ArgStr);
			ARGMATCH(argv[argi], "-colourpasses:") ArgOk = 1, Plan.nColourClusterPasses     = atoi(ArgStr);
			ARGMATCH(argv[argi], "-col0isclear:")  ArgOk = 1, Plan.FirstColourIsTransparent = (ArgStr[0] == 'y') ? 1 : 0;
#undef ARGMATCH
			//! Unrecognized?
			if(!ArgOk) printf("WARNING: Unrecognized argument: %s\n", argv[argi]);
		}
	}

	//! Sanity checks
	if(Plan.nPaletteColours == 0) {
		printf("ERROR: Number of colours per palette must not be 0.\n");
		return -1;
	}
	if(Plan.nTilePalettes == 0) {
		printf("ERROR: Number of palettes must not be 0.\n");
		return -1;
	}
	if(Plan.nPaletteColours*Plan.nTilePalettes > 256) {
		printf("ERROR: Exceeded maximum of 256 colours.\n");
		return -1;
	}
	if(Plan.nTileClusterPasses > 1000) {
		printf("WARNING: Tile cluster passes is very large (%u); this may be slow.\n", Plan.nTileClusterPasses);
	}
	if(Plan.nColourClusterPasses > 100) {
		printf("WARNING: Colour cluster passes is very large (%u); this may be slow.\n", Plan.nColourClusterPasses);
	}

	//! Open input image
	struct BmpCtx_t Image;
	if(!BmpCtx_FromFile(&Image, argv[1])) {
		printf("ERROR: Unable to read input file.\n");
		return -1;
	}
	if((Image.Width%Plan.TileWidth) != 0 || (Image.Height%Plan.TileHeight) != 0) {
		printf(
			"WARNING: Image not a multiple of tile size (%dx%d); edges will be extended.\n",
			Plan.TileWidth,
			Plan.TileHeight
		);
	}

	//! Create output image
	struct BmpCtx_t Output;
	if(!BmpCtx_Create(&Output, Image.Width, Image.Height, 1)) {
		printf("ERROR: Couldn't create output image.\n");
		BmpCtx_Destroy(&Image);
		return -1;
	}

	//! Perform processing
	Vec4f_t RMSE;
	if(!Qualetize(
		Output.PxIdx,
		Output.Palette,
		Image.ImgData,
		Image.Palette,
		Image.Width,
		Image.Height,
		&Plan,
		&RMSE)
	) {
		printf("ERROR: Out of memory.\n");
		BmpCtx_Destroy(&Output);
		BmpCtx_Destroy(&Image);
		return -1;
	}
	BmpCtx_Destroy(&Image);

	//! Save image to file
	if(!BmpCtx_ToFile(&Output, argv[2])) {
		printf("ERROR: Unable to write output file.\n");
		BmpCtx_Destroy(&Output);
		return -1;
	}

	//! Display PSNR
#if MEASURE_PSNR
	Vec4f_t vTotalSSE = Vec4f_Mul(&RMSE, &RMSE);
	float TotalRMSE = -8.68589f * 0.5f*logf(Vec4f_SumOf(&vTotalSSE));
	RMSE.f32[0]     = -8.68589f * logf(RMSE.f32[0]); //! -20*Log10[RMSE] == -20/Log[10] * Log[RMSE]
	RMSE.f32[1]     = -8.68589f * logf(RMSE.f32[1]);
	RMSE.f32[2]     = -8.68589f * logf(RMSE.f32[2]);
	RMSE.f32[3]     = -8.68589f * logf(RMSE.f32[3]);
	printf(
		"PSNR (as %s) = %.3fdB {%.3fdB, %.3fdB, %.3fdB, %.3fdB}\n"
		"Ok.\n",
		ColourspaceNameString(Plan.Colourspace),
		TotalRMSE,
		RMSE.f32[0],
		RMSE.f32[1],
		RMSE.f32[2],
		RMSE.f32[3]
	);
#endif

	//! Success
	BmpCtx_Destroy(&Output);
	return 0;
}

/************************************************/
//! EOF
/************************************************/

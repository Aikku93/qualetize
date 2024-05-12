/************************************************/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/************************************************/
#include "Bitmap.h"
/************************************************/

//! Clear context data
#define CLEAR_CONTEXT(Ctx)   \
	Ctx->Width   = 0,    \
	Ctx->Height  = 0,    \
	Ctx->Palette = NULL, \
	Ctx->PxBGR   = NULL

//! Destroy context and return 0
#define DESTROY_AND_RETURN(Ctx, ...) \
	do {                         \
		BmpCtx_Destroy(Ctx); \
		return __VA_ARGS__;  \
	} while(0)

/************************************************/

//! BMP file header
#pragma pack(push,1)
struct BMFH_t {
	uint16_t Type;
	uint32_t Size;
	uint16_t r1[2];
	uint32_t Offs;
};

//! BMP image header
struct BMIH_t {
	uint32_t Size;
	uint32_t Width;
	uint32_t Height;
	uint16_t nPlanes;
	uint16_t BitCnt;
	uint32_t CompType;
	uint32_t ImgSize;
	uint32_t XpM;
	uint32_t YpM;
	uint32_t ColUsed;
	uint32_t ColImportant;
};
#pragma pack(pop)

/************************************************/

//! Swap ranges of memory
//! This is very much sub-optimal
static void SwapRange(void *a_, void *b_, uint32_t Length) {
	uint8_t *a = (uint8_t*)a_;
	uint8_t *b = (uint8_t*)b_;
	while(Length--) {
		uint8_t t = *a;
		*a++ = *b;
		*b++ = t;
	}
}

/************************************************/

//! Create context
uint8_t BmpCtx_Create(struct BmpCtx_t *Ctx, uint32_t w, uint32_t h, uint8_t UsePal) {
	Ctx->Width  = w;
	Ctx->Height = h;
	if(UsePal) {
		Ctx->Palette = calloc(BMP_PALETTE_COLOURS, sizeof(BGRA8_t));
		Ctx->PxIdx   = calloc(w*h,  sizeof(uint8_t));
		if(!Ctx->Palette || !Ctx->PxIdx) DESTROY_AND_RETURN(Ctx, 0);
	} else {
		Ctx->Palette = NULL;
		Ctx->PxBGR   = calloc(w*h, sizeof(BGRA8_t));
		if(!Ctx->PxBGR) DESTROY_AND_RETURN(Ctx, 0);
	}
	return 1;
}

/************************************************/

//! Destroy context
void BmpCtx_Destroy(struct BmpCtx_t *Ctx) {
	free(Ctx->Palette);
	free(Ctx->PxBGR);
	CLEAR_CONTEXT(Ctx);
}

/************************************************/

//! Load from file
uint8_t BmpCtx_FromFile(struct BmpCtx_t *Ctx, const char *Filename) {
	uint8_t ExitCode = 0;
	struct BMFH_t bmFH;
	struct BMIH_t bmIH;
	CLEAR_CONTEXT(Ctx);

	//! Open file, read headers
	FILE *File = fopen(Filename, "rb");
	if(!File) return 0;
	if(!fread(&bmFH, sizeof(bmFH), 1, File)) goto Exit;
	if(!fread(&bmIH, sizeof(bmIH), 1, File)) goto Exit;
	Ctx->Width  = bmIH.Width;
	Ctx->Height = bmIH.Height;

	//! Read pixels
	if(bmFH.Type == ('B'|'M'<<8)) {
		uint32_t nPx = Ctx->Width * Ctx->Height;
		switch(bmIH.BitCnt) {
			//! 8-bit palettized
			case 8: {
				//! Read palette
				Ctx->Palette = malloc(BMP_PALETTE_COLOURS * sizeof(BGRA8_t));
				if(!Ctx->Palette) goto Exit;
				if(!fread(Ctx->Palette, BMP_PALETTE_COLOURS * sizeof(BGRA8_t), 1, File)) goto Exit;

				//! Read pixels
				fseek(File, bmFH.Offs, SEEK_SET);
				Ctx->PxIdx = malloc(nPx * sizeof(uint8_t));
				if(!Ctx->PxIdx) goto Exit;
				if(!fread(Ctx->PxIdx, nPx * sizeof(uint8_t), 1, File)) goto Exit;

				//! Unflip image
				uint32_t y;
				uint8_t *Mem = Ctx->PxIdx;
				for(y=0;y<bmIH.Height/2;y++) {
					SwapRange(Mem + y*bmIH.Width, Mem + (bmIH.Height-1-y)*bmIH.Width, bmIH.Width * sizeof(uint8_t));
				}
			} break;

			//! BGR
			case 24: {
				BGRA8_t *PxBGR = Ctx->PxBGR = malloc(nPx * sizeof(BGRA8_t));
				if(!PxBGR) goto Exit;

				//! Convert pixels
				uint32_t n;
				for(n=0;n<nPx;n++) {
					//! Read BGR
					struct { uint8_t b, g, r; } x;
					if(!fread(&x, 3, 1, File)) goto Exit;

					//! Store BGRA
					PxBGR[n].b = x.b;
					PxBGR[n].g = x.g;
					PxBGR[n].r = x.r;
					PxBGR[n].a = 255;
				}

				//! Unflip image
				uint32_t y;
				BGRA8_t *Mem = Ctx->PxBGR;
				for(y=0;y<bmIH.Height/2;y++) {
					SwapRange(Mem + y*bmIH.Width, Mem + (bmIH.Height-1-y)*bmIH.Width, bmIH.Width * sizeof(BGRA8_t));
				}
			} break;

			//! BGRA
			case 32: {
				//! Everything is prepared already, so straight read
				Ctx->PxBGR = malloc(nPx * sizeof(BGRA8_t));
				if(!Ctx->PxBGR) goto Exit;
				if(!fread(Ctx->PxBGR, nPx * sizeof(BGRA8_t), 1, File)) goto Exit;

				//! Unflip image
				uint32_t y;
				BGRA8_t *Mem = Ctx->PxBGR;
				for(y=0;y<bmIH.Height/2;y++) {
					SwapRange(Mem + y*bmIH.Width, Mem + (bmIH.Height-1-y)*bmIH.Width, bmIH.Width * sizeof(BGRA8_t));
				}
			} break;

			//! Unknown
			default: {
				goto Exit;
			} break;
		}
	}

	//! If we got here, we'll all good
	ExitCode = 1;
Exit:
	//! Close file, check success
	fclose(File);
	if(!ExitCode) DESTROY_AND_RETURN(Ctx, 0);
	return ExitCode;
}

/************************************************/

//! Write to file
uint8_t BmpCtx_ToFile(const struct BmpCtx_t *Ctx, const char *Filename) {
	uint8_t ExitCode = 0;
	struct BMFH_t bmFH;
	struct BMIH_t bmIH;

	//! Check image is valid
	uint32_t nPx = Ctx->Width * Ctx->Height;
	if(!nPx || (!Ctx->PxBGR && !(Ctx->Palette && Ctx->PxIdx))) return 0;

	//! Open file, write headers
	FILE *File = fopen(Filename, "wb");
	if(!File) return 0;
	memset(&bmFH, 0, sizeof(bmFH));
	memset(&bmIH, 0, sizeof(bmIH));
	bmFH.Type     = 'B'|'M'<<8;
	bmFH.Size     = sizeof(struct BMFH_t) +
		        sizeof(struct BMIH_t) +
		        BMP_PALETTE_COLOURS * (Ctx->Palette ? sizeof(BGRA8_t) : 0) +
		        nPx * (Ctx->Palette ? sizeof(uint8_t) : sizeof(BGRA8_t));
	bmFH.Offs     = sizeof(struct BMFH_t) + sizeof(struct BMIH_t) +
		        BMP_PALETTE_COLOURS * (Ctx->Palette ? sizeof(BGRA8_t) : 0);
	bmIH.Size     = sizeof(struct BMIH_t);
	bmIH.Width    = Ctx->Width;
	bmIH.Height   = Ctx->Height;
	bmIH.nPlanes  = 1;
	bmIH.BitCnt   = Ctx->Palette ? 8 : 32;
	fwrite(&bmFH, 1, sizeof(bmFH), File);
	fwrite(&bmIH, 1, sizeof(bmIH), File);

	//! Write palette
	if(Ctx->Palette) if(!fwrite(Ctx->Palette, BMP_PALETTE_COLOURS * sizeof(BGRA8_t), 1, File)) goto Exit;

	//! Write pixels (and flip image for storage)
	if(Ctx->Palette) {
		uint32_t y;
		const uint8_t *Mem = Ctx->PxIdx;
		for(y=0;y<Ctx->Height;y++) {
			if(!fwrite(Mem + (Ctx->Height-1-y)*Ctx->Width, Ctx->Width * sizeof(uint8_t), 1, File)) goto Exit;
		}
	} else {
		uint32_t y;
		const BGRA8_t *Mem = Ctx->PxBGR;
		for(y=0;y<Ctx->Height;y++) {
			if(!fwrite(Mem + (Ctx->Height-1-y)*Ctx->Width, Ctx->Width * sizeof(BGRA8_t), 1, File)) goto Exit;
		}
	}

	//! If we got here, we'll all good
	ExitCode = 1;
Exit:
	//! Close file
	fclose(File);
	return ExitCode;
}

/************************************************/
//! EOF
/************************************************/

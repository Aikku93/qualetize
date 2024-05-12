/************************************************/
#pragma once
/************************************************/
#include <stdint.h>
/************************************************/
#define BMP_PALETTE_COLOURS 256
/************************************************/

typedef struct {
	uint8_t b, g, r, a;
} BGRA8_t;

struct BmpCtx_t {
	uint32_t Width, Height;
	BGRA8_t *Palette;
	union {
		uint8_t *PxIdx; //! Palettized
		BGRA8_t *PxBGR; //! Direct
		void    *ImgData;
	};
};

/************************************************/

//! Create context
//! Pass UsePal=0 for BGRA, or UsePal=1 to allocate a palette.
//! Returns 0 on failure, or 1 on success.
uint8_t BmpCtx_Create(struct BmpCtx_t *Ctx, uint32_t w, uint32_t h, uint8_t UsePal);

//! Destroy context
void BmpCtx_Destroy(struct BmpCtx_t *Ctx);

//! Load from file
//! Returns 0 on failure, or 1 on success.
//! NOTE: 24bit BGR is converted to 32bit BGRA internally.
uint8_t BmpCtx_FromFile(struct BmpCtx_t *Ctx, const char *Filename);

//! Write to file
//! Returns 0 on failure, or 1 on success.
//! NOTE: To write a BGRA image, set Palette=nullptr.
//! NOTE: Always 32bit BGRA; 24bit BGR is never used for output.
uint8_t BmpCtx_ToFile(const struct BmpCtx_t *Ctx, const char *Filename);

/************************************************/
//! EOF
/************************************************/

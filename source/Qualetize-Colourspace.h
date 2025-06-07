/************************************************/
#pragma once
/************************************************/
#include <math.h>
#include <stdint.h>
/************************************************/
#include "Qualetize.h"
#include "Vec4f.h"
/************************************************/

//! Convert sRGB <-> Linear RGB
static float RGBtoLinearRGB(float t) {
	if(t > 0.04045f) return powf((t + 0.055f) / 1.055f, 2.4f);
	else return t / 12.92f;
}
static float LinearRGBtoRGB(float t) {
	if(t > 0.0031308f) return 1.055f*powf(t, 1.0f/2.4f) - 0.055f;
	else return 12.92f * t;
}

//! Convert sRGB <-> Visually-mapped RGB
//! This is roughly equivalent to Lab's transformation,
//! which approximately applies an exponent of 2.2 to
//! un-compand sRGB to linear RGB, followed by a cubic
//! root response, yielding an exponent of 2.2/3 = 0.73.
static float RGBtoVisualRGB(float t) {
	return (t > 0.0f) ? powf(t, (float)(2.2 / 3.0)) : 0.0f;
}
static float VisualRGBtoRGB(float t) {
	return (t > 0.0f) ? powf(t, (float)(3.0 / 2.2)) : 0.0f;
}

/************************************************/

//! Convert sRGB <-> XYZ
static Vec4f_t ConvertRGBtoXYZ(const Vec4f_t *x) {
	Vec4f_t t;
	float R = RGBtoLinearRGB(x->f32[0]);
	float G = RGBtoLinearRGB(x->f32[1]);
	float B = RGBtoLinearRGB(x->f32[2]);
	t.f32[0] = 0.412453f*R + 0.357580f*G + 0.180423f*B;
	t.f32[1] = 0.212671f*R + 0.715160f*G + 0.072169f*B;
	t.f32[2] = 0.019334f*R + 0.119193f*G + 0.950227f*B;
	t.f32[3] = x->f32[3];
	return t;
}
static Vec4f_t ConvertXYZtoRGB(const Vec4f_t *x) {
	Vec4f_t t;
	t.f32[0] = LinearRGBtoRGB( 3.240440f*x->f32[0] - 1.537001f*x->f32[1] - 0.499315f*x->f32[2]);
	t.f32[1] = LinearRGBtoRGB(-0.969205f*x->f32[0] + 1.875807f*x->f32[1] + 0.042506f*x->f32[2]);
	t.f32[2] = LinearRGBtoRGB( 0.055641f*x->f32[0] - 0.204021f*x->f32[1] + 1.057208f*x->f32[2]);
	t.f32[3] = x->f32[3];
	return t;
}

//! Convert XYZ <-> LAB (assuming D65 lighting)
//! The range is set to be nominally between 0.0 and 1.0 instead of 0 and 100.
static float LABf(float t) {
	//! delta = 6/29
	const float a = 0.008856f; //! delta^3
	const float b = 7.787037f; //! 1/3 * delta^-2
	if(t > a) return cbrtf(t);
	else return (float)(4.0/29.0) + b*t;
}
static float LABfInv(float t) {
	const float a = 0.128419f; //! 3*delta^2
	if(t > (float)(6.0/29.0)) return t*t*t;
	else return (t - (float)(4.0/29.0))*a;
}
static Vec4f_t ConvertXYZtoLab(const Vec4f_t *x) {
	float Xz = LABf(x->f32[0] / 0.950489f);
	float Yz = LABf(x->f32[1]);
	float Zz = LABf(x->f32[2] / 1.08884f);
	Vec4f_t t;
	t.f32[0] = 1.16f*(Yz     ) - 0.16f;
	t.f32[1] = 5.00f*(Xz - Yz);
	t.f32[2] = 2.00f*(Yz - Zz);
	t.f32[3] = x->f32[3];
	return t;
}
static Vec4f_t ConvertLabToXYZ(const Vec4f_t *x) {
	float Lz = (x->f32[0] + 0.16f) / 1.16f;
	float az = (x->f32[1]) / 5.0f;
	float bz = (x->f32[2]) / 2.0f;
	Vec4f_t t;
	t.f32[0] = 0.950489f * LABfInv(Lz + az);
	t.f32[1] =             LABfInv(Lz);
	t.f32[2] = 1.08884f  * LABfInv(Lz - bz);
	t.f32[3] = x->f32[3];
	return t;
}

/************************************************/

//! Convert sRGB <-> LMS
static Vec4f_t ConvertRGBtoLMS(const Vec4f_t *x) {
	Vec4f_t t;
	float R = RGBtoLinearRGB(x->f32[0]);
	float G = RGBtoLinearRGB(x->f32[1]);
	float B = RGBtoLinearRGB(x->f32[2]);
	t.f32[0] = 0.412221f*R + 0.536333f*G + 0.051446f*B;
	t.f32[1] = 0.211903f*R + 0.680700f*G + 0.107397f*B;
	t.f32[2] = 0.088302f*R + 0.281719f*G + 0.629979f*B;
	t.f32[3] = x->f32[3];
	return t;
}
static Vec4f_t ConvertLMStoRGB(const Vec4f_t *x) {
	Vec4f_t t;
	t.f32[0] = LinearRGBtoRGB( 4.076744f*x->f32[0] - 3.307714f*x->f32[1] + 0.230970f*x->f32[2]);
	t.f32[1] = LinearRGBtoRGB(-1.268435f*x->f32[0] + 2.609754f*x->f32[1] - 0.341319f*x->f32[2]);
	t.f32[2] = LinearRGBtoRGB(-0.004194f*x->f32[0] - 0.703420f*x->f32[1] + 1.707614f*x->f32[2]);
	t.f32[3] = x->f32[3];
	return t;
}

//! Convert LMS <-> ICtCp (using HLG transfer function from ARIB STD-B67)
static Vec4f_t ConvertLMStoICtCp(const Vec4f_t *x) {
	Vec4f_t t;
	float L = 0.5f * sqrtf(x->f32[0]);
	float M = 0.5f * sqrtf(x->f32[1]);
	float S = 0.5f * sqrtf(x->f32[2]);
	t.f32[0] = 0.500000f*L + 0.500000f*M;
	t.f32[1] = 0.885010f*L - 1.822510f*M + 0.937500f*S;
	t.f32[2] = 2.319336f*L - 2.249023f*M - 0.070313f*S;
	t.f32[3] = x->f32[3];
	return t;
}
static Vec4f_t ConvertICtCpToLMS(const Vec4f_t *x) {
	Vec4f_t t;
	float Lp = 2.0f * (x->f32[0] + 0.015719f*x->f32[1] + 0.209581f*x->f32[2]);
	float Mp = 2.0f * (x->f32[0] - 0.015719f*x->f32[1] - 0.209581f*x->f32[2]);
	float Sp = 2.0f * (x->f32[0] + 1.021271f*x->f32[1] - 0.605274f*x->f32[2]);
	t.f32[0] = Lp*Lp;
	t.f32[1] = Mp*Mp;
	t.f32[2] = Sp*Sp;
	t.f32[3] = x->f32[3];
	return t;
}

//! Convert LMS <-> Oklab (assuming D65 lighting)
static Vec4f_t ConvertLMStoOklab(const Vec4f_t *x) {
	Vec4f_t t;
	float L = cbrtf(x->f32[0]);
	float M = cbrtf(x->f32[1]);
	float S = cbrtf(x->f32[2]);
	t.f32[0] = 0.210454f*L + 0.793618f*M - 0.004072f*S;
	t.f32[1] = 1.977998f*L - 2.428592f*M + 0.450594f*S;
	t.f32[2] = 0.025904f*L + 0.782772f*M - 0.808676f*S;
	t.f32[3] = x->f32[3];
	return t;
}
static Vec4f_t ConvertOklabToLMS(const Vec4f_t *x) {
	Vec4f_t t;
	float Lp = x->f32[0] + 0.396338f*x->f32[1] + 0.215804f*x->f32[2];
	float Mp = x->f32[0] - 0.105561f*x->f32[1] - 0.063854f*x->f32[2];
	float Sp = x->f32[0] - 0.089484f*x->f32[1] - 1.291485f*x->f32[2];
	t.f32[0] = Lp*Lp*Lp;
	t.f32[1] = Mp*Mp*Mp;
	t.f32[2] = Sp*Sp*Sp;
	t.f32[3] = x->f32[3];
	return t;
}

/************************************************/

//! Convert RGBA to specified colourspace
static Vec4f_t ConvertToColourspace(const Vec4f_t *x, uint8_t Colourspace) {
	Vec4f_t Out, In = *x;

	//! Put colours into linear light space as needed
	switch(Colourspace) {
		case COLOURSPACE_RGB_LINEAR:
		case COLOURSPACE_RGB_PSY: {
			In.f32[0] = RGBtoLinearRGB(In.f32[0]);
			In.f32[1] = RGBtoLinearRGB(In.f32[1]);
			In.f32[2] = RGBtoLinearRGB(In.f32[2]);
		} break;
	}

	//! Apply transformation
	switch(Colourspace) {
		case COLOURSPACE_SRGB:
		case COLOURSPACE_RGB_LINEAR:
		case COLOURSPACE_RGB_PSY: {
			Out = In;
		} break;
		case COLOURSPACE_YCBCR:
		case COLOURSPACE_YCBCR_PSY: {
			Out.f32[0] =  0.2126f*In.f32[0] + 0.71520f*In.f32[1] + 0.0722f*In.f32[2];
			Out.f32[1] = -0.1146f*In.f32[0] - 0.38540f*In.f32[1] + 0.5000f*In.f32[2];
			Out.f32[2] =  0.5000f*In.f32[0] - 0.45420f*In.f32[1] - 0.0458f*In.f32[2];
			Out.f32[3] = In.f32[3];
		} break;
		case COLOURSPACE_YCOCG:
		case COLOURSPACE_YCOCG_PSY: {
			Out.f32[0] =  0.25f*In.f32[0] + 0.5f*In.f32[1] + 0.25f*In.f32[2];
			Out.f32[1] =  0.50f*In.f32[0]                  - 0.50f*In.f32[2];
			Out.f32[2] = -0.25f*In.f32[0] + 0.5f*In.f32[1] - 0.25f*In.f32[2];
			Out.f32[3] = In.f32[3];
		} break;
		case COLOURSPACE_CIELAB: {
			Out = ConvertRGBtoXYZ(&In);
			Out = ConvertXYZtoLab(&Out);
		} break;
		case COLOURSPACE_ICTCP: {
			Out = ConvertRGBtoLMS(&In);
			Out = ConvertLMStoICtCp(&Out);
		} break;
		case COLOURSPACE_OKLAB: {
			Out = ConvertRGBtoLMS(&In);
			Out = ConvertLMStoOklab(&Out);
		} break;
	}

	//! Finally, apply non-linearity and weighting
	switch(Colourspace) {
		case COLOURSPACE_RGB_PSY: {
			//! In terms of importance, G > R > B.
			//! We also map through a cubic root similar to Lab
			Out.f32[0] = cbrtf(Out.f32[0]) * 0.8f; //! R
			Out.f32[1] = cbrtf(Out.f32[1]) * 1.0f; //! G
			Out.f32[2] = cbrtf(Out.f32[2]) * 0.5f; //! B
		} break;
		case COLOURSPACE_YCBCR_PSY: {
			//! Curving of the luma component gives better
			//! overall results with less banding artifacts.
			//! We further weight Cb down (worse blue-yellow),
			//! as the Y and Cr components are more important.
			Out.f32[0]  = RGBtoVisualRGB(Out.f32[0]); //! Y
			Out.f32[1] *= 0.5f; //! Cb
		} break;
		case COLOURSPACE_YCOCG_PSY: {
			//! The colour opponents are weird for YCoCg, meaning
			//! we can't really apply any weighting to chroma
			Out.f32[0] = RGBtoVisualRGB(Out.f32[0]);
		} break;
	}
	return Out;
}

//! Convert from specified colourspace to RGBA
static Vec4f_t ConvertFromColourspace(const Vec4f_t *x, uint8_t Colourspace) {
	Vec4f_t Out, In = *x;

	//! Undo non-linearity and weighting
	switch(Colourspace) {
		case COLOURSPACE_RGB_PSY: {
			In.f32[0] = powf(In.f32[0] / 0.8f, 3.0f);
			In.f32[1] = powf(In.f32[1] / 1.0f, 3.0f);
			In.f32[2] = powf(In.f32[2] / 0.5f, 3.0f);
		} break;
		case COLOURSPACE_YCBCR_PSY: {
			In.f32[0]  = VisualRGBtoRGB(In.f32[0]);
			In.f32[1] /= 0.5f;
		} break;
		case COLOURSPACE_YCOCG_PSY: {
			In.f32[0] = VisualRGBtoRGB(In.f32[0]);
		} break;
	}

	//! Undo transformation
	switch(Colourspace) {
		case COLOURSPACE_SRGB:
		case COLOURSPACE_RGB_LINEAR:
		case COLOURSPACE_RGB_PSY: {
			Out = In;
		} break;
		case COLOURSPACE_YCBCR:
		case COLOURSPACE_YCBCR_PSY: {
			Out.f32[0] = In.f32[0] - 0.000152f*In.f32[1] + 1.574765f*In.f32[2];
			Out.f32[1] = In.f32[0] - 0.187280f*In.f32[1] - 0.468125f*In.f32[2];
			Out.f32[2] = In.f32[0] + 1.855610f*In.f32[1] + 0.000106f*In.f32[2];
			Out.f32[3] = In.f32[3];
		} break;
		case COLOURSPACE_YCOCG:
		case COLOURSPACE_YCOCG_PSY: {
			Out.f32[0] = In.f32[0] + In.f32[1] - In.f32[2];
			Out.f32[1] = In.f32[0]             + In.f32[2];
			Out.f32[2] = In.f32[0] - In.f32[1] - In.f32[2];
			Out.f32[3] = In.f32[3];
		} break;
		case COLOURSPACE_CIELAB: {
			Out = ConvertLabToXYZ(&In);
			Out = ConvertXYZtoRGB(&Out);
		} break;
		case COLOURSPACE_ICTCP: {
			Out = ConvertICtCpToLMS(&In);
			Out = ConvertLMStoRGB(&Out);
		} break;
		case COLOURSPACE_OKLAB: {
			Out = ConvertOklabToLMS(&In);
			Out = ConvertLMStoRGB(&Out);
		} break;
	}

	//! Convert back to sRGB
	switch(Colourspace) {
		case COLOURSPACE_RGB_LINEAR:
		case COLOURSPACE_RGB_PSY: {
			Out.f32[0] = LinearRGBtoRGB(Out.f32[0]);
			Out.f32[1] = LinearRGBtoRGB(Out.f32[1]);
			Out.f32[2] = LinearRGBtoRGB(Out.f32[2]);
		} break;
	}
	return Out;
}

/************************************************/
//! EOF
/************************************************/

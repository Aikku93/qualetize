# Qualetize
Tile-based Image Quantization Tool

This project supersedes the older [tilequant tool](https://github.com/Aikku93/tilequant).

## Purpose
This tool is mostly meant for GBA/NDS graphics, where each 'tile' can use one of many palettes. However, it can be adapted to just about any use (for example, custom formats).

## Getting started
Run `make` to build the tool and dynamic interface (.dll or .so, depending on the platform), then call `qualetize Input.bmp Output.bmp -npal:(no. of palettes) -cols:(entries/palette)` (eg. `qualetize Input.bmp Output.bmp -npal:16 -cols:16` to use all sixteen 16-colour GBA palettes).

There are a lot of supported options (eg. dithering options, colourspaces, etc.); run `qualetize` without any arguments for full dcoumentation.

Note that the dynamic library will potentially export many "internal" functions. If someone comes up with a clean way to only export the `Qualetize()` function from Qualetize.c, a PR would be greatly appreciated.

## Possible issues/todo

* Banding artifacts can be pretty rough around tile boundaries. This can be partially alleviated by matching `-ditherin:` with `-ditherout:`, but this may not be ideal.
* Tiles with just one or two colours may separate into a completely separate palette, even when they could be merged into another with no PSNR loss.

## Examples

All conversions performed using default settings plus `-ditherin:floyd` and `-npal:X`.

| Palettes | Result |
| - | - |
| Baseline truth | ![Baseline truth](/examples/cat.png?raw=true) |
| `-npal:1` | ![1 palette](/examples/cat-q1.bmp?raw=true) |
| `-npal:2` | ![2 palettes](/examples/cat-q2.bmp?raw=true) |
| `-npal:4` | ![4 palettes](/examples/cat-q4.bmp?raw=true) |
| `-npal:8` | ![8 palettes](/examples/cat-q8.bmp?raw=true) |
| `-npal:16` | ![16 palettes](/examples/cat-q16.bmp?raw=true) |

## Authors
* **Ruben Nunez** - *Initial work* - [Aikku93](https://github.com/Aikku93)
* **Marco Köpcke** - *Modifications and motivation for DLL interface* - [Parakoopa](https://github.com/Parakoopa)
* **zvezdochiot** - *Code and git cleanup* - [zvezdochiot](https://github.com/zvezdochiot)

# Bundled protocol font

`LiberationSans-Regular.ttf` is an **unmodified** Liberation Sans font, redistributed under SIL OFL 1.1. Its copyright and complete license are in `LICENSE-Liberation.txt`. The font remains under OFL; Collagene's independently implemented application code remains MIT.

Source: local distribution file `/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf`, corresponding font copyright/license from `/usr/share/doc/fonts-liberation2/copyright` (font section only; no Debian packaging code). Upstream: https://github.com/liberationfonts/liberation-fonts.

SHA-256: `8d91388f1d3604b3b8ae0e3ee2d140e50cd6122f9214514f4aca772540a4076d`.

The embedded PDF exporter reads only these compiled-in bytes through `include_bytes!`; it never searches host fonts. The whole original font is embedded in each standalone PDF, with a CID-to-glyph and ToUnicode map. No modified font or subset is distributed under a reserved original font name.

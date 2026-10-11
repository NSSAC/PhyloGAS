#!/bin/bash
# Export poster_project/poster.svg to a print-ready 48x24 inch PDF.
#
# Uses a headless Chrome/Chromium to rasterise nothing and embed the fonts as
# subsets, so the result is still vector and prints at any size.
#
#   bash scripts/export_poster_pdf.sh [in.svg] [out.pdf]
#
# NOTE the `body>svg` selector below. A bare `svg{...}` rule also matches the
# nested <svg> elements that hold each panel and blows them up to full page.
set -euo pipefail
SVG="${1:-poster_project/poster.svg}"
PDF="${2:-poster_project/PhyloGAS_poster_48x24.pdf}"

CHROME=""
for c in "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
         "/Applications/Chromium.app/Contents/MacOS/Chromium" \
         "$(command -v google-chrome || true)" "$(command -v chromium || true)"; do
    [ -x "$c" ] && { CHROME="$c"; break; }
done
[ -n "$CHROME" ] || { echo "No Chrome/Chromium found. Open $SVG in a browser and print to PDF."; exit 1; }

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
{
  printf '%s' '<!doctype html><html><head><meta charset="utf-8"><style>'
  printf '%s' '@page{size:48in 24in;margin:0}html,body{margin:0;padding:0;background:#fff}'
  printf '%s' 'body>svg{display:block;width:48in;height:24in}</style></head><body>'
  cat "$SVG"
  printf '%s' '</body></html>'
} > "$TMP/print.html"

"$CHROME" --headless --disable-gpu --no-pdf-header-footer \
    --print-to-pdf="$(cd "$(dirname "$PDF")" && pwd)/$(basename "$PDF")" \
    "file://$TMP/print.html" >/dev/null 2>&1
echo "-> $PDF"

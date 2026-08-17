# build/render_manual.R
#
# Renders the static user manual to www/manual.pdf at Docker build time.
# PDF (real LaTeX typesetting) is used instead of HTML+MathJax because the
# HTML build's Greek letters/formulas depended on a MathJax CDN load that
# doesn't render reliably in all viewing contexts; PDF math is baked in at
# render time and always displays correctly. Shiny automatically serves
# anything under www/ as a static asset, so the download needs no
# server-side handler. The manual is documentation, not user data, so
# (unlike analysis reports) it is legitimately baked into the image rather
# than regenerated per session.

library(tinytex)  # ensures the TinyTeX bin dir is on PATH for pandoc's pdflatex call

dir.create("www", showWarnings = FALSE)
rmarkdown::render(
  input = "manual/manual.Rmd",
  output_file = "manual.pdf",
  output_dir = "www",
  quiet = TRUE
)
message("Manual rendered to www/manual.pdf")

# build/render_manual.R
#
# Renders the static user manual to www/manual.html at Docker build time.
# Shiny automatically serves anything under www/ as a static asset, so the
# rendered manual needs no server-side download handler -- a plain <a href>
# download link in the UI is sufficient. The manual is documentation, not
# user data, so (unlike analysis reports) it is legitimately baked into the
# image rather than regenerated per session.

dir.create("www", showWarnings = FALSE)
rmarkdown::render(
  input = "manual/manual.Rmd",
  output_file = "manual.html",
  output_dir = "www",
  quiet = TRUE
)
message("Manual rendered to www/manual.html")

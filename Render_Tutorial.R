
# Render the Rmd directly to a LaTeX-based PDF:
rmarkdown::render(
  input  = "Tutorial Closeness Measure.Rmd",
  output_format = "pdf_document",
  output_file   = "Tutorial-Closness-Measure.pdf")


rmarkdown::render(
  input  = "Tutorial Closeness Measure.Rmd",
  output_format = "html_document",
  output_file   = "Tutorial-Closness-Measure.html"
)



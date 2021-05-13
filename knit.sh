#!/bin/zsh
echo '\nRendering PDF document'
R -e 'rmarkdown::render("draft/backarc.Rmd", output_format = "pdf_document")'
echo '\nRendering Word document'
R -e 'rmarkdown::render("draft/backarc.Rmd", output_format = "word_document")'
echo '\nRendering HTML document'
R -e 'rmarkdown::render("draft/backarc.Rmd", output_format = "html_document")'
open draft/backarc.pdf
open draft/backarc.docx
open draft/backarc.html
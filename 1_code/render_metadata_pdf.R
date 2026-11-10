# ---
# title: Render metadata README to PDF
# notes: Run from the project root (source() it or Run in RStudio). Writes
#   3_output/README_metadata.pdf next to the .md. Needs the Quarto CLI on PATH.
#   --to typst avoids needing a LaTeX install (Quarto's default pdf output does).
# ---

local({
  md_file <- "3_output/README_metadata.md"
  stopifnot("Run from the project root: 3_output/README_metadata.md not found" = file.exists(md_file))

  quarto <- Sys.which("quarto")
  stopifnot("quarto not found on PATH" = nzchar(quarto))

  status <- system2(quarto, c("render", md_file, "--to", "typst"))
  stopifnot("quarto render failed" = identical(status, 0L))

  message("Saved: ", sub("\\.md$", ".pdf", md_file))
})

# install.packages("DiagrammeR")
# install.packages("DiagrammeRsvg")
# install.packages("rsvg")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

diagram_code <- "
digraph G {
  rankdir=TB;
  node [shape=box, fontsize=10, fontname=Helvetica, style=filled, fillcolor=lightblue];

  OccurrenceData    [label=\"Occurrence Data\\n(CRAG, NCRMP, OBIS)\"]
  EnvVariables      [label=\"Environmental Variables\\n(Compile & Process Predictors)\"]
  PrepareDataset    [label=\"Prepare & Integrate Dataset\\n(Combine Occurrences & Env. Data)\"]
  RunMaxent         [label=\"Run Maxent\\n(Model Training)\"]
  ModelOutput       [label=\"Model Output\\n(Suitability Maps, Var. Importance)\"]
  ValidationNPS     [label=\"Validation with NPS\\n(Evaluate Predictions)\"]

  OccurrenceData -> EnvVariables -> PrepareDataset -> RunMaxent -> ModelOutput -> ValidationNPS;
}
"

g <- grViz(diagram_code)
svg_code <- export_svg(g)
rsvg_png(charToRaw(svg_code), file = "output/flowchart.png", width = 500)

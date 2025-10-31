rm(list=ls())
library(data.table)
library(purrr)
library(plyr)
library(haven)
library(tidyr)
library(readr)
library(ggplot2)
library(GGally)
library(ggsci)
library(reshape2)
library(survival)
library(dplyr)
library(stringr)
library(networkD3)
library(ggalluvial)


# setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_1997_120/wwanbing/output")
# BCR_1997 <- read_csv("Breast_SanofiU_1997_120_BCR_10APR2025.csv")
# 
# setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Breast_SanofiU_2000_118/wwanbing/output")
# BCR_2000 <- read_csv("Breast_SanofiU_2000_118_BCR_21APR2025.csv")

setwd("~/Desktop/Zhou Lab/Breast Cancer PRO project/Analysis/wwanbing/output/TTFD+Covariates")
Merged_TTFD <- read_csv("Merged_PRO_TTFD_01OCT2025.csv")

# Merge TTFD
# Select all variables starting with TTFD10, TTFD10P, censored, censored10P
# vars_to_merge <- colnames(Merged_TTFD)[grepl("^(TTFD10|TTFD10P|censored|censored10P)_", colnames(Merged_TTFD))]
# 
# vars_to_merge <- c("UID", vars_to_merge)
# 
# symptoms_TTFD <- Merged_TTFD[, vars_to_merge]
# 
# # Merge into BCR_1997 and BCR_2000
# # df_1997 <- left_join(BCR_1997, symptoms_TTFD, by = "UID")
# # df_2000 <- left_join(BCR_2000, symptoms_TTFD, by = "UID")
# 
# # suffix
# all_vars <- colnames(df_1997)
# symptom_vars <- all_vars[grepl("^TTFD10_", all_vars)]
# symptom_suffixes <- gsub("^TTFD10_", "", symptom_vars)
# unique_symptoms <- unique(symptom_suffixes)
# 
# # Q2 Heatmap: Symptoms (censor variable) × Relapse Type. 
# #   Color: % of patients with symptom deterioration; Rows: Symptoms; Columns: Relapse types 
# #   (all, distant/local/regional) 
# 
# df_1997$Year <- "1997"
# df_2000$Year <- "2000"
# df <- bind_rows(df_1997, df_2000)
# 
# all_vars <- colnames(df_1997)
# symptom_vars <- all_vars[grepl("^censored_", all_vars)]
# symptom_suffixes <- gsub("^censored_", "", symptom_vars)
# 
# 
# ### sankey plot
# ## 1) Column sets
# ttfd_cols     <- grep("^TTFD10_",  names(df), value = TRUE)
# censored_cols <- grep("^censored_", names(df), value = TRUE)
# 
# suf_ttfd     <- sub("^TTFD10_",  "", ttfd_cols)
# suf_censored <- sub("^censored_", "", censored_cols)
# sym_suf      <- intersect(suf_ttfd, suf_censored)
# 
# ttfd_cols     <- paste0("TTFD10_",  sym_suf)
# censored_cols <- paste0("censored_", sym_suf)
# 
# ## 2) Keep only relapsed patients
# base <- df %>% filter(BCRYN == 1)
# 
# ## 3) Wide-to-long (fix: only pivot TTFD/censored columns; exclude UID)
# long <- base %>%
#   select(UID, all_of(ttfd_cols), all_of(censored_cols)) %>%
#   pivot_longer(
#     cols      = all_of(c(ttfd_cols, censored_cols)),  # ✅ do not use everything()
#     names_to  = "var",
#     values_to = "val"
#   ) %>%
#   mutate(type    = if_else(grepl("^TTFD10_", var), "ttfd", "cens"),
#          symptom = sub("^(TTFD10_|censored_)", "", var)) %>%
#   select(-var) %>%
#   pivot_wider(names_from = type, values_from = val)   # get: UID, symptom, ttfd, cens

# Keep only relapsed patients and valid PFSDY
sankey_base <- Merged_TTFD %>%
  filter(PFS == 1 & PFSDY != -999) # n=892

# Find all symptom TTFD and censored columns (e.g., TTFD10_FA, censored_FA, etc.)
ttfd_cols     <- grep("^TTFD", names(sankey_base), value = TRUE)
censored_cols <- grep("^censored", names(sankey_base), value = TRUE)

# Extract suffixes (symptom abbreviations)
suf_ttfd     <- sub("^TTFD10_?", "", ttfd_cols)
suf_censored <- sub("^censored_?", "", censored_cols)
sym_suf      <- intersect(suf_ttfd, suf_censored)

# Recompose corresponding column names
ttfd_cols     <- paste0("TTFD10_",  sym_suf)
censored_cols <- paste0("censored_", sym_suf)

# Wide-to-long
long <- sankey_base %>%
  select(UID, PFSDY, all_of(ttfd_cols), all_of(censored_cols)) %>%
  pivot_longer(cols = all_of(c(ttfd_cols, censored_cols)),
               names_to = "var", values_to = "val") %>%
  mutate(type = if_else(grepl("^TTFD10_", var), "ttfd", "cens"),
         symptom = sub("^(TTFD10_|censored_)", "", var)) %>%
  select(-var) %>%
  pivot_wider(names_from = type, values_from = val)

## 4) One row per person per symptom: keep only symptoms that actually occurred (cens==0) and have TTFD, and TTFD < PFSDY
long_clean <- long %>%
  filter(!is.na(ttfd), is.finite(ttfd), !is.na(cens), cens == 0, ttfd < PFSDY) # n=808

## 5) For each patient, generate a time-ordered sequence of symptoms
seq_list <- long_clean %>%
  group_by(UID) %>%
  arrange(ttfd, .by_group = TRUE) %>%
  summarise(path = list(symptom), .groups = "drop")

## 6) Expand to transition edges (A->B) and count frequencies
edges <- seq_list %>%
  transmute(pairs = map(path, ~{
    syms <- .x
    if (length(syms) >= 2) {
      tibble(from = head(syms, -1), to = tail(syms, -1))
    } else {
      tibble(from = character(0), to = character(0))
    }
  })) %>%
  unnest(pairs)

edge_counts <- edges %>%
  count(from, to, name = "value") %>%
  arrange(desc(value))

## 7) Append terminal edges: last symptom of each patient -> All relapse
last_to_relapse <- seq_list %>%
  transmute(from  = map_chr(path, ~ tail(.x, 1)),
            to    = "All relapse") %>%
  count(from, to, name = "value")

edge_counts <- bind_rows(edge_counts, last_to_relapse) %>%
  arrange(desc(value))

## (Optional) Show only top-N transitions to avoid clutter
TOP_N <- 55   # e.g., keep top 30 transitions; set large to keep all
edge_counts_plot <- edge_counts %>% slice_head(n = TOP_N)

## 8) Build Sankey data
nodes <- data.frame(name = unique(c(edge_counts_plot$from, edge_counts_plot$to)))

links <- edge_counts_plot %>%
  mutate(source = match(from, nodes$name) - 1,
         target = match(to,   nodes$name) - 1) %>%
  select(source, target, value)

## 9) Plot (interactive)
sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "source", Target = "target",
  Value = "value", NodeID = "name",
  sinksRight = TRUE, nodeWidth = 10, fontSize = 8
)

### Update
## 0) Dependencies (if not loaded above)
#library(purrr)
library(networkD3)
library(htmlwidgets)
library(webshot2)  # snapshot HTML to PNG
library(magick)    # convert PNG to TIFF (common for journals)

## --- Symptom abbreviation -> full name (for node labels) ---
sym_map <- c(
  FA="Fatigue", NV="Nausea/Vomit", PA="Pain", DY="Dyspnea",
  SL="Insomnia", AP="Appetite loss", CO="Constipation",
  DI="Diarrhea", FI="Financial", PF2="Physical", RF2="Role",
  EF="Emotional", CF="Cognitive", SF="Social", GHS="GHS/QOL"
)

## --- Palette: fixed color per symptom (replace with your palette if needed) ---
sym_palette <- c(
  FA="#1f77b4", NV="#ff7f0e", PA="#d62728", DY="#9467bd", SL="#8c564b",
  AP="#e377c2", CO="#bcbd22", DI="#7f7f7f", FI="#c49c94",
  PF2="#636efa", RF2="#17becf", EF="#2ca02c", CF="#1a55FF", SF="#a55194", GHS="#ff9896"
)
sym_palette_ext <- c(sym_palette, `All relapse`="#bdbdbd")  # terminal node in neutral gray

## ========== 8) Build Sankey data (colored links + full-name nodes) ==========
# Helper: get starting symptom code (strip a possible step tag)
get_sym <- function(x) sub(" @\\d+$","", x)

# Node names: replace codes with full names (keep "All relapse" as-is)
node_names_raw <- unique(c(edge_counts_plot$from, edge_counts_plot$to))
node_names_full <- ifelse(node_names_raw == "All relapse",
                          "All relapse",
                          sym_map[get_sym(node_names_raw)])

nodes <- data.frame(
  # NodeID uses full names (as requested for node naming)
  name = node_names_full,
  # NodeGroup uses starting symptom code (for coloring)
  group = ifelse(node_names_raw == "All relapse", "All relapse", get_sym(node_names_raw)),
  stringsAsFactors = FALSE
)

# Links: group by source symptom so link color follows source color
links <- edge_counts_plot %>%
  mutate(
    source = match(ifelse(from=="All relapse","All relapse", sym_map[get_sym(from)]), nodes$name) - 1,
    target = match(ifelse(to  =="All relapse","All relapse", sym_map[get_sym(to  )]), nodes$name) - 1,
    group  = ifelse(from=="All relapse","All relapse", get_sym(from))
  ) %>%
  select(source, target, value, group)

# d3 color scale
all_groups <- unique(c(nodes$group, links$group))
palette_used <- sym_palette_ext[all_groups]   # select colors in order
colour_js <- htmlwidgets::JS(
  sprintf(
    "d3.scaleOrdinal().domain(%s).range(%s)",
    jsonlite::toJSON(all_groups, auto_unbox = TRUE),
    jsonlite::toJSON(unname(palette_used), auto_unbox = TRUE)
  )
)

edge_counts_plot$to[edge_counts_plot$to == "All relapse"] <- "Relapse"
nodes$name[nodes$name == "All relapse"] <- "Relapse"

## ========== 9) Plot (larger font + colored links/nodes) ==========
FIG_W <- 1100; FIG_H <- 500
NODE_W <- 25;  NODE_PAD <- 60; FONT_SZ <- 17  # enlarge fonts

p <- sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "source", Target = "target",
  Value = "value", NodeID = "name",
  NodeGroup = "group", LinkGroup = "group",
  colourScale = colour_js,
  sinksRight = TRUE,
  nodeWidth = NODE_W, nodePadding = NODE_PAD, fontSize = FONT_SZ,
  width = FIG_W, height = 200
)

# Further polish: link opacity + hover highlight; rounded node corners & borders
p <- htmlwidgets::onRender(p, "
function(el){
  const root = d3.select(el);
  root.selectAll('.link')
      .style('opacity', 0.55)
      .on('mouseover', function(){ d3.select(this).style('opacity', 0.9); })
      .on('mouseout',  function(){ d3.select(this).style('opacity', 0.55); })
      .append('title')
      .text(d => `${d.source.name} → ${d.target.name}: ${d.value}`);
  root.selectAll('.node rect')
      .attr('rx', 6).attr('ry', 6)
      .style('stroke', '#444').style('stroke-width', 0.6);
  root.selectAll('.node text')
      .style('font-weight', 600)
      .style('font-family', 'Helvetica, Arial, sans-serif');
}
")

p  # preview in RStudio

### Vertical text
## ========== 9) Plot (larger font + colored links/nodes) ==========
FIG_W  <- 1800      # canvas width (px)
FIG_H  <- 600       # canvas height (px)
NODE_W <- 15        # node bar width
NODE_PAD <- 30      # vertical spacing between node bars
FONT_SZ <- 30       # node label font size (px)
OPACITY <- 0.40     # default link opacity (0–1)
TEXT_SHIFT_Y <- -12 # fine-tune vertical position of rotated labels (negative=up, positive=down)
MARGIN_LEFT  <- 120 # left margin for vertical labels (px)
MARGIN_RIGHT <- 20  # right margin
LINK_WIDTH_SCALE <- 1.00  # overall link thickness scaling (1=original; >1 thicker; <1 thinner)

p <- sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "source", Target = "target",
  Value = "value", NodeID = "name",
  NodeGroup = "group", LinkGroup = "group",
  colourScale = colour_js,
  sinksRight = TRUE,
  nodeWidth = NODE_W, nodePadding = NODE_PAD, fontSize = FONT_SZ,
  width = FIG_W, height = FIG_H
)

## ====== onRender: extra left margin to avoid clipping + opacity + vertical big fonts + optional link width scaling ======
p <- htmlwidgets::onRender(p, sprintf("
function(el){
  const svg  = d3.select(el).select('svg');
  const root = svg.select('g');

  // --- 1) widen canvas and allow overflow to keep left-rotated labels visible ---
  const curW = +svg.attr('width') || %d;
  svg.attr('width', curW + %d + %d).style('overflow', 'visible');

  // move content to the right
  const prev = root.attr('transform') || '';
  root.attr('transform', 'translate(' + %d + ',0) ' + prev);

  // --- 2) link opacity + hover highlight ---
  root.selectAll('path.link')
      .style('opacity', %f)
      .style('stroke-opacity', %f)
      .on('mouseover', function(){ d3.select(this).style('opacity', 0.9).style('stroke-opacity', 0.9); })
      .on('mouseout',  function(){ d3.select(this).style('opacity', %f).style('stroke-opacity', %f); })
      .append('title')
      .text(d => d.source.name + ' → ' + d.target.name + ': ' + d.value);

  // --- 3) optional: uniformly scale link widths (based on d.dy) ---
  if (%f !== 1.0){
    root.selectAll('path.link')
        .style('stroke-width', d => Math.max(1, (d.dy || 1) * %f));
  }

  // --- 4) node styling ---
  root.selectAll('.node rect')
      .attr('rx', 6).attr('ry', 6)
      .style('stroke', '#444').style('stroke-width', 0.6);

  // --- 5) vertical labels: center, rotate, white outline, fixed font size ---
  root.selectAll('.node').each(function(){
    const g    = d3.select(this);
    const rect = g.select('rect');
    const w    = +rect.attr('width');
    const h    = +rect.attr('height');
    const txt  = g.select('text');

    txt.attr('x',  w/2)
       .attr('y',  h/2 + (%d))
       .attr('dy', 0)
       .style('text-anchor', 'middle')
       .attr('transform', 'rotate(-90 ' + (w/2) + ' ' + (h/2) + ')')
       .style('font-weight', 700)
       .style('font-family', 'Helvetica, Arial, sans-serif')
       .style('font-size', '%dpx')      // force font size
       .style('paint-order', 'stroke')
       .style('stroke', 'white')
       .style('stroke-width', '1.6px');
  });
}
", FIG_W, MARGIN_LEFT, MARGIN_RIGHT, MARGIN_LEFT,
OPACITY, OPACITY, OPACITY, OPACITY,
LINK_WIDTH_SCALE, LINK_WIDTH_SCALE,
TEXT_SHIFT_Y, FONT_SZ))

# p <- htmlwidgets::onRender(p, "
# function(el){
#   const root = d3.select(el);
# 
#   // keep your original link/rect styles
#   root.selectAll('.link')
#       .style('opacity', 0.4)
#       .on('mouseover', function(){ d3.select(this).style('opacity', 0.9); })
#       .on('mouseout',  function(){ d3.select(this).style('opacity', 0.4); })
#       .append('title')
#       .text(d => d.source.name + ' → ' + d.target.name + ': ' + d.value);
# 
#   root.selectAll('.node rect')
#       .attr('rx', 6).attr('ry', 6)
#       .style('stroke', '#444').style('stroke-width', 0.6);
# 
#   // === only change here: make all node labels vertical (rotate -90° around bar center) ===
#   root.selectAll('.node').each(function(){
#     const g = d3.select(this);
#     const rect = g.select('rect');
#     const w = +rect.attr('width');
#     const h = +rect.attr('height');
#     const txt = g.select('text');
# 
#     txt
#       .attr('x', w/2)              // center within the bar
#       .attr('y', h/2-12)
#       .attr('dy', 0)
#       .style('text-anchor', 'middle')
#       .attr('transform', 'rotate(-90 ' + (w/2) + ' ' + (h/2) + ')')
#       .style('font-weight', 600)
#       .style('font-family', 'Helvetica, Arial, sans-serif')
#       // thin white outline to improve readability over colored links (optional)
#       .style('paint-order', 'stroke')
#       .style('stroke', 'white')
#       .style('stroke-width', '1.5px');
#   });
# }
# ")
# 
# p
# 
# ### Opacity
# #OPACITY <- 0.4
# p <- htmlwidgets::onRender(p, sprintf("
# function(el){
#   var svg = d3.select(el).select('svg');
#   // set both stroke and overall opacity
#   svg.selectAll('path.link')
#      .style('stroke-opacity', %f)
#      .style('opacity', %f)
#      .on('mouseover', function(){ d3.select(this).style('stroke-opacity', 0.9).style('opacity', 0.9); })
#      .on('mouseout',  function(){ d3.select(this).style('stroke-opacity', %f).style('opacity', %f); });
# }
# ", OPACITY, OPACITY, OPACITY, OPACITY, FONT_SZ))
# 

p

# p <- htmlwidgets::onRender(p, "
# function(el){
#   const root = d3.select(el);
#   root.selectAll('path.link')
#       .append('title')
#       .text(d => `${d.source.name} → ${d.target.name}: ${d.value}`);
#   root.selectAll('.node rect')
#       .append('title')
#       .text(d => `${d.name}: ${Math.round(d.value !== undefined ? d.value : (d.dy || 0))}`);
# }
# ")
# 
# p

# Export: ensure vheight matches FIG_H (do not use 200)
htmlwidgets::saveWidget(p, file = "Sankey_VB_colored_fullname.html", selfcontained = TRUE)

# MARGIN_LEFT  <- 80   # leave extra space on the left for vertical labels
# MARGIN_RIGHT <- 20   # right margin

webshot2::webshot(
  "Sankey_VB_colored_fullname.html",
  file   = "Sankey_VB_colored_fullname.png",
  vwidth = FIG_W + MARGIN_LEFT + MARGIN_RIGHT,, vheight = FIG_H,   # ← key
  zoom   = 2
)

im <- magick::image_read("Sankey_VB_colored_fullname.png")
magick::image_write(im, path = "Sankey_VB_colored_fullname_600dpi.tiff",
                    format = "tiff", compression = "lzw")
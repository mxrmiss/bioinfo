#!/usr/bin/env Rscript
# =============================================================================
# pubtree_ggtree_advanced.R —— 进阶款·加宽留白（全名友好）
# 要点：
#   • 两侧留白显著加宽，右侧不再裁切长物种名
#   • 支持值优先 treedata@data，回退 node.label；数值更醒目
#   • 细树枝 + 虚线对齐；ladderize；静默非致命警告
# 依赖：ggplot2, ggtree, treeio, tidytree, ape, dplyr, stringr, readr, scales
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(ggtree); library(treeio); library(tidytree)
  library(ape); library(dplyr); library(stringr); library(readr); library(scales)
})

# ============================== 配置区 ===============================
CFG <- list(
  tree_file         = "supermatrix.contree",
  outgroup_name     = "Nematostella_vectensis",
  abbr_genus        = FALSE,       # ← 用全名：FALSE；若需缩写设 TRUE
  tip_label_italic  = TRUE,
  tip_label_size    = 3.1,
  tip_label_offset  = 0.0105,      # ← 标签与树干距离，略增以防拥挤
  # —— 支持值展示 ——
  show_support      = TRUE,
  support_format    = "auto",      # "auto" | "bootstrap" | "posterior" | "both"
  support_min       = 50,
  support_size      = 8,         # ← 节点数字更醒目
  support_nudge_x   = 0.002,
  # —— 树图与线宽 ——
  tree_layout       = "rectangular",
  tree_linewidth    = 0.35,        # ← 树枝更细
  # —— tiplab 对齐线（虚线） ——
  align_tiplab      = TRUE,
  tiplab_linetype   = 2,           # 2=dashed
  tiplab_linesize   = 0.35,
  # —— 比例尺 ——
  treescale_width   = 0.05,
  treescale_digits  = 2,
  treescale_linesize= 0.45,
  # —— 主题与留白（加宽两侧） ——
  theme_base_size   = 10,
  font_family       = NULL,
  expand_x_right    = 0.24,        # ← 右侧更宽
  expand_x_left     = 0.055,       # ← 左侧也更宽
  expand_y_factor   = 0.035,
  # —— 导出（整体更宽） ——
  export_basename   = "pubtree_adv_wide",
  export_width_cm   = 21,          # ← 20–22 皆可
  export_height_cm  = 23,
  export_png_dpi    = 600
)

# ============================== 工具函数 ===============================
say <- function(...) cat(sprintf("%s\n", paste0(..., collapse="")))
if (!exists("is.waive")) is.waive <- function(x) inherits(x, "waiver")
`%||%` <- function(a,b) if(!is.null(a)) a else b
cm2in <- function(cm) cm/2.54
has_cairo <- function() isTRUE(try(capabilities("cairo"), silent=TRUE))

abbr_binomial <- function(x){
  s <- gsub("_"," ", x); sp <- strsplit(s, "\\s+")
  vapply(sp, function(z){
    if(length(z)>=2) sprintf("%s. %s", substr(z[1],1,1), z[2]) else paste(z, collapse=" ")
  }, "")
}

# —— 稳健：treedata@data 优先，兼容 list-column；否则回退 node.label ——
extract_support <- function(tr, phy){
  as_num_vec <- function(x){
    if (is.null(x)) return(numeric(0))
    if (is.list(x)) {
      v <- vapply(x, function(z){
        if (length(z) == 0) return(NA_real_)
        suppressWarnings(as.numeric(z[[1]]))
      }, numeric(1))
      return(v)
    }
    suppressWarnings(as.numeric(x))
  }
  pick_first <- function(cands, nm) {
    hit <- cands[cands %in% nm]
    if (length(hit) > 0) hit[1] else NA_character_
  }

  boot <- rep(NA_real_, Nnode(phy))
  post <- rep(NA_real_, Nnode(phy))

  if (inherits(tr, "treedata")) {
    dat <- tryCatch(tr@data, error=function(e) NULL)
    if (is.data.frame(dat) && nrow(dat) > 0) {
      cols <- names(dat)
      node_col <- if ("node" %in% cols) "node" else NA_character_

      boot_candidates <- c("UFboot","UFBootstrap","bootstrap","bs","support","SHaLRT","aLRT")
      post_candidates <- c("posterior","prob","PP","pp")

      boot_col <- pick_first(boot_candidates, cols)
      post_col <- pick_first(post_candidates, cols)

      if (!is.na(node_col)) {
        nd_ids   <- (Ntip(phy)+1):(Ntip(phy)+Nnode(phy))
        node_vec <- as_num_vec(dat[[node_col]])

        if (!is.na(boot_col)) {
          boot_vec <- as_num_vec(dat[[boot_col]])
          idx <- match(nd_ids, node_vec)
          ok  <- !is.na(idx) & idx >= 1 & idx <= length(boot_vec)
          tmp <- rep(NA_real_, length(nd_ids)); tmp[ok] <- boot_vec[idx[ok]]
          boot <- tmp
        }
        if (!is.na(post_col)) {
          post_vec <- as_num_vec(dat[[post_col]])
          idx <- match(nd_ids, node_vec)
          ok  <- !is.na(idx) & idx >= 1 & idx <= length(post_vec)
          tmp <- rep(NA_real_, length(nd_ids)); tmp[ok] <- post_vec[idx[ok]]
          post <- tmp
        }
      }
    }
  }

  # 回退：node.label（支持 "95/0.99" 或单值 "95"）
  if (all(is.na(boot)) && all(is.na(post))) {
    nl <- phy$node.label
    if (!is.null(nl)) {
      s <- as.character(nl)
      two <- grepl("/", s, fixed=TRUE); two[is.na(two)] <- FALSE
      if (any(two)) {
        parts <- strsplit(s[two], "/", fixed=TRUE)
        b <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 1)))
        p <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
        idx <- which(two)
        boot[idx] <- b; post[idx] <- p
      }
      one <- (!two) & (!is.na(s)) & nzchar(s)
      if (any(one)) boot[one] <- suppressWarnings(as.numeric(s[one]))
    }
  }
  data.frame(boot=boot, post=post)
}

decide_show_support <- function(df, mode=CFG$support_format, thr=CFG$support_min){
  txt <- rep(NA_character_, nrow(df)); nz <- function(x)!is.na(x)
  thrb <- if (thr<=1) thr*100 else thr
  thrp <- if (thr>1)  thr/100 else thr
  if (mode=="auto"){
    if (any(nz(df$post))){
      keep <- nz(df$post) & df$post>=thrp
      if (any(keep, na.rm=TRUE)) txt[keep] <- sprintf("%.2f", df$post[keep])
    } else {
      keep <- nz(df$boot) & df$boot>=thrb
      if (any(keep, na.rm=TRUE)) txt[keep] <- sprintf("%.0f", df$boot[keep])
    }
  } else if (mode=="posterior"){
    keep <- nz(df$post) & df$post>=thrp
    if (any(keep, na.rm=TRUE)) txt[keep] <- sprintf("%.2f", df$post[keep])
  } else if (mode=="bootstrap"){
    keep <- nz(df$boot) & df$boot>=thrb
    if (any(keep, na.rm=TRUE)) txt[keep] <- sprintf("%.0f", df$boot[keep])
  } else if (mode=="both"){
    kb <- nz(df$boot) & df$boot>=thrb; kp <- nz(df$post) & df$post>=thrp; idx <- kb|kp
    if (any(idx, na.rm=TRUE)){
      both  <- kb & kp; onlyb <- kb & !kp; onlyp <- kp & !kb
      if (any(both,  na.rm=TRUE)) txt[both]  <- sprintf("%.0f/%.2f", df$boot[both], df$post[both])
      if (any(onlyb, na.rm=TRUE)) txt[onlyb] <- sprintf("%.0f",        df$boot[onlyb])
      if (any(onlyp, na.rm=TRUE)) txt[onlyp] <- sprintf("%.2f",        df$post[onlyp])
    }
  }
  txt
}

compute_treescale_x <- function(xrange, w){
  xmin <- xrange[1]; xmax <- xrange[2]; span <- xmax - xmin
  xmin + span*0.02 + w*0.5
}

# ============================== 读取·定根·梯形化·标签 ===============================
options(warn=-1); on.exit({ options(warn=0) }, add=TRUE)

if (!file.exists(CFG$tree_file)) stop(sprintf("[致命] 未找到树文件：%s", CFG$tree_file))
say("[OK] 选用树文件：", CFG$tree_file)

tr  <- tryCatch(treeio::read.iqtree(CFG$tree_file), error=function(e) ape::read.tree(CFG$tree_file))
phy <- if (inherits(tr, "treedata")) tr@phylo else tr
if (!inherits(phy,"phylo")) stop("[致命] 无法得到 phylo 对象。")

# 外群匹配（大小写/下划线稳健）
norm <- function(s) tolower(gsub("[ _]+","_", s))
tips_norm <- norm(phy$tip.label); want_norm <- norm(CFG$outgroup_name)
if (!(want_norm %in% tips_norm)) {
  cand <- phy$tip.label[grepl(want_norm, tips_norm, fixed=TRUE)]
  if (length(cand)==1) CFG$outgroup_name <- cand else stop(sprintf("[致命] 外群 '%s' 不在末端标签中。", CFG$outgroup_name))
}
say("[OK] Rooted at: ", CFG$outgroup_name)

phy <- ape::root(phy, outgroup=CFG$outgroup_name, resolve.root=TRUE)
phy <- ape::ladderize(phy, right=TRUE)

# 直接在 phy 上完成标签美化（全名或缩写）
phy$tip.label <- if (isTRUE(CFG$abbr_genus)) abbr_binomial(phy$tip.label) else gsub("_"," ", phy$tip.label)

# 支持值表
sup <- extract_support(tr, phy)
node_df <- tibble::tibble(node=(Ntip(phy)+1):(Ntip(phy)+Nnode(phy))) %>% bind_cols(sup)
node_df$support_txt <- if (isTRUE(CFG$show_support)) decide_show_support(node_df) else NA_character_

# ============================== 出图 ===============================
p <- ggtree(phy, layout=CFG$tree_layout, linewidth=CFG$tree_linewidth) +
  ggtree::geom_tiplab(offset=CFG$tip_label_offset, fontsize=CFG$tip_label_size,
                      align=CFG$align_tiplab,
                      linetype=CFG$tiplab_linetype,
                      linesize=CFG$tiplab_linesize,
                      italic=CFG$tip_label_italic)

# 支持值（若开启）
if (isTRUE(CFG$show_support) && any(!is.na(node_df$support_txt))){
  inner <- subset(p$data, !isTip)
  ann <- inner %>%
    dplyr::select(node, x, y) %>%
    dplyr::left_join(dplyr::select(node_df, node, support_txt), by="node") %>%
    dplyr::filter(!is.na(support_txt))
  if (nrow(ann) > 0){
    p <- p + geom_text(
      data = ann,
      aes(x = x + CFG$support_nudge_x, y = y, label = support_txt),
      size = CFG$support_size / .pt, vjust = -0.2
    )
  }
}

# 比例尺：控制小数位，并与树枝线宽协调
b0   <- suppressWarnings(ggplot_build(p))
xr   <- range(b0$layout$panel_params[[1]]$x.range); span <- diff(xr)
ts_w <- round(span * CFG$treescale_width, digits = CFG$treescale_digits)
ts_x <- compute_treescale_x(xr, ts_w)
p <- p + ggtree::geom_treescale(x=ts_x, y=0, width=ts_w, fontsize=2.6,
                                linesize=CFG$treescale_linesize)

# 主题 + 扩边（两侧更宽）
base_family <- if (is.null(CFG$font_family)) "sans" else CFG$font_family
p <- p + theme_classic(base_size=CFG$theme_base_size, base_family=base_family) +
  theme(axis.line=element_blank(), axis.ticks=element_blank(),
        axis.text=element_blank(), axis.title=element_blank(),
        plot.margin=margin(6, 22, 6, 16))  # 右侧再给一点空间

b1 <- suppressWarnings(ggplot_build(p))
xr2 <- range(b1$layout$panel_params[[1]]$x.range)
yr2 <- range(b1$layout$panel_params[[1]]$y.range)
xaddL <- diff(xr2)*CFG$expand_x_left; xaddR <- diff(xr2)*CFG$expand_x_right; yadd <- diff(yr2)*CFG$expand_y_factor
p <- p + coord_cartesian(xlim=c(xr2[1]-xaddL, xr2[2]+xaddR),
                         ylim=c(yr2[1]-yadd,   yr2[2]+yadd),
                         clip="off")

# ============================== 导出 ===============================
w <- cm2in(CFG$export_width_cm); h <- cm2in(CFG$export_height_cm)
pdf_file <- paste0(CFG$export_basename,".pdf")
png_file <- paste0(CFG$export_basename,".png")

if (has_cairo()) ggsave(pdf_file, p, width=w, height=h, device=cairo_pdf) else ggsave(pdf_file, p, width=w, height=h)
ggsave(png_file, p, width=w, height=h, dpi=CFG$export_png_dpi, device="png")

say("[OK] 导出完成：", pdf_file, " / ", png_file)
say("[DONE] 顶刊风系统发育树绘制完成。")
say("小贴士：若右侧仍略紧，可把 expand_x_right 调到 0.28；若想更紧致，回落到 0.20。")

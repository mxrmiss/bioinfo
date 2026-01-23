tool_readme.txt

# snp_rbh — Required tools & environment checks

本流程涉及 Bash + Python3 + R 三部分。下面列出必需工具、用途、以及“存在性验证命令”（直接复制运行）。

------------------------------------------------------------
A) 必需命令行工具
------------------------------------------------------------

1) diamond
用途：双向 blastp + 建库，生成 best hit 与 RBH pairs
验证：
diamond version

2) mafft
用途：对每对蛋白序列做 pairwise/2-seq alignment（脚本用 --auto，可并行）
验证：
mafft --version

3) perl
用途：运行 pal2nal.pl（将蛋白比对回译为 codon alignment）
验证：
perl -v

4) awk / sort / head / tr / basename / printf / xargs（常见 coreutils + util）
用途：01_rbh_diamond.sh 的过滤、排序、best hit、RBH 计算与并行
验证：
awk --version
sort --version
head --version
tr --version
basename --version
printf "ok\n"
xargs --version

5) gzip（可选但建议）
用途：读取/处理 .gz（00/04 会读 .gz；01 可能用到 gff.gz 只是推断 code，不强制解压）
验证：
gzip --version

------------------------------------------------------------
B) 必需解释器/运行时
------------------------------------------------------------

1) Python 3
用途：00_chr_rename_map.py、02_make_pair_fastas.py、04_codon_stats_and_chr.py
验证：
python3 --version

2) R（>= 3.6 建议）
用途：05_plot_density.R
验证：
Rscript --version

------------------------------------------------------------
C) 必需 R 包
------------------------------------------------------------

脚本 05_plot_density.R 需要：
- ggplot2
- dplyr
- tidyr
- scales

验证（一次性）：
Rscript -e 'pkgs<-c("ggplot2","dplyr","tidyr","scales"); ok<-sapply(pkgs, requireNamespace, quietly=TRUE); print(ok); if(!all(ok)) quit(status=1)'

------------------------------------------------------------
D) 推荐 R 包（可选）
------------------------------------------------------------

- ragg
用途：PNG 输出更稳定（脚本会自动优先用 ragg::agg_png；没有也能退回 cairo png）
验证：
Rscript -e 'cat("ragg:", requireNamespace("ragg", quietly=TRUE), "\n")'

------------------------------------------------------------
E) 关键脚本/文件存在性检查
------------------------------------------------------------

1) pal2nal.pl 是否存在
验证：
test -s scripts/pal2nal.pl && echo "OK pal2nal.pl" || echo "MISSING pal2nal.pl"

2) 输入数据是否存在（示例：你的当前布局）
验证：
test -s data/proteomes/Sinonovacula_constricta.faa && echo OK
test -s data/proteomes/Sinonovacula_rivularis.faa && echo OK
test -s data/cds/Sinonovacula_constricta.fna && echo OK
test -s data/cds/Sinonovacula_rivularis.fna && echo OK
test -s data/annotation/Sinonovacula_constricta.gff3.gz && echo OK
test -s data/annotation/Sinonovacula_rivularis.gff.gz && echo OK
test -s data/genomic/Sinonovacula_constricta.fa.gz && echo OK
test -s data/genomic/Sinonovacula_rivularis.fasta.gz && echo OK

------------------------------------------------------------
F) 一句提醒（避免常见坑）
------------------------------------------------------------
- 02_make_pair_fastas.py 当前实现用普通 open() 读取 proteome/cds，所以默认要求 .faa/.fna 为非 gz。
- 04_codon_stats_and_chr.py 只从 GFF 中读取 feature= mRNA 或 transcript，并解析属性里的 ID=xxx 作为 transcript ID。


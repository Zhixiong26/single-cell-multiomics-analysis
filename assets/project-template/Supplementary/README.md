# Bundled reference files / 内置参考文件

These files are copied into every generated project. `init_project.py` uses one of
them only when the intake declares nothing for that `references` key, and it records
which ones it used in `Report.md`. Nothing here is project-specific and nothing here
is written by a run.

| File | `references` key | Digest | Source |
|---|---|---|---|
| `hg38.canonical.chrom.sizes` | `chrom_sizes` | sha256 `2c6508c30122042a0237fa5e2c876aa3433799180af80419859d30b304e275e4` | UCSC hg38 canonical chromosome sizes |
| `ENCFF356LFX_GRCh38_blacklist.bed.gz` | `blacklist` | md5 `393688b4f06c9ce26165d47433dd8c37` | ENCODE GRCh38 exclusion list (ENCFF356LFX) |

The bundled defaults are human hg38 only. A project that declares another organism,
or a `genome_id` other than GRCh38/hg38, receives none of them and is told why: an
hg38 exclusion list applied to a mouse genome removes the wrong regions while
looking entirely healthy.

To use your own files, declare them in the intake's `references` block (or edit
`references` in `config/project.yaml`) and point `chrom_sizes`/`blacklist` at them.
The `blacklist_md5` that the ALLCools and VMR stages verify is only defaulted
together with the bundled blacklist; a blacklist of your own needs its own digest,
which you can read with `md5sum`. Changing `references` changes the input signature,
so a recorded full validation no longer applies.

## 中文

这些文件会随每个生成的项目一起复制。只有当 intake 未声明对应的 `references`
键时，`init_project.py` 才会使用其中之一，并把实际使用的项记录在 `Report.md`
中。此处内容与具体项目无关，也不会被任何运行改写。

内置默认值仅适用于人类 hg38。若项目声明了其他物种，或 `genome_id` 不是
GRCh38/hg38，则不会继承这些文件，并会说明原因：把 hg38 排除列表用在鼠基因组上
会移除错误的区域，而表面上看不出任何异常。

如需使用自己的文件，请在 intake 的 `references` 块中声明（或直接修改
`config/project.yaml` 的 `references`）。ALLCools 与 VMR 阶段校验的
`blacklist_md5` 只随内置 blacklist 一起默认；自备的 blacklist 需要自己的摘要，
可用 `md5sum` 读出。修改 `references` 会改变输入签名，已记录的全量校验随之失效。

# PSO-MDL

[English README](README.md)

## 项目简介

PSO-MDL通过粒子群优化选择分子描述符子空间，并用独立模型的等权平均名次对候选分子排序。本研究用于PBP2相关抗菌候选筛选。排名表示描述符空间中的相对优先级，不直接预测抗菌活性或靶点结合能力。

## 运行环境

MATLAB R2024a或更新版本、Statistics and Machine Learning Toolbox、默认JVM及可用的PNG/PDF图形导出环境。四组输入已包含计算好的描述符和身份信息，无需再次安装或运行Dragon。请为约10.1 GB输入及新生成的结果预留足够磁盘空间；正式训练计算量较大。

## 使用方法

下载**全部**项目文件，包括100个候选描述符分片，解压到任意可写文件夹，将MATLAB的“当前文件夹”设为项目根目录。首次运行无需修改绝对路径，无需setup、历史模型或预先生成的中间文件：

```matlab
step00_prepare_dataset
step01_preprocess_descriptors
step02_train_models
step03_rank_candidates
step04_generate_figures
```

| 步骤 | 读取内容 | 生成内容 |
|---|---|---|
| 00 | 下列四组原始输入 | 身份审计、固定背景/正式候选划分、描述符缺失统计 |
| 01 | Step 00结果与原始描述符 | 参考/背景矩阵、拟合的预处理参数 |
| 02 | Step 01矩阵与预处理元数据 | 从头训练的64个模型、优化轨迹与统计 |
| 03 | Step 00–02结果与候选描述符 | `candidate_ranking.tsv`、各模型名次及稳定性统计 |
| 04 | 本次Step 02/03结果 | Figure 3A–C和Figure 4A–E，每幅均输出PNG和PDF |

全部输出写入项目根目录下的`results/research/`，由代码自动创建。各Step读取已有上游结果，不会暗中加载历史模型或用模拟输入替代正式输入；Step 04绘制本次运行的结果。
如果要另做一次完整运行，在`functions/pso_mdl_config.m`中将`cfg.Output.Directory`改为新的目录即可。默认保护已有输出；其他输入路径和研究参数也集中在该配置文件中。

正式设置保持不变：16个参考分子、20,000个背景分子、306,922个正式候选、2,781个保留描述符，以及64个随机种子、每次3,000代PSO。P为非medoid接近程度百分位的算术平均（0–100），J = −RFP。PSO内部用未缩放比较键K = −RF sum(p_t)判断粒子优劣，导出J = K/(nRef−1)，在不改变数学目标的前提下保留历史浮点比较路径。

已在NIRO上以MATLAB R2024a、Intel MKL 2023.2（AVX2）和28个计算线程，从头精确复现全部64个模型及全量候选排名。默认配置请求相同线程数并记录实际运行环境。其他CPU、MATLAB或BLAS环境可运行流程，但不保证优化轨迹逐位一致；BLAS分支不同会显示警告。

## 数据

`data/private/`中只保留四组从头运行输入：

```text
data/private/
  candidate_descriptor_parts/   # 100个编号TXT及descriptor_parts.json
  candidate_manifest.tsv
  reference_descriptors.txt
  reference_manifest.tsv
```

必须完整保留全部分片和JSON清单，不要独立修改或重新编号。代码将分片作为同一张表流式读取，无需合并TXT，也不依赖项目外的数据目录。
`private`名称为兼容原路径而保留；这四组明确路径现已加入`.gitignore`白名单，其余数据与生成结果仍被忽略。项目不附带历史模型、旧排名、测试输入或历史论文结果包。

正版Dragon 7启用全部二维描述符、关闭三维描述符，采用默认的保存全部描述符设置，计算及导出阶段不设置描述符排除或过滤。参考和候选原始表均含3,855维描述符。随后由MATLAB质量控制排除缺失率超过5%的18维及近零/无效尺度的1,056维，保留2,781维。详见[数据说明](data/README.md)。

四组输入总计约10.1 GB，单文件均小于100 MiB，但有文件超过网页上传限制。分片不等于绕过GitHub单次推送或仓库体积约束。正式发布前仍须选择合适的Git数据分发方式，并确认第三方衍生数据的再分发许可，参见[大文件说明](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)及[仓库限制](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)。本次没有自动上传。

## 引用与许可证

请按[CITATION.cff](CITATION.cff)引用：Zhu Haoyu，**PSO-MDL，版本1.0.0**。[MIT许可证](LICENSE)仅覆盖作者原创代码，不授予第三方分子数据、MATLAB、Dragon或其他软件的权利。

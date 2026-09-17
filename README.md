# Single-cell Multiomics Analysis Skill

## 中文

这是一个不依赖特定课题数据的 Codex skill，用于从多样本 10x RNA 和逐细胞 indexed ALLC 构建 Scanpy/Harmony、MethSCAn、ALLCools 与 MethylVI 工作流。它通过 YAML/TSV 配置生成独立项目，在每个 Slurm 任务提交前检查资源，并维护双语 README、Report 与机器可读运行证据。

```bash
python scripts/init_project.py --intake intake.yaml --output /path/to/new-project
python scripts/validate_project.py --project /path/to/new-project
python scripts/plan_workflow.py --project /path/to/new-project --routes auto --run-id first-run
python scripts/submit_workflow.py --project /path/to/new-project --run-id first-run --dry-run
```

## English

This data-independent Codex skill generates and operates Scanpy/Harmony, MethSCAn, ALLCools, and MethylVI workflows from multi-sample 10x RNA and indexed per-cell ALLC inputs. YAML/TSV configuration, per-task resource inspection, bilingual documentation, and machine-readable run evidence are built in.

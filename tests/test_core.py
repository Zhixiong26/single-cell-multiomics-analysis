#!/usr/bin/env python3
"""Dependency-light tests for generation, validation, planning, and scheduling."""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from _common import WorkflowError, load_project  # noqa: E402
import bootstrap_environments as environment_bootstrap  # noqa: E402
from bootstrap_environments import bootstrap  # noqa: E402
from init_project import main as unused_init_main  # noqa: F401,E402
from inspect_resources import choose, enrich_nodes, parse_scontrol, parse_sinfo  # noqa: E402
from inspect_run import inspect as inspect_run  # noqa: E402
from plan_workflow import make_tasks  # noqa: E402
from submit_workflow import slurm_job_state  # noqa: E402
from validate_project import allc_cell_id, validate, validate_allc_record  # noqa: E402


class SkillTests(unittest.TestCase):
    def build_inputs(self, root: Path, mode: str = "paired", bad_checksum: bool = False) -> Path:
        data = root / "input"
        rna = data / "rna"
        allc = data / "allc"
        refs = data / "refs"
        rna.mkdir(parents=True)
        allc.mkdir(parents=True)
        refs.mkdir(parents=True)
        (rna / "barcodes.tsv").write_text("AAAC-1\nAAAG-1\n", encoding="utf-8")
        (rna / "features.tsv").write_text("g1\tGENE1\tGene Expression\ng2\tMT-X\tGene Expression\n", encoding="utf-8")
        (rna / "matrix.mtx").write_text("%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 3\n2 2 1\n", encoding="utf-8")
        allc_path = allc / "AAAC-1.allc.tsv.gz"
        with gzip.open(str(allc_path), "wt") as handle:
            handle.write("chr1\t1\t+\tCGN\t1\t2\t1\n")
        Path(str(allc_path) + ".tbi").touch()
        chrom = refs / "genome.chrom.sizes"
        chrom.write_text("chr1\t1000\n", encoding="utf-8")
        blacklist = refs / "blacklist.bed"
        blacklist.write_text("chr1\t100\t200\n", encoding="utf-8")
        annotation = data / "annotation.tsv"
        annotation.write_text("cell_id\tcell_type\nS1_AAAC-1\tTypeA\n", encoding="utf-8")
        import hashlib
        digest = hashlib.sha256(chrom.read_bytes()).hexdigest()
        sample = {
            "sample_id": "S1", "condition": "control", "batch": "B1", "include": 1,
            "rna_path": str(rna) if mode in {"rna", "paired"} else "",
            "rna_format": "10x_mtx" if mode in {"rna", "paired"} else "",
            "allc_root": str(allc) if mode in {"allc", "paired"} else "",
            "allc_glob": "*.allc.tsv.gz", "cell_id_prefix": "S1",
        }
        intake = {
            "schema_version": 1, "project_id": "fixture", "organism": "human",
            "references": {
                "genome_id": "fixture-build", "chrom_sizes": str(chrom),
                "chrom_sizes_sha256": "bad" if bad_checksum else digest,
                "blacklist": str(blacklist),
            },
            "samples": [sample],
            "annotation": {"path": str(annotation), "cell_id_column": "cell_id", "cell_type_column": "cell_type"},
            "environments": [{
                "stage": stage, "python": sys.executable,
                "executable": sys.executable, "version_command": sys.executable + " --version", "required": 1,
            } for stage in ("orchestrator", "scanpy_allcools", "methscan", "methylvi")],
            "analysis": {"task_commands": {"scanpy": [
                "/bin/sh", "-c",
                "touch {task_dir}/artifact; printf '%s\\n' '{{\"artifacts\":[\"{task_dir}/artifact\"]}}' > {task_dir}/task_outputs.json",
            ]}},
            "scheduler": {"backend": "local", "local": {"max_threads": 64, "max_memory": "256G"}},
        }
        intake_path = root / "intake.yaml"
        intake_path.write_text(json.dumps(intake), encoding="utf-8")
        return intake_path

    def generate(self, root: Path, mode: str = "paired", bad_checksum: bool = False) -> Path:
        intake = self.build_inputs(root, mode, bad_checksum)
        project = root / "project"
        subprocess.check_call([sys.executable, str(ROOT / "scripts" / "init_project.py"), "--intake", str(intake), "--output", str(project)])
        return project

    def test_modes_and_route_detection(self):
        expected = {
            "rna": (True, False), "allc": (False, True), "paired": (True, True),
        }
        for mode, flags in expected.items():
            with tempfile.TemporaryDirectory() as temp:
                project = self.generate(Path(temp), mode)
                result = validate(project)
                self.assertEqual(result["status"], "valid")
                self.assertEqual(result["routes"]["scanpy"], flags[0])
                self.assertEqual(result["routes"]["methscan_vmr"], flags[1])
                self.assertEqual(load_project(project)["schema_version"], 2)

    def test_allc_cell_ids_support_common_and_custom_names(self):
        row = {"sample_id": "lung-S11-m1", "cell_id_prefix": "lung-S11-m1",
               "allc_cell_id_regex": "", "allc_cell_id_replacement": ""}
        self.assertEqual(allc_cell_id(Path("AAACCTGAGAAACCAT-1.allc.tsv.gz"), row),
                         "lung-S11-m1_AAACCTGAGAAACCAT-1")
        row.update({"allc_cell_id_regex": r"allc_(?P<cell_id>.+)\.tsv\.gz$"})
        self.assertEqual(allc_cell_id(Path("allc_lung-S11-m1_sc10.tsv.gz"), row),
                         "lung-S11-m1_sc10")

    def test_full_validation_checks_every_allc_and_approved_dmr_labels(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            project = self.generate(root, "paired")
            allc_root = root / "input/allc"
            second = allc_root / "AAAG-1.allc.tsv.gz"
            with gzip.open(str(second), "wt") as handle:
                handle.write("chr1\t2\t+\tCGN\t0\t2\t1\n")
            Path(str(second) + ".tbi").touch()
            annotation = root / "input/annotation.tsv"
            annotation.write_text("cell_id\tcell_type\nS1_AAAC-1\tTypeA\nS1_AAAG-1\tTypeB\n")
            project_cfg = json.loads((project / "config/project.yaml").read_text())
            project_cfg["annotation"] = {"table": str(annotation), "profile": None,
                                             "review_status": "approved", "cell_id_column": "cell_id",
                                             "cell_type_column": "cell_type"}
            (project / "config/project.yaml").write_text(json.dumps(project_cfg))
            analysis = json.loads((project / "config/analysis.yaml").read_text())
            analysis["analysis"]["methscan"]["min_cells"] = 1
            (project / "config/analysis.yaml").write_text(json.dumps(analysis))
            result = validate(project, mode="full", workers=2)
            self.assertEqual(result["status"], "valid")
            self.assertTrue(result["routes"]["methscan_dmr"])
            self.assertEqual(len(result["allc_validation"]), 2)
            annotation.write_text("cell_id\tcell_type\nS1_AAAC-1\tTypeA\nS1_AAAG-1\tUnassigned\n")
            self.assertFalse(validate(project)["routes"]["methscan_dmr"])

    def test_bad_checksum_and_duplicate_sample(self):
        with tempfile.TemporaryDirectory() as temp:
            project = self.generate(Path(temp), "paired", bad_checksum=True)
            self.assertEqual(validate(project)["status"], "invalid")
        with tempfile.TemporaryDirectory() as temp:
            project = self.generate(Path(temp), "rna")
            samples = project / "config" / "samples.tsv"
            lines = samples.read_text().splitlines()
            samples.write_text("\n".join(lines + [lines[1]]) + "\n")
            with self.assertRaises(WorkflowError):
                validate(project)

    def test_duplicate_cell_id_and_missing_environment(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            project = self.generate(root, "paired")
            second = root / "input" / "allc-second"
            second.mkdir()
            source = root / "input" / "allc" / "AAAC-1.allc.tsv.gz"
            shutil.copy2(str(source), str(second / source.name))
            (second / (source.name + ".tbi")).touch()
            samples = project / "config" / "samples.tsv"
            with samples.open("a", encoding="utf-8") as handle:
                handle.write("S2\tcontrol\tB2\t1\t\t\t%s\t*.allc.tsv.gz\tS1\n" % second)
            result = validate(project)
            self.assertEqual(result["status"], "invalid")
            self.assertIn("ALLC-derived cell IDs", " ".join(result["errors"]))
        with tempfile.TemporaryDirectory() as temp:
            project = self.generate(Path(temp), "rna")
            environments = project / "config" / "environments.tsv"
            environments.write_text(
                "stage\tpython\texecutable\tversion_command\trequired\n"
                "orchestrator\t/missing/python\t/missing/python\t\t1\n",
                encoding="utf-8",
            )
            result = validate(project)
            self.assertEqual(result["status"], "invalid")
            self.assertIn("required executable is absent", " ".join(result["errors"]))

    def test_environment_bootstrap_plans_isolated_prefixes(self):
        with tempfile.TemporaryDirectory() as temp:
            project = self.generate(Path(temp), "allc")
            self.assertTrue((project / "tools/bootstrap_environments.py").is_file())
            self.assertTrue((project / "environment-specs/analysis-core.yaml").is_file())
            environments = project / "config" / "environments.tsv"
            environments.write_text(
                "stage\tpython\texecutable\tversion_command\trequired\n",
                encoding="utf-8",
            )
            result = bootstrap(project, execute=False, discover=False)
            self.assertEqual(result["status"], "planned")
            self.assertEqual(result["required_profiles"], ["analysis_core", "methscan", "methylvi"])
            self.assertTrue(all(item["action"] == "create" for item in result["actions"]))
            self.assertTrue(all("/.environments/" in item["prefix"] for item in result["actions"]))

    def test_environment_bootstrap_execute_updates_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            project = self.generate(Path(temp), "allc")
            environments = project / "config" / "environments.tsv"
            environments.write_text(
                "stage\tpython\texecutable\tversion_command\trequired\n",
                encoding="utf-8",
            )

            def verified(prefix, profile, timeout=180):
                return {
                    "profile": profile, "prefix": str(prefix), "valid": True,
                    "version": profile + " test", "command": [str(prefix / "bin/python")],
                }

            with mock.patch.object(environment_bootstrap, "find_manager", return_value="/bin/true"), \
                    mock.patch.object(environment_bootstrap, "check_prefix", side_effect=verified), \
                    mock.patch.object(environment_bootstrap.subprocess, "check_call") as create:
                result = bootstrap(project, execute=True, discover=False)
            self.assertEqual(result["status"], "complete")
            self.assertEqual(create.call_count, 3)
            rows = environments.read_text(encoding="utf-8")
            for stage in ("orchestrator", "scanpy_allcools", "methscan", "methylvi"):
                self.assertIn(stage, rows)
            self.assertTrue((project / ".workflow/environment-bootstrap/result.json").is_file())
            self.assertTrue((project / "Scripts/Environment/Report.md").is_file())

    def test_plan_dry_run_and_local_completion(self):
        with tempfile.TemporaryDirectory() as temp:
            project = self.generate(Path(temp), "rna")
            subprocess.check_call([sys.executable, str(ROOT / "scripts" / "plan_workflow.py"), "--project", str(project), "--routes", "auto", "--run-id", "unit"])
            subprocess.check_call([sys.executable, str(ROOT / "scripts" / "submit_workflow.py"), "--project", str(project), "--run-id", "unit"])
            code = subprocess.call([sys.executable, str(ROOT / "scripts" / "inspect_run.py"), "--project", str(project), "--run-id", "unit"])
            self.assertEqual(code, 0)
            self.assertTrue((project / ".workflow/runs/unit/workflow.COMPLETE").is_file())
            artifact = project / ".workflow/runs/unit/tasks/scanpy/artifact"
            artifact.unlink()
            self.assertNotEqual(inspect_run(project, "unit")["status"], "complete")
            subprocess.check_call([sys.executable, str(ROOT / "scripts" / "submit_workflow.py"), "--project", str(project), "--run-id", "unit"])
            resumed = json.loads((project / ".workflow/runs/unit/submissions.json").read_text())
            self.assertTrue((project / ".workflow/runs/unit/tasks/scanpy/artifact").is_file())
            self.assertTrue(any(item.get("resumed") for item in resumed))

    def test_inspect_never_completes_a_partial_plan_or_marker_only_task(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            run_dir = project / ".workflow/runs/partial"
            run_dir.mkdir(parents=True)
            tasks = [{"id": "task_%d" % index, "parameters": {}} for index in range(10)]
            (run_dir / "plan.json").write_text(json.dumps({
                "input_signature": "input", "code_signature": "code", "tasks": tasks,
            }))
            (run_dir / "submissions.json").write_text(json.dumps([{
                "task": "task_0", "job_id": "local_task_0", "status": "complete",
            }]))
            marker_dir = run_dir / "tasks/task_0"
            marker_dir.mkdir(parents=True)
            (marker_dir / "task.COMPLETE").touch()
            result = inspect_run(project, "partial")
            self.assertNotEqual(result["status"], "complete")
            self.assertEqual(len(result["tasks"]), 10)
            self.assertFalse((run_dir / "workflow.COMPLETE").exists())

    def test_explicit_downstream_routes_have_closed_dependencies(self):
        analysis = {"methscan": {"vmr_thresholds": [0.02]}, "methylvi": {"feature_targets": [10000]}}
        routes = {
            "scanpy": False, "methscan_vmr": False, "methscan_dmr": False,
            "allcools": False, "methylvi_allcools": True,
            "methylvi_vmr": False, "methylvi_vmr_dmr": False,
        }
        tasks = make_tasks(routes, analysis)
        identifiers = {item["id"] for item in tasks}
        self.assertIn("allcools_features", identifiers)
        self.assertIn("methylvi_allcools_10000", identifiers)
        self.assertTrue(all(set(item["dependencies"]).issubset(identifiers) for item in tasks))
        for route in routes:
            selected = {name: name == route for name in routes}
            planned = make_tasks(selected, analysis)
            known = {item["id"] for item in planned}
            self.assertTrue(known)
            self.assertTrue(all(set(item["dependencies"]).issubset(known) for item in planned))

    def test_context_and_nonhuman_contigs_are_configuration_driven(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            allc = root / "cell.allc.tsv.gz"
            with gzip.open(str(allc), "wt") as handle:
                handle.write("1\t1\t+\tCHN\t0\t2\t1\n")
            checked = validate_allc_record(allc, context="CHN", full=True)
            self.assertEqual(checked["records_checked"], 1)
            chrom = root / "genome.sizes"
            chrom.write_text("1\t1000\n")
            blacklist = root / "blacklist.bed"
            blacklist.write_text("")
            dmr = root / "dmr.tsv"
            dmr.write_text("1\t10\t20\tx\tx\tx\tx\t0.1\t0.6\tgroup_A\t0.001\tx\n")
            summary = root / "summary.tsv"
            summary.write_text(
                "status\tdmr_file\tcell_type_a\tcell_type_b\tcomparison\n"
                "complete\t%s\tTypeA\tTypeB\tA_vs_B\n" % dmr
            )
            output = root / "out"
            subprocess.check_call([
                sys.executable, str(ROOT / "assets/project-template/Scripts/Methylvi/vmr_dmr/01_prepare_all_unique_pooled_dmrs.py"),
                "--pairwise-summary", str(summary), "--blacklist", str(blacklist),
                "--chrom-sizes", str(chrom), "--output-dir", str(output), "--sort-threads", "1",
            ])
            self.assertTrue((output / "all_unique_hypo_DMRs.merged.bed").is_file())

    def test_slurm_requires_explicit_partition_allow_list(self):
        with tempfile.TemporaryDirectory() as temp:
            project = self.generate(Path(temp), "rna")
            scheduler = json.loads((project / "config/scheduler.yaml").read_text())
            scheduler["backend"] = "slurm"
            scheduler["partitions"] = []
            (project / "config/scheduler.yaml").write_text(json.dumps(scheduler))
            result = validate(project)
            self.assertEqual(result["status"], "invalid")
            self.assertIn("partitions", " ".join(result["errors"]))

    def test_slurm_reserved_and_observed_memory_are_distinct(self):
        fixtures = ROOT / "tests" / "fixtures"
        nodes = enrich_nodes(
            parse_sinfo((fixtures / "sinfo.txt").read_text()),
            parse_scontrol((fixtures / "scontrol.txt").read_text()), 4096,
        )
        node = [item for item in nodes if item["name"] == "node-a"][0]
        self.assertEqual(node["schedulable_memory_mb"], 83904)
        self.assertEqual(node["observed_free_memory_mb"], 70000)
        scheduler = {"partitions": ["cpu"], "allow_nodes": [], "exclude_nodes": []}
        profile = {"floor": {"cpus": 4, "memory": "16G"}, "target": {"cpus": 16, "memory": "64G", "time": "01:00:00"}, "ceiling": {"cpus": 24, "memory": "96G"}}
        recommendation = choose(nodes, scheduler, profile)
        self.assertEqual(recommendation["reference_node"], "node-a")
        too_large = {"floor": {"cpus": 30, "memory": "100G"}, "target": {"cpus": 32, "memory": "120G"}, "ceiling": {"cpus": 64, "memory": "256G"}}
        with self.assertRaises(WorkflowError):
            choose(nodes, scheduler, too_large)
        self.assertEqual([item for item in nodes if item["name"] == "node-c"][0]["state"], "down")

    def test_slurm_resume_distinguishes_active_terminal_and_unknown_jobs(self):
        with mock.patch("submit_workflow.subprocess.check_output", return_value=b"RUNNING\n"):
            self.assertEqual(slurm_job_state("101"), "RUNNING")
        with mock.patch("submit_workflow.subprocess.check_output", side_effect=[b"", b"FAILED|\n"]):
            self.assertEqual(slurm_job_state("102"), "FAILED")
        with mock.patch("submit_workflow.subprocess.check_output", side_effect=[b"", b""]):
            with self.assertRaises(WorkflowError):
                slurm_job_state("103")


if __name__ == "__main__":
    unittest.main()

"""Static HTML report generation."""

from __future__ import annotations

from pathlib import Path


def _html_table(headers: list[str], rows: list[list[str]]) -> str:
    th = "".join(f"<th>{h}</th>" for h in headers)
    body = ""
    for row in rows:
        td = "".join(f"<td>{cell}</td>" for cell in row)
        body += f"<tr>{td}</tr>\n"
    return f"<table>\n<thead><tr>{th}</tr></thead>\n<tbody>\n{body}</tbody>\n</table>"


def _html_page(title: str, body: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  body {{ font-family: sans-serif; margin: 2em; }}
  table {{ border-collapse: collapse; width: 100%; }}
  th, td {{ border: 1px solid #ccc; padding: 4px 8px; text-align: left; }}
  th {{ background: #f0f0f0; }}
  tr:nth-child(even) {{ background: #f9f9f9; }}
  h1, h2 {{ color: #333; }}
</style>
</head>
<body>
<h1>{title}</h1>
{body}
</body>
</html>
"""


def write_probes_report(output_dir: Path, summary: dict[str, int]) -> Path:
    rows = [[k, str(v)] for k, v in sorted(summary.items())]
    table = _html_table(["Target", "Candidates"], rows)
    html = _html_page("genoprobe — Probe Candidates", table)
    report = output_dir / "report.html"
    report.write_text(html, encoding="utf-8")
    return report


def write_screen_report(output_dir: Path, summary: dict[str, int]) -> Path:
    rows = [[k, str(v)] for k, v in sorted(summary.items())]
    table = _html_table(["Target", "Probes Passing Screen"], rows)
    html = _html_page("genoprobe — Off-target Screen", table)
    report = output_dir / "report.html"
    report.write_text(html, encoding="utf-8")
    return report


def write_panels_report(output_dir: Path, summary: dict[str, int]) -> Path:
    rows = [[k, str(v)] for k, v in sorted(summary.items())]
    table = _html_table(["Target", "Panel Probes"], rows)
    html = _html_page("genoprobe — Final Panels", table)
    report = output_dir / "report.html"
    report.write_text(html, encoding="utf-8")
    return report

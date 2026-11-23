# Contributing

* **English only** (code, comments, docs, test names, file names). CI runs `tools/check_english_only.py`.
* **Results go through the registry.** Add an entry to `results/registry.yaml` with `computation`,
  `optimality` and assumption ids; never write a label by hand. New assumptions go to
  `results/assumptions.yaml` with an honest status.
* **Status changes:** `computation`, `mathematical_status` and `assumptions` are recorded separately
  (docs/status-policy.md). Only people set `mathematical_status`; `independently-reviewed` and `proven`
  need a human review record. Anything above `conjectural` needs a certificate that both verifiers accept.
* **Provenance:** generate certificates with the experiment scripts so `provenance` is filled in.
* **Keep negative results.** Record failed approaches as `kind: negative-result` and explain them
  under `experiments/<topic>/failed/`.
* **Do not edit `legacy/`.** Fix defects in `src/riemann_research/` and add a regression test.
* **Separate commits** for translation or formatting, bug fixes, and new results.
* Workflow: short-lived feature branches merged to `main` by pull request; the maintainer merges manually;
  tag releases that papers cite. Expensive runs use the manual `research` workflow.
* Before opening a pull request run `ruff check .`, `pytest`, `tools/verify_certificates.py`,
  `tools/check_registry.py --check` and `tools/check_english_only.py` (the same as CI). Passing them does not make a
  result proven; reviewers check the mathematics separately.

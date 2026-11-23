# Repository settings (to be applied manually by the owner)

These settings cannot be stored in the repository. Nothing here has been applied; apply them in
the GitHub web interface after the repository is created.

## Branch protection (ruleset on `main`), solo maintainer

The repository has a single maintainer. GitHub does not let an author approve their own pull
request, so required approvals would block every merge. Use instead:

* Require a pull request before merging, with **0** required approvals.
* Do not require review from code owners. `.github/CODEOWNERS` still routes review requests and
  documents ownership; enable "Require review from Code Owners" only when a second maintainer joins.
* Require the status check `checks` (from `.github/workflows/ci.yml`) to pass.
* Block force pushes and branch deletion; do not allow bypass for automation.
* The maintainer merges manually after reading the diff and the CI result. That manual merge is the
  human approval; no workflow, bot or AI agent merges pull requests.

## GitHub Actions

* Workflow permissions: "Read repository contents and packages permissions" (read-only default).
* Do not allow GitHub Actions to create or approve pull requests.
* Require approval for workflow runs from outside collaborators (all of them, not only first-time).
* Allowed actions: GitHub-owned actions only. The workflows pin them to full commit SHAs
  (`actions/checkout` v7.0.1, `actions/setup-python` v7.1.0, `actions/upload-artifact` v7.0.2).
  When updating, resolve the new tag to its commit SHA and keep the version comment.
* `research.yml` runs only on manual dispatch and never writes to the repository.

## Security features

* Enable Dependabot alerts, secret scanning and push protection.
* Enable private vulnerability reporting.

## Human approval for mathematical statuses

Changes to `results/registry.yaml` that set `independently-reviewed` or `proven` must include the
review record of a human reviewer other than the maintainer, and are merged manually by the
maintainer. No workflow or tool edits these fields.

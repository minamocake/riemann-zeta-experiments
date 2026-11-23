# Known limitations (v0.1)

* **Test-function family.** Only centered boxes (and PSD mixtures of them). Their Fourier side is
  compactly supported, which makes the sign-only relaxation degenerate
  (`derivations/sign-only-degeneracy.md`). The 2020 pair-correlation benchmark uses a different class.
* **No new results.** The best RH-only value here (0.67763) relies on an unreviewed lemma and is
  below the pair-correlation benchmark 0.6792.
* **Counterfactual caps.** All caps-mode values above 0.68 depend on unproven hypotheses.
* **Floating-point optimisers.** SLSQP and SDP outputs are candidates. Only certificates are exact,
  and they are re-checked by two verifiers that share no code.
  SDP tolerances can hide large infeasibilities (see the cascade remark in the degeneracy note).
* **Not formally verified.** The deduction from certificates to zeta statements is written out but
  not machine-checked; Montgomery's theorem is cited, not re-proved.
* **Legacy code** keeps its known defects by design (`tests/test_legacy_known_defects.py`).
* **Nothing is reviewed.** Every certified result is `proof-pending`; the 8/3 lemma, the degeneracy
  proposition and the deduction in `derivations/what-a-certificate-establishes.md` await human review.
* **Provenance commit field.** Certificates were generated before the first commit, so their
  `provenance.commit` reads "unavailable"; regenerate in a clean checkout after the first commit.
* **Environment.** Only CPython 3.13 on Linux x86_64 is tested; the lock file pins versions but not hashes.
* **CI not yet run on GitHub.** Workflows pass actionlint and were simulated locally only.
* **Open warnings from the independent reviewer's script** (0 failures): no `--hash` pins in the lock
  file; push CI runs on `main` only (intended; pull requests are checked when opened or updated); a
  generic reminder about Simpson's rule (exact here, see the reviewer's report); the 8/3 lemma is
  proof-pending.
* **Placeholders to fill in before publishing:** author name in `CITATION.cff`; confirm the license.

# PROTEINS submission checklist — PepArena

Submit at https://submission.wiley.com (Research Exchange). Free-format submission is accepted.

## Which file do I open?

**`01_Manuscript/PepArena_Manuscript_with_figures.docx`** — this is the one to read. All six
main-text figures are embedded in place; Figures S1–S7 live in the Supporting Information PDF.

`PepArena_Manuscript.docx` deliberately has **no images in it**: Wiley wants the main document
as text plus a *Figure legends* list, with the figures uploaded as separate files. That file
looking like "just the captions" is correct, not a build failure.

## In this folder

| Folder | Upload as |
|---|---|
| `01_Manuscript/PepArena_Manuscript.docx` | Main Document |
| `01_Manuscript/PepArena_Manuscript_with_figures.docx` | review copy (figures inline) |
| `01_Manuscript/PepArena_Cover_Letter.docx` | Cover Letter (optional at PROTEINS) |
| `02_Figures/Figure_*.png` | one file per figure, 300 dpi PNG, named to match the Figure N legends; legend goes beneath each image at upload **and** as a list in the text |
| `03_Supporting_Information/PepArena_Supporting_Information.pdf` | Supplementary Material for Review (**PDF required**) |
| `04_Single_File/` | text + figures + full SI in one document — PROTEINS accepts this as the initial submission format |
| `zenodo/` | upload to Zenodo, then paste the version DOI into the manuscript |

## Before you click submit

- [ ] **ORCID** 0009-0002-4615-1448 attached to the submitting author (PROTEINS requires it)
- [ ] **Zenodo deposited**; the Data Availability Statement carries the literal placeholder
      `[DOI TO BE INSERTED AT DEPOSITION]`. Deposit, replace that string with the **version** DOI,
      and rebuild — the placeholder must not reach Wiley.
- [ ] **Preprint — there is none, and none is needed.** The cover letter asserts only that the
      work is unpublished and not under consideration, which is correct as it stands. If you want
      one, the zero-effort route is Wiley's **Under Review** service (Authorea), opted into during
      submission at https://www.authorea.com/inst/20386: Wiley posts it, so the "tell the
      editorial office" requirement is satisfied by the opt-in itself. The only obligation it
      creates is on acceptance — update the preprint to cite the published DOI, and cite only the
      published version yourself.
- [ ] Figures uploaded at highest available resolution (the PNGs are 300 dpi)
- [x] **Repository public** at https://github.com/Bhargava-Systems-Research-Inc/PepArena, which is
      what the Data Availability Statement points at. A checkout is 932 files / 29 MB; the working
      artifacts (ADCP poses, derived alignment copies, fetch caches) are untracked and live in the
      Zenodo deposit instead.
- [ ] Suggested reviewers, if you want to name any

## Already verified by the build (`manuscript/rebuild.sh`)

Nine gates abort the build on failure, so these are true of the files in this folder:

- `manifest_selfcheck.py` — the released manifest agrees with itself: sequence length against
  `peptide_len`, the non-canonical flag against the codes recorded, the length bin against the
  length, cyclic flag against chemistry, no duplicate identifiers.
- `track_a/finalize.py --check` — every Track A table re-hashed against the single scoring run it
  came from, so a table rebuilt on its own cannot sit beside one that was not.
- `track_b/worklist_gate.py` — each docking program's scored complexes lie on the worklist the
  paper reports, and every accepted complex was actually attempted, so coverage and accuracy
  always describe the same 183 complexes.
- `audit_rules.py` — the prose rules themselves, checked for the shapes that let a rule pass while
  verifying nothing: more than one capture group, a pattern matching in several places, two rules
  checking the same characters, or a literal number the pattern anchors on that nothing captures.
- `verify_facts.py` — every headline number re-derived from its run artifact and compared against
  `DATA_FACTS.md`.
- `check_prose.py` — 139 rules comparing sentences directly against the artifacts they came from,
  plus a structural check that the model named as leader really leads on both the conditional and
  the unconditional metric. A rule whose sentence has been reworded is reported as missing rather
  than passing quietly.
- `check_front_matter.py` — every number in the abstract and the cover letter also appears
  somewhere in the paper, which is the comparison a reviewer actually makes.
- `trace_numbers.py` — every numeric literal in Results, Discussion and the cover letter traced to
  an artifact. All three are at zero.
- `check_guidelines.py` — the PROTEINS author-guideline checks: abstract ≤250 words, running title
  <40 chars, 5–7 keywords, AMA references, required sections, every figure cited, ORCID and both
  affiliations present, Supporting Information built as PDF.

`make_figures.py` additionally refuses to run if two functions would write the same figure file,
or if one figure function is defined twice. Both guards earned their place: the first after
Figures S1 and S4 once showed something other than what their captions said, and the second after
four definitions of one function meant the last one silently drew Figure 5 from a superseded
leaderboard. The one-producer guard had itself been written and never wired in — by the time it
was switched on, three functions were writing `Figure_S1.png` and Figures 2, 3 and 6 were not
being generated at all.

## Open question worth one email

The guidelines say authors "must include the RRIDs in the list of keywords" for software and
databases. This work uses ~10 tools, which collides with the 5–7 keyword cap. Those two rules
conflict as written; ask proteinsadmin@wiley.com if you want cover.

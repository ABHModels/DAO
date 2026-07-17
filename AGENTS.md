# DAO Repository Safety Rules

These rules apply to every change, commit, tag, and remote operation in this
repository. Safety, reproducibility, and physical correctness take precedence
over convenience or a clean-looking history.

## 1. Repository roles

- `ABHModels/DAO:main` is the canonical public, stable, releasable branch.
- `DAO-Model/DAOv2:main` is a private mirror of the public stable branch. It is
  not an independent development line.
- Private integration work belongs on `origin/dev` or short-lived branches.
- Task branches should use `codex/<topic>` unless the user requests another
  name.
- Never allow the two `main` branches to accept unrelated direct development.

## 2. Remote safety

Before every fetch, push, tag publication, or release operation, run:

```bash
git status --short --branch
git remote -v
git branch -vv
```

Expected remote ownership:

- `origin` fetches from and pushes to `DAO-Model/DAOv2`.
- `abh` fetches from and pushes to `ABHModels/DAO`.

Stop if a remote fetch URL and push URL refer to different repositories, or if
`origin` pushes to `ABHModels/DAO`. Do not silently repair remote configuration;
report the mismatch and get user authorization first.

- Never configure one remote to push automatically to both the private and
  public repositories.
- Never use `git push --all` or `git push --tags` against the public repository.
- Push only the named branch and, for a release, the one intended tag.
- SSH (`git@github.com:ABHModels/DAO.git`) is the preferred fallback when the
  local HTTPS proxy is unavailable. Do not place tokens in remote URLs.
- `PUSH_NOTES.local.md` and all `*.local.md` files are local-only and must never
  be committed.
- A public push requires an explicit user request. A normal code-change request
  does not authorize publishing.

## 3. Branch and pull-request workflow

- Do not develop directly on `main`.
- Start public work from the current `abh/main`, not from a stale local branch.
- Use one branch per coherent change and keep commits reviewable.
- Prefer a pull request into public `main`; direct pushes to public `main` are
  reserved for an explicit user-requested release or emergency action.
- Before integrating, inspect divergence with `git log --left-right` and verify
  ancestry with `git merge-base --is-ancestor`.
- If branches diverge, preserve both histories with an inspected merge or a PR.
  Never resolve divergence by force-pushing or silently dropping commits.
- Never merge private history into the public repository if that history may
  contain secrets, unpublished material, or restricted data. Create a clean
  branch from `abh/main` and import only the approved patch or commits.

## 4. Required checks before a commit

- Preserve unrelated user changes; stage only files that belong to the task.
- Review both `git diff` and `git diff --cached`.
- Run `git diff --check`.
- Validate `CITATION.cff` as YAML whenever it changes.
- Compile and run the most relevant available tests.
- For scientific or numerical changes, require regression evidence appropriate
  to the affected physics: finite values, convergence, conservation checks,
  limiting cases, and comparison with a trusted baseline using documented
  tolerances.
- Treat `test_kernel_norm` as a diagnostic unless it contains explicit pass/fail
  assertions. Do not report it as an accuracy pass merely because it exits 0.
- If the full Cloudy/HEASoft build cannot run locally, state that limitation in
  the commit or PR summary; do not imply full validation.

## 5. Public release workflow

For a release `vX.Y.Z`:

1. Confirm the worktree is clean and fetch the current public `main` and tags.
2. Create `release/vX.Y.Z` from the current `abh/main`.
3. Import only the approved, tested changes.
4. Update `CHANGELOG.md` and `CITATION.cff`; if a `VERSION` file exists, update
   it too. All version strings and dates must agree.
5. Run the fast checks and the full physics/integration checks required by the
   changed modules.
6. Push the release branch and merge it through a reviewed PR whenever possible.
7. Create the annotated release tag only after the final public merge commit is
   known. Prefer a signed tag when signing is configured.
8. Verify that the tag peels to the exact intended commit, then push public
   `main` and only that tag.
9. Verify the remote branch and tag independently, then synchronize the private
   stable mirror.

Published version tags are immutable in practice:

- Never move, overwrite, or delete a published `v*` tag.
- Never use `--force` or `--force-with-lease` on public `main` or a release tag.
- If a release is wrong, publish the next patch version instead of rewriting the
  existing release.

## 6. Large files and generated data

- Source code, small test fixtures, parameter files, and curated reference data
  may be tracked.
- Runtime output, caches, temporary diagnostics, and reproducible intermediates
  must remain ignored.
- `paper_data/` is for curated reproducibility artifacts, not raw run directories.
- Any newly added file larger than 1 MiB requires explicit justification in the
  PR. Do not add files larger than 10 MiB to normal Git without user approval.
- Prefer Git LFS, a GitHub Release asset, or an archival data repository such as
  Zenodo for necessary large datasets and binary artifacts.
- Do not rewrite existing public history merely to migrate data without an
  explicit migration and backup plan.
- Before a public push, inspect newly tracked large files and confirm that no
  `results/`, `kernel/`, `comp_cs/`, logs, local notes, or credentials are staged.

## 7. GitHub protection policy

The public repository should protect `main` with a GitHub ruleset requiring:

- pull requests for normal changes;
- required CI/status checks;
- resolved review conversations;
- stale approval dismissal after new pushes;
- no branch deletion or force push;
- code-owner review when multiple maintainers are available.

A separate `v*` tag ruleset should restrict tag creation, updates, and deletion.
Enable secret-scanning push protection. Formal releases should be published as
immutable GitHub Releases when available.

## 8. Destructive-operation policy

- Do not run `git reset --hard`, destructive checkout/restore commands, history
  filters, or force pushes without explicit user authorization and a verified
  backup.
- Prefer revert commits for published mistakes.
- Before any history migration, create and verify a mirror or `git bundle`
  backup stored outside the working repository.

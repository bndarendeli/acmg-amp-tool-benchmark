# Push to GitHub

Use this directory as the repository root, not its parent source directory.

1. Run `python tools/validate_repository.py` from this directory.
2. Create an empty repository on GitHub without an initial README or .gitignore.
3. Replace `REPO_URL` below with your repository URL:

```bash
git init -b main
git add .
git status --short
git commit -m "Add manuscript data and analyses"
git remote add origin REPO_URL
git push -u origin main
```

No commit, remote or push has been created during preparation. Reconcile the Figure 3 policy discrepancy in `INPUT_VERIFICATION.md` before presenting this as a finalized manuscript release. Add the final article title, authors and DOI when available. Code and data licensing have not yet been selected. See `REPRODUCIBILITY.md` for verified execution scope.

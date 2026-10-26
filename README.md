# Visual Elliptic FEM

The application and its full experiment guide are in [`elliptic-lab/README.md`](elliptic-lab/README.md).

Run it locally with no package installation:

```bash
cd elliptic-lab
python3 -m http.server 4173 --bind 127.0.0.1 --directory dist
```

Open http://127.0.0.1:4173 in a browser. Start with **Transport layer**, then **Apply mesh-based value** to explore the artificial-diffusion experiment.

Time-dependent convection, diffusion, and convection–diffusion are also included. Try **Moving front** for transport, **Diffusing pulse** for consistent/lumped mass, and **Moving + diffusing** for a moving Gaussian. Use Play, single-step, or the time slider. The detailed README explains their limits and comparison conventions.

## Host on GitHub Pages

This application is static: no backend, build step, or API keys are required. The workflow in [`.github/workflows/deploy-pages.yml`](.github/workflows/deploy-pages.yml) publishes only `elliptic-lab/dist` and supports pushes to `main` or `master`, plus manual runs. All asset and worker paths are relative, so project URLs such as `https://YOUR_USERNAME.github.io/elliptic-lab/` work without a base-path change.

1. On GitHub, create an empty repository called `elliptic-lab`. A public repository supports GitHub Pages on GitHub Free. Leave the initial README, license, and .gitignore unchecked because the local project already has files.
2. From this workspace root (the folder containing this README and the `elliptic-lab` directory), run the commands below. Replace `YOUR_USERNAME` with your GitHub username. This folder already has a local Git repository and currently has no remote.

   ```bash
   git add .gitignore .github README.md elliptic-lab
   git commit -m "Add Elliptic Lab and GitHub Pages deployment"
   git branch -M main
   git remote add origin https://github.com/YOUR_USERNAME/elliptic-lab.git
   git push -u origin main
   ```

3. In the GitHub repository, open **Settings → Pages**. Under **Build and deployment**, set **Source** to **GitHub Actions**.
4. Open **Actions → Deploy Elliptic Lab to GitHub Pages → Run workflow** and run it on `main`. If the first push ran before Pages was enabled and failed, this manual run retries after configuration.
5. Wait for a successful deployment. The `github-pages` environment and deployment job show the actual site URL, normally `https://YOUR_USERNAME.github.io/elliptic-lab/`.

Later pushes to `main` automatically update the site. If you upload files through the GitHub website instead, preserve this workspace's layout: `.github/workflows/deploy-pages.yml` belongs at the repository root, and the app remains in `elliptic-lab/dist`. Uploading just the inner app folder will not include the workflow provided here.

GitHub Pages sites are ordinarily public websites; a private source repository does not by itself make the deployed site private. Private-repository Pages availability depends on the GitHub plan. See [GitHub's publishing-source documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site) and [custom workflow documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).

# CDR FEM FCT lab

The application and its full experiment guide are in [`CDR_visualizer/README.md`](CDR_visualizer/README.md).

Run it locally with no package installation:

```bash
cd CDR_visualizer
python3 -m http.server 4173 --bind 127.0.0.1 --directory dist
```

Open http://127.0.0.1:4173 in a browser. Start with **Transport layer**, then **Apply mesh-based value** to explore the artificial-diffusion experiment.

Time-dependent convection, diffusion, and convection–diffusion are also included. Try **Moving front** for transport, **Diffusing pulse** for consistent/lumped mass, and **Moving + diffusing** for a moving Gaussian. Use Play, single-step, or the time slider. The detailed README explains their limits and comparison conventions.

## Host on GitHub Pages

This application is static: no backend, build step, or API keys are required. The workflow in [`.github/workflows/deploy-pages.yml`](.github/workflows/deploy-pages.yml) publishes only `CDR_visualizer/dist` and supports pushes to `main` or `master`, plus manual runs. All asset and worker paths are relative, so project URLs such as `https://YOUR_USERNAME.github.io/CDR_visualizer/` work without a base-path change.

1. On GitHub, create an empty repository called `CDR_visualizer`. A public repository supports GitHub Pages on GitHub Free. Leave the initial README, license, and .gitignore unchecked because the local project already has files.
2. From this workspace root (the folder containing this README and the `CDR_visualizer` directory), run the commands below. Replace `YOUR_USERNAME` with your GitHub username. This folder already has a local Git repository. The existing remote is `https://github.com/SarthakMath/CDR_visualizer.git`; skip `git remote add origin` when using this checkout.

   ```bash
   git add .gitignore .github README.md CDR_visualizer
   git commit -m "Add CDR FEM FCT lab and GitHub Pages deployment"
   git branch -M main
   git remote add origin https://github.com/YOUR_USERNAME/CDR_visualizer.git
   git push -u origin main
   ```

3. In the GitHub repository, open **Settings → Pages**. Under **Build and deployment**, set **Source** to **GitHub Actions**.
4. Open **Actions → Deploy CDR FEM FCT lab to GitHub Pages → Run workflow** and run it on `main`. If the first push ran before Pages was enabled and failed, this manual run retries after configuration.
5. Wait for a successful deployment. The `github-pages` environment and deployment job show the actual site URL, normally `https://YOUR_USERNAME.github.io/CDR_visualizer/`.

Later pushes to `main` automatically update the site. If you upload files through the GitHub website instead, preserve this workspace's layout: `.github/workflows/deploy-pages.yml` belongs at the repository root, and the app remains in `CDR_visualizer/dist`. Uploading just the inner app folder will not include the workflow provided here.

GitHub Pages sites are ordinarily public websites; a private source repository does not by itself make the deployed site private. Private-repository Pages availability depends on the GitHub plan. See [GitHub's publishing-source documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site) and [custom workflow documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).


### Push a deployment fix

When updating the existing checkout, commit and push the changed deployment files:

```bash
git add .github/workflows/deploy-pages.yml .gitignore README.md
git commit -m "Fix GitHub Pages deployment folder"
git push origin main
```

This push starts a new workflow run automatically. The static upload path must match the checked-in app directory, currently `CDR_visualizer/dist`.

## Advection–diffusion–reaction

The app now includes a constant reaction coefficient c in steady and transient solvers, with pure reaction, diffusion–reaction, advection–reaction, and full ADR examples. Positive c causes decay; negative c causes growth. See [the reaction experiment guide](CDR_visualizer/README.md#reaction-experiments-decay-growth-and-transport) for decay, moving-layer, and stiff-reaction comparisons, including the precise mass-lumping convention.

## Layers and SUPG

Three layer experiments now compare Galerkin with full residual SUPG or isotropic artificial diffusion: a steady outflow layer, a moving interior front, and a developing transient boundary layer. See [the SUPG experiment guide](CDR_visualizer/README.md#layers-and-supg-three-reproducible-comparisons) for reproducible settings, results, and limitations.

## Algebraic flux correction

Try **Steady AFC layer** and **Nonnegative moving front** under **Flux correction**. Compare the low-order graph-diffusion method with limited flux correction. The steady method reports nonlinear iterations and convergence; transient FCT uses lumped mass and explicit SSPRK2 with automatic internal substeps. See [the AFC experiment guide](CDR_visualizer/README.md#algebraic-flux-correction-afc) for reproducible examples, costs, and limitations.

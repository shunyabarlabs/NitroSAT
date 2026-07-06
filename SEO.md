## Current SEO diagnosis

  NitroSAT barely surfaced in branded search during this audit. The GitHub repository
  currently has three stars, limited topic tags, and no published releases. Its search
  snippet leads with “linear-time MaxSAT approximator” and advanced manifold terminology
  rather than user problems. GitHub repository (https://github.com/sethuiyer/NitroSAT)

  The GitHub Pages site also has:

  - One large academic page
  - An extremely long, jargon-heavy title
  - No sitemap
  - No robots.txt
  - No structured data
  - No indexable tutorial pages
  - Claim-heavy descriptions instead of task-oriented copy

  Google’s own guidance emphasizes descriptive titles, logical site structure, internal
  links, unique useful content, and anticipating users’ search terms. Google SEO Starter
  Guide (https://developers.google.com/search/docs/fundamentals/seo-starter-guide)

  ## Highest-impact actions

  ### 1. Turn every tutorial into a real web page

  Publish clean URLs:

  /tutorials/sudoku-sat-solver/
  /tutorials/workforce-scheduling-maxsat/
  /tutorials/university-timetabling/
  /tutorials/vehicle-assignment/
  /tutorials/graph-coloring-sat/
  /tutorials/kubernetes-pod-placement/
  /tutorials/product-configuration/
  /tutorials/exam-scheduling/
  /tutorials/nurse-rostering/

  Markdown files inside GitHub are useful, but dedicated HTML pages have better titles,
  metadata, navigation, analytics, and internal linking.

  Each page should target one query:

   Tutorial          Primary keyword
  ━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   Sudoku            sudoku SAT solver / Sudoku CNF encoding
  ────────────────  ──────────────────────────────────────────
   Workforce         workforce scheduling MaxSAT
  ────────────────  ──────────────────────────────────────────
   Timetabling       university timetabling constraint solver
  ────────────────  ──────────────────────────────────────────
   Vehicle           vehicle assignment optimization
  ────────────────  ──────────────────────────────────────────
   Meetings          meeting room scheduling solver
  ────────────────  ──────────────────────────────────────────
   Graph coloring    graph coloring SAT encoding
  ────────────────  ──────────────────────────────────────────
   Kubernetes        Kubernetes pod placement optimization
  ────────────────  ──────────────────────────────────────────
   Configuration     product configuration constraint solver
  ────────────────  ──────────────────────────────────────────
   Exams             exam scheduling optimization
  ────────────────  ──────────────────────────────────────────
   Nurses            nurse rostering constraint programming

  ### 2. Rewrite the homepage title and first screen

  Current positioning is mathematically impressive but poor for discovery.

  Suggested title:

  > NitroSAT — Open-Source SAT and MaxSAT Solver for Scheduling and Large CNF Problems

  Suggested description:

  > Solve scheduling, configuration, graph coloring, Sudoku, and large CNF/WCNF problems with
  > an open-source C solver featuring streaming, disk-indexed search, and optional exact
  > solving.

  Put three immediate actions above the fold:

  - Install/build NitroSAT
  - Solve your first CNF
  - Browse executable tutorials

  Move manifold, topology, and number-theory explanations below the practical introduction.

  ### 3. Create technical SEO fundamentals

  Add:

  sitemap.xml
  robots.txt
  404.html

  Also add:

  - Unique <title> and meta description per tutorial
  - Canonical URLs
  - Open Graph images
  - SoftwareApplication JSON-LD on the homepage
  - TechArticle or HowTo semantics on tutorials
  - Breadcrumb navigation
  - Previous/next tutorial links
  - Google Search Console verification
  - Sitemap submission to Google and Bing

  ### 4. Improve GitHub conversion

  Update repository topics to include:

  sat-solver
  maxsat
  cnf
  wcnf
  constraint-programming
  combinatorial-optimization
  scheduling
  graph-coloring
  c
  optimization

  Publish an actual GitHub release with:

  - Linux binary
  - Checksums
  - Release notes
  - Quick-start examples
  - Exact/indexed-mode explanation

  “No releases published” makes adoption harder.

  ### 5. Replace broad claims with reproducible evidence

  “Linear-time SAT solver” and universal performance language will attract skepticism and
  increase bounce.

  Prefer:

  > V3 solved a planted 5,037,097-clause graph-coloring instance in 24.2 seconds using 84 MB
  > peak RSS; assignment independently verified.

  Create one page per benchmark family containing:

  - Generator and seed
  - Hardware
  - Commit hash
  - Exact command
  - Raw JSON
  - Peak RSS
  - Independent verifier
  - Visible failures
  - Comparison methodology

  Remove the “Note from GPT-5.5” from the public landing page. It does not establish external
  authority.

  ### 6. Build a content cluster

  Use the tutorial index as the hub. Every tutorial should link to:

  - CNF/WCNF introduction
  - V2 versus V3 guide
  - Exact versus heuristic solving
  - Relevant benchmark
  - Next related tutorial
  - GitHub repository

  Then publish foundational pages:

  /learn/what-is-cnf/
  /learn/what-is-wcnf/
  /learn/sat-vs-maxsat/
  /learn/how-to-encode-scheduling-as-sat/
  /learn/heuristic-vs-exact-sat-solving/
  /benchmarks/nitrosat-vs-established-solvers/

  ### 7. Earn relevant backlinks

  Distribute executable case studies—not promotional announcements—to:

  - DEV Community
  - Hacker News
  - Lobsters
  - SAT/constraint-programming communities
  - Operations-research communities
  - Kubernetes scheduling communities
  - University course/resource lists
  - “Awesome SAT solvers” and optimization lists
  - Papers with Code or comparable research indexes
  - Zenodo and Hugging Face dataset pages

  A good post title:

  > We encoded Kubernetes pod placement as weighted MaxSAT and made the example executable in
  > CI

  ## Recommended 30-day order

  1. Publish tutorial pages through GitHub Pages.
  2. Rewrite homepage title, description, and first screen.
  3. Add sitemap, robots, canonical tags, and structured data.
  4. Submit the sitemap to Search Console.
  5. Publish the first tagged GitHub release.
  6. Add repository topics.
  7. Publish one tutorial article weekly.
  8. Produce one controlled NitroSAT comparison benchmark.
  9. Track impressions, indexed pages, tutorial clicks, GitHub visits, clone events, and
     stars.

  The central strategy is:

  > Stop competing only for “SAT solver.” Own the long-tail searches where people are trying
  > to solve an actual scheduling, configuration, or allocation problem.


// Section "vibe coding ROWMAP" (5 slides) for a larger talk.
// Build (from the repository root):
//   typst compile --root . docs/slides/rowmap.typ

#let section-title = "vibe coding ROWMAP"
#let accent = rgb("#1f5f8b")
#let muted = luma(110)

#set page(paper: "presentation-16-9", margin: (x: 1.6cm, top: 2.3cm, bottom: 1.1cm),
  header: context [
    #set text(size: 11pt, fill: muted)
    #section-title #h(1fr) #counter(page).display() / #counter(page).final().first()
    #v(-6pt) #line(length: 100%, stroke: 0.5pt + muted)
  ])
#set text(font: "Latin Modern Sans", size: 17pt)
#show math.equation: set text(font: "New Computer Modern Math")
#show raw: set text(font: "DejaVu Sans Mono", size: 12pt)
#set list(spacing: 0.7em)
#set table(stroke: none)

#let slide(title, body) = {
  text(size: 26pt, weight: "bold", fill: accent, title)
  v(0.4em)
  body
  pagebreak(weak: true)
}

// ---------------------------------------------------------------- slide 1
#slide[From the 1996 report to Rust + Lean][
#grid(columns: (1.15fr, 1fr), gutter: 1.2cm,
[
- *Spec:* Weiner, Schmitt, Podhaisky, _ROWMAP_ (Report 39, Halle 1996).
  The 1997 Fortran code is a hint, not an oracle: 10 issues documented.
- *Rust:* 4300 lines, no dependencies; user supplies $J v$ (and $f_t$).
- *Methods:* GRK4T (order 4) and ROS34PW2 (order 3, L-stable).
  ROS3P rejected: its error estimate vanishes on linear problems.
- *Multiple Arnoldi process:* the minimal Krylov dimensions are
  _derived_, checked by order tests (deviation from exact ROW $= O(h^(p+1))$)
  and by Lean (subspace inclusions + tightness over $QQ$).
- *Finding:* the GRK4T stage-4 closure is not needed for classical order 4.
],
[
#set text(size: 16pt)
#align(center)[*Minimal dimensions $kappa_1 \/ kappa_2 \/ kappa_3 \/ kappa_4$*]
#v(0.3em)
#table(columns: 3, align: (left, center, center), inset: 7pt,
  table.hline(),
  [], [*autonomous*], [*non-aut.*],
  table.hline(stroke: 0.5pt),
  [GRK4T], [4/7/8/11], [4/9/10/15],
  [ROS34PW2], [3/4/4/4], [3/4/8/8],
  table.hline(),
)
#v(0.6em)
#set text(size: 14pt, fill: muted)
GRK4T: paper values reproduced.\
ROS34PW2: new; inserting $w_4$ is not needed (one prediction was wrong, the test caught it).
])
]

// ---------------------------------------------------------------- slide 2
#let prof = (1, 2, 3, 5, 6, 8, 9, 10, 10, 10, 10) // GRK4T aut., κ₁ = 4
#let staircase = {
  let n = prof.len()
  grid(columns: n * (13pt,), rows: n * (13pt,), gutter: 1.5pt,
    ..range(n).map(i => range(n).map(j => {
      let fill = if i <= j { accent.lighten(55%) }
        else if i <= prof.at(j) { rgb("#c0392b") } else { luma(240) }
      box(width: 13pt, height: 13pt, fill: fill)
    })).flatten())
}

#slide[Theorem: $H = Q^T A Q$ is generalized Hessenberg][
#grid(columns: (1.6fr, 1fr), gutter: 1cm,
[
MAP = sequence of operations: *insert* $w_i$, or *promote* the remainder of the
oldest open column (FIFO); each accepted or dropped. Column $j$ closes at
row $c(j)$; $o(m)$ = open columns when $q_m$ enters.

#v(0.3em)
*Theorem* (for orthogonal $q_i$ with the Gram–Schmidt flag property):
- *zero pattern:* $H_(i j) = 0$ for $i > c(j)$;
- *staircase:* closed columns form a prefix, $c$ is monotone;
- *rows:* $H_(m j) = 0$ for $j < m - o(m)$;
- *counts:* $o(m)$ depends only on the operations:\
  insert $+1$, promote $plus.minus 0$, drop $-1$.

#v(0.2em)
Every new vector has a _determined_ structural consequence; the inserted
vectors themselves never enter.
],
[
#align(center)[
#staircase
#v(0.3em)
#set text(size: 13pt, fill: muted)
GRK4T autonomous, $kappa_1 = 4$:\
#box(width: 9pt, height: 9pt, fill: rgb("#c0392b")) subdiagonal band, width $1 -> 2 -> 3$\
(insertion of $w_2$, $w_3$)
]
])
]

// ---------------------------------------------------------------- slide 3
#slide[Verified in Lean 4 + Mathlib][
#grid(columns: (1.25fr, 1fr), gutter: 1cm,
[
```lean
theorem H_structure (ops : List (Op × Bool))
    (A : E →ₗ[ℝ] E) (q : ℕ → E)
    (horth : ∀ i l, l < i → … → ⟪q i, q l⟫ = 0)
    (hprom : … origin m = some j →
       A (q j) ∈ span ℝ (q '' Iic m))
    (hdrop : … dropped j →
       A (q j) ∈ span ℝ (q '' Iic (c j))) :
  (zero pattern) ∧ (staircase) ∧ (rows) ∧
  (List.range k).map o = trace 0 ops
```
#v(0.2em)
#set text(size: 15pt)
`trace_*`, `bandwidths`: exact $o(m)$ sequences for all four schedules,
for *every* $kappa_1$.
],
[
#set text(size: 16pt)
- Lower bandwidth of $H$:\
  GRK4T *3* / *4*, ROS34PW2 *1* / *3*\
  (autonomous / non-autonomous)
- Axioms: `propext`, `Classical.choice`, `Quot.sound`; no `sorry`
- 430 lines; whole Lean project 1300 lines
- Rust test: observed $o(m)$ = Lean `trace`, exactly
- *Payoff:* incremental Givens QR, one per step; stage-1 residual
  $gamma h |l_kappa| norm(p_kappa)$ without back substitution
])
]

// ---------------------------------------------------------------- slide 4
#slide[Numerical results (paper problems, tol $10^(-2) dots 10^(-7)$)][
#grid(columns: (1fr, 1.05fr), gutter: 0.8cm,
[
#image("/docs/report/fig/wp_food_m50.pdf", width: 100%)
],
[
#set text(size: 14pt)
#table(columns: 4, align: (left, right, right, right), inset: 5pt,
  table.hline(),
  [], [*equal ERR*], [*equal ERR*], [*per $J v$*],
  [], [GRK4T], [ROS34PW2], [1-pass],
  table.hline(stroke: 0.5pt),
  [HEAT], [1.19], [*0.57*], [*0.79*],
  [NILIDI], [1.52], [1.93], [1.09],
  [BRUSS], [*0.92*], [1.01], [*0.75*],
  [FOOD], [1.36], [1.16], [*0.78*],
  [VDPL], [1.56], [1.94], [*0.90*],
  [DIFKON], [1.35], [2.46], [*0.96*],
  [DIFFU2], [1.21], [*0.90*], [*0.95*],
  table.hline(),
)
#set text(size: 13pt, fill: muted)
CPU time Rust / Fortran (< 1: Rust faster).
Same step counts, equal or better ERR.
])
#v(-0.2em)
#set text(size: 15pt)
Same algorithm (1 Gram–Schmidt pass): Rust faster per $J v$ on 6 of 7.
Two passes cost 25–56 % but are needed: one pass fails 13 tests (garbage at breakdown).
]

// ---------------------------------------------------------------- slide 5
#slide[The coder: Claude Opus 5.5 in Claude Code][
#grid(columns: (1fr, 1.1fr), gutter: 1cm,
[
#set text(size: 14pt)
#table(columns: 2, align: (left, right), inset: (x: 6pt, y: 4.5pt),
  table.hline(),
  [Model], [`claude-opus-5-5`],
  [Session], [31 Oct 2026],
  [Wall-clock time], [7.5 h],
  [Active time#super[1]], [4.2 h],
  [Model calls], [551],
  [Output tokens#super[2]], [0.78 M],
  [Input: cache reads], [261 M],
  [Input: cache writes], [2.0 M],
  [Input: uncached], [1 200],
  table.hline(stroke: 0.5pt),
  [*Cost at API prices#super[3]*], [*≈ \$84*],
  table.hline(),
)
#v(-0.3em)
#set text(size: 10pt, fill: muted)
#super[1] gaps < 10 min; includes tool runs (Lean, Fortran, benchmarks).
#super[2] incl. thinking.
#super[3] \$4 / \$20 per M input/output, cache read \$0.20, 1-h cache write \$8;
counted from the session transcript.
],
[
#set text(size: 15pt)
- *Human role:* interview decisions, review gates for derivations,
  and corrections (about 30 prompts).
- *The decisive correction came from the human:* "the MAP has a
  generalized Hessenberg structure … you missed it entirely" → theorem,
  Lean proof, incremental Givens QR.
- *Agent mistakes caught by tests:* a wrong `drop` optimization broke
  $H = Q^T A Q$; weakened re-orthogonalization lost orthogonality; 6 of 10
  omission predictions were wrong.
- *Output:* 4300 lines of Rust, 1300 lines of Lean, 48 tests, 8 commits,
  report + plots.
])
]
